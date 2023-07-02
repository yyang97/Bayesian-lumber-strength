# the original model 
# y* = y* I(c*x > l) + alpha * y * I(c*x < l)
# first maybe simplified as 
# y* = y* I(x > l) + alpha * y * I(x < l)

#
# yseq <- seq(from = 1, to = 100, by =.1)
# ystarseq <- ((1/0.2)*yseq + 32/.8 - 1/0.2*32/.8)*ifelse(.8*yseq<32 & yseq > 32,1,0) + yseq*ifelse(.8*yseq>=32 | yseq < 32,1,0)
# 
# plot(yseq,ystarseq, type = 'l', ylim =c(0,100))
# # the dash line is no damaged plot 
# abline(yseq, yseq, lty = 2)


x <- seq(from = -10, to= 10, length.out = 500)
ind_y <- ifelse(x>0,1,0)
plot(x,ind_y,type = "l")
lines(x,sapply(x,smooth_ind, s = 1),type = "l",col = "red")
lines(x,sapply(x,smooth_ind, s = 10),type = "l",col = "blue")
lines(x,sapply(x,smooth_ind, s = 100),type = "l",col = "green")
legend("topleft", legend=c("Indicator", "s = 1", "s = 10", "s = 100"),
       col=c("black","red", "blue","green"), lty=1, cex=0.8)


### new 

#' The smooth function for the indicator function, which is also called "sigmoid function".
#' 
#' @param x, the variable 
#' @param s, the hyper parameter, higher s means closer the indicator function
#' @returns the smoothed value.  

smooth_ind <- function(x,s){
  return(1/(1+exp(-x*s)))
}



#' The damage model, smoothed version of 
#' "y* = y* I(c*x > l) + alpha * y * I(c*x < l)"
#' 
#' @param y, the original strength
#' @param alpha, the damage parameter
#' @param l, the proof loading level
#' @param c, the threshold parameter. Damage effects happen exceeding c.
#' @param s, the temperature parameter for s
#' @return, the weakened y*, y*<y.
dmg_model <- function(y,alpha,l,c,s){
  return(y*smooth_ind(c*y-l,s)+ alpha*y*smooth_ind(l-c*y,s))
}

#' The damage-inverse model, given a weakened ystar, find the original y
#' 
#' @param ystar, the weakened strength
#' @param alpha, the damage parameter
#' @param l, the proof loading level
#' @param c, the threshold parameter. Damage effects happen exceeding c.
#' @param s, the temperature parameter for s
#' @return, the original y, y>y*.
dmg_inverse <- function(ystar,alpha,l,c,s){
  uniroot((function (x) dmg_model(x,alpha,l,c,s) - ystar), lower = 0, upper = 1000)$root
}
# the abs gradient at damage-inverse
dmg_inverse_grad <- function(ystar,alpha,l,c,s){
  abs(numDeriv::grad(func = (function(x) dmg_inverse(x,alpha,l,c,s)),
                     x  = ystar))
}  
# dmg_inverse_grad(42,a lpha,l,c,s)
# ystar <- y*smooth_ind(c*y-l)+ alpha*smooth_ind(l-c*y)

dmglik <- function(ystar,alpha,l,c,s,mu,sigma){
  y <- dmg_inverse(ystar,alpha,l,c,s)
  dnorm(y, mean = mu, sd = sigma, log = TRUE)+ log(dmg_inverse_grad(ystar,alpha,l,c,s))
}

dmglik(30,alpha = 0.4,l,c,s,mu,sigma)


# y <- seq(from  = 0.01, to = 100, length = 200)
# ystar <- dmg_model(y,alpha,l,c,s)
# plot(y,ystar,type = 'l',lty = 1)
# lines(y,y,lty = 2)
# legend("topleft", legend = c("ystar","y"),lty=1:2)
# abline(v = l)
# abline(v = l/c)

# 
# 
# ystar <- seq(from  = 0.1, to = 100, length = 200)
# y <- sapply(ystar,dmg_inverse,alpha = alpha,l = l,c = c ,s = s)
# plot(ystar,y,type = 'l',lty = 1)
# lines(ystar,ystar,lty = 2)
# legend("topleft", legend = c("y","ystar"),lty=1:2)
# abline(v = l)
# abline(v = l/c)







## plot 
set.seed(7777)
N <-300
mu <- 48
sigma <- 19
l <-  32
alpha <- 0.45
c <- 0.65
s <- 1

## data 

y <- rnorm(N, mean = mu, sd = sigma )
y <- y[y>0]
y_obs <- list()
# group 1, y^* <y < l
y_obs$group1 <- y[which(y < l)]
# group 2, y^*<l<y
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23,alpha,l,c,s)
y_obs$group2 <- length(group23[which(group23 > l &group23_star < l)])
y_obs$group3 <- group23_star[which(group23_star >l)]

negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik2-lik3)
  #return (-lik3)
}

theta0 <- c(mu,sigma,alpha,c)


optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
      lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimCheck::optim_proj(optimout$par,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))



# the only alpha plot
s <- 1
negdmglik_alpha  <- function(alpha){
  mu <- mu
  sigma <- sigma
  c <- c
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  # return(-lik1-lik2-lik3)
  return (-lik2-lik3)
}

alpha_seq <- seq(from = 0.01, to = 0.99, by = 0.01)
negloglik <- sapply(alpha_seq,negdmglik_alpha)

plot(alpha_seq,negloglik, type = "l")

plot(alpha_seq[1:30],negloglik[1:30],type = "l")


# when alpha is small, why it has little impact on the likelihood?


## group 2 likelihood 
negdmglik_alpha_g2  <- function(alpha){
  mu <- mu
  sigma <- sigma
  c <- c
    lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )

  return (-lik2)
}
likg2 <- sapply(alpha_seq,negdmglik_alpha_g2)

## group 3 likelihood 
negdmglik_alpha_g3  <- function(alpha){
  mu <- mu
  sigma <- sigma
  c <- c
  

  #lik2 <- 0
  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  # return(-lik1-lik2-lik3)
  return (-lik3)
}
likg3 <- sapply(alpha_seq,negdmglik_alpha_g3)


plot(alpha_seq,likg2, type = "l")
plot(alpha_seq,likg3, type = "l")

# 
# negdmglik_g1g3 <- function(theta){
#   mu <- theta[1]
#   sigma <- theta[2]
#   alpha <- theta[3]
#   c <- theta[4]
# 
#   lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
# 
#   lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
#   return(-lik1-lik3)
# }
# 
# 
# optimCheck::optim_proj(theta0,
#                        negdmglik_g1g3,xrng = .5,
#                        xnames = c("mu","sigma","alpha","c"))
