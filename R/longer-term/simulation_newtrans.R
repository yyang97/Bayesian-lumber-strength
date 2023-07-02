
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






## plot 
set.seed(8888)
# the orignal sample size 
N <- 300
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
# alpha <- 0.9
c <- 0.8
alpha <- 0.9
# alpha <- 0.24
# c <- 0.65
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
optimCheck::optim_proj(theta0,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))


optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
                  lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimCheck::optim_proj(optimout$par,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))