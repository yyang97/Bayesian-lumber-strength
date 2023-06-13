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



# y <- seq(from  = 0.01, to = 100, length = 200)
# ystar <- dmg_model(y,alpha,l,c,s)
# plot(y,ystar,type = 'l',lty = 1)
# lines(y,y,lty = 2)
# legend("topleft", legend = c("ystar","y"),lty=1:2)
# abline(v = l)
# abline(v = l/c)



ystar <- seq(from  = 0.1, to = 100, length = 200)
y <- sapply(ystar,dmg_inverse,alpha = alpha,l = l,c = c ,s = s)
plot(ystar,y,type = 'l',lty = 1)
lines(ystar,ystar,lty = 2)
legend("topleft", legend = c("y","ystar"),lty=1:2)
abline(v = l)
abline(v = l/c)


dmglik(2,alpha = alpha,l = l,c = c,s = 1,mu = mu,sigma = sigma)
alpha

dmglik(42,alpha = alpha,l = l,c = c,s = 1,mu = mu,sigma = sigma)
dmglik(42,alpha = .4,l = l,c = c,s = 1,mu = mu,sigma = sigma)





## plot 
set.seed(8888)
N <- 1000
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
alpha <- 0.9
c <- 0.8
s <- 1

## data 

y <- rnorm(N, mean = mu, sd = sigma )

y_obs <- list()
# group 1, y^* <y < l
y_obs$group1 <- y[which(y < l)]
# group 2, y^*<l<y
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23,alpha,l,c,s)
y_obs$group2 <- length(group23[which(group23 > l &group23_star < l)])
y_obs$group3 <- group23[which(group23_star >l)]

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
}

theta0 <- c(mu,sigma,alpha,c)
#theta0 <- c(48,19,.7,.8)
negdmglik_model(theta0)
optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
      lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimCheck::optim_proj(theta0,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))


negdmglik_g1g3 <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))

  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik3)
}


optimCheck::optim_proj(theta0,
                       negdmglik_g1g3,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))
