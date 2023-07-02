# simplied version 
# x ~ N(1,5)
# if x >c , observe the uncensored dataset
# if x < c , censored dataset 

set.seed(8888)
mu <- 3
sigma <- 5 
N <- 300
x <- rnorm(N,mean = mu, sd = sigma)
c <- 1
x_ob <- x[x>c]
theta0 <- c(mu,sigma)
neglik_norm <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  lik_unc <- sum(dnorm(x_ob, mean = mu, sd = sigma, log = TRUE))
  lik_cen <- (N - length(x_ob))*log(pnorm(c,mean = mu,sd = sigma))
  return(-lik_cen-lik_unc)
}

optim(theta0,neglik_norm)

optimCheck::optim_proj(theta0,
                       neglik_norm,
                       xrng = 1,
                       xnames = c("mu","sigma"))

# sigma_seq <- seq(from = 0.1, to = 10, length = 200)
# neglik_simple_cen <- function(sigma){
#   # lik_unc <- sum(dnorm(x_ob, mean = mu, sd = sigma, log = TRUE))
#   lik_cen <- (N - length(x_ob))*log(pnorm(c,mean = mu,sd = sigma))
#   return(-lik_cen)
# }
# 
# neglik_simple_unc <- function(sigma){
#   lik_unc <- sum(dnorm(x_ob, mean = mu, sd = sigma, log = TRUE))
#   #lik_cen <- (N - length(x_ob))*log(pnorm(c,mean = mu,sd = sigma))
#   return(-lik_unc)
# }
# 
# plot(sigma_seq,sapply(sigma_seq,neglik_simple_cen),type = "l")
# plot(sigma_seq,sapply(sigma_seq,neglik_simple_unc),type = "l")


# Now consider a left-censoring problem in exponential distribution
# x ~ exp(5),
# if x >c , observe the uncensored dataset
# if x < c , censored dataset


lambda <- 3
x <- rexp(N,rate = lambda)
x_ob <- x[x>c]

neglik_exp <- function(lambda){
  lik_unc <- sum(dexp(x_ob, rate = lambda, log = TRUE))
  lik_cen <- (N - length(x_ob))*log(pexp(c,rate = lambda))
  return(-lik_cen-lik_unc)
}
optimize(neglik_exp,c(0,6))

optimCheck::optim_proj(lambda,
                       neglik_exp,
                       xrng = 1,
                       xnames = c("lambda"))



# Now consider a left-censoring problem in weilbull distribution
shape <- 2
scale <- 5
x <- rweibull(n = N,shape = shape, scale = scale)
x_ob <- x[x>c]
theta0<- c(shape,scale)
neglik_weibull <- function(theta){
  shape <- theta[1]
  scale <- theta[2]
  lik_unc <- sum(dweibull(x_ob, shape = shape, scale = scale, log = TRUE))
  lik_cen <- (N - length(x_ob))*log(pweibull(c,shape = shape, scale = scale))
  return(-lik_cen-lik_unc)
}



optim(theta0,neglik_weibull)

optimCheck::optim_proj(theta0,
                       neglik_weibull,
                       xrng = .75,
                       xnames = c("Shape","Scale"))
