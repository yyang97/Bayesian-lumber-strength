bending <- read.csv("bending-pl.csv", header = T)

bending[,1] <- bending[,1]/1000*6.895
l <- 4500/1000*6.895
# bending_mar <- bending[1:195,]
# bending_dmg <- bending[195:341,]

R100_data <- bending[bending[,2] == "R100",1]
R20_data <- bending[bending[,2] == "R20",1]
R20R100_data <- bending[bending[,2] == "R20R100",1]


# 97 pieces breaking in the loading process.
g2 <- 97

# the total sample size 
#N <- length(R20_data) + length(R20R100_data) + g2
mu <- 48 
sigma <- 19

s <- 1


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

logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}

negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  # alpha <- theta[3]
  alpha <- expit(theta[3])
  c <- theta[4]
  
  lik1 <- sum(dnorm(R20_data,mean = mu, sd = sigma, log = TRUE))
  lik2 <- g2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(R20R100_data,dmglik,alpha,l,c,s,mu,sigma))
  lik4 <- sum(dnorm(R100_data,mean = mu, sd = sigma, log = TRUE))
  return(-lik1-lik2-lik3-lik4)
}

negdmglik_model_grad <- function(theta){
  numDeriv::grad(negdmglik_model,theta_test)
}
theta_test <- c(mu,sigma,logit(0.5),.6)
negdmglik_model(theta_test)


theta_test <- c(mu,sigma,logit(0.5),.6)
theta_0 <- c(47.9985218, 18.9804222,logit(0.4548490),  0.6515348)
negdmglik_model(theta_0)
optimout <- optim(par = theta_0,
                  fn = negdmglik_model,
                  gr = negdmglik_model_grad,
                  method = "BFGS")
optimout$par

expit(optimout$par[3])

# [1] 1678.031
#theta_test <- c(47.9985218, 18.9804222,0.4548490,  0.6515348)





s <- 1
negdmglik_model_unconstraint <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  # alpha <- theta[3]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(R20_data,mean = mu, sd = sigma, log = TRUE))
  lik2 <- g2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(R20R100_data,dmglik,alpha,l,c,s,mu,sigma))
  lik4 <- sum(dnorm(R100_data,mean = mu, sd = sigma, log = TRUE))
  return(-lik1-lik2-lik3-lik4)
}

theta_opt <- optimout$par
theta_opt[3] <- expit(theta_opt[3])
optimCheck::optim_proj(theta_opt,
                       negdmglik_model_unconstraint ,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))
# 
## the only plot for alpha

negdmglik_alpha <- function(alpha){
  mu <- 47.9985218
  sigma <- 18.9804222
  # alpha <- theta[3]
  c <- 0.6515348

  lik1 <- sum(dnorm(R20_data,mean = mu, sd = sigma, log = TRUE))
  lik2 <- g2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(R20R100_data,dmglik,alpha,l,c,s,mu,sigma))
  lik4 <- sum(dnorm(R100_data,mean = mu, sd = sigma, log = TRUE))
  return(-lik1-lik2-lik3-lik4)
}

s <- 1
alpha_opt <- optimize(negdmglik_alpha,c(.01,.99))


alpha_seq <- seq(from = 0.01, to = 0.5, length =200)
plot(alpha_seq, sapply(alpha_seq,negdmglik_alpha),type = "l", ylab = "neglik")
# 
# 
# 
# 
# negdmglik_alpha <- function(alpha){
#   mu <- 47.9985218
#   sigma <- 18.9804222
#   alpha <- expit(alpha)
#   c <- 0.6515348
#   
#   lik1 <- sum(dnorm(R20_data,mean = mu, sd = sigma, log = TRUE))
#   lik2 <- g2 * log(
#     pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
#       pnorm(l, mean = mu, sd = sigma)
#   )
#   #lik2 <- 0
#   lik3 <- sum(sapply(R20R100_data,dmglik,alpha,l,c,s,mu,sigma))
#   lik4 <- sum(dnorm(R100_data,mean = mu, sd = sigma, log = TRUE))
#   return(-lik1-lik2-lik3-lik4)
# }
# 
# 
# alpha <- seq(from = -1, to = 1, length = 200)
# plot(alpha, sapply(alpha,negdmglik_alpha),type = "l", ylab = "neglik")



##-----------test  the result--------



test_data <- cbind(R20R100_data,
                   alpha45 = sapply(R20R100_data, dmg_inverse,0.45,l,0.65,s),
                   alpha25 = sapply(R20R100_data, dmg_inverse,0.25,l,0.65,s))
test_data 
test_data[,1] < test_data[,2]
