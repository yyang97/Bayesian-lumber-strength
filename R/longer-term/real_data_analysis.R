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



negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
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

theta_test <- c(mu,sigma,0.5,.6)
negdmglik_model(theta_test)



s <- 1
optimout <- optim(par = theta_test,
                  fn = negdmglik_model,
                  method = "L-BFGS-B",
                  lower = c(-Inf,-Inf,.1,-Inf),
                  upper = c(Inf,Inf,.9,Inf))
optimout$par


# [1] 1678.031
#theta_test <- c(47.9985218, 18.9804222,0.4548490,  0.6515348)
optimCheck::optim_proj(optimout$par,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))


