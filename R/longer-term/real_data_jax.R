#setwd("~/Downloads/github/Bayesian-lumber-strength/R/longer-term")
bending <- read.csv("bending-pl.csv", header = T)
# convert psi to Mpa
#  1 thousand psi = 6.895 MPa
bending[,1] <- bending[,1]/1000*6.895
l <- 4500/1000*6.895
# bending_mar <- bending[1:195,]
# bending_dmg <- bending[195:341,]

R100_data <- bending[bending[,2] == "R100",1]
R20_data <- bending[bending[,2] == "R20",1]
R20R100_data <- bending[bending[,2] == "R20R100",1]

g2 <- 97


logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}
negdmglik_theta <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  #alpha <- expit(theta[3])
  c <- theta[4]
  
  
  lik1 <- sum(dnorm(R20_data,mean = mu, sd = sigma, log = TRUE))
  lik2 <- g2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- dmglik_jax(R20R100_data,alpha,l,c,s,mu,sigma)
  lik4 <- sum(dnorm(R100_data,mean = mu, sd = sigma, log = TRUE))
  return(-lik1-lik2-lik3-lik4)
}
negdmglik_theta(theta0)

optimout <- optim(theta0,negdmglik_theta,method = "L-BFGS-B",
                  lower = rep(0.1,4),upper = c(Inf,Inf,4,Inf))


# optimout <- optim(theta0,negdmglik_theta,method = "L-BFGS-B",
#                   lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))
# par_optim = optimout$par
# par_optim[3] = logit(optim)