setwd("~/Downloads/github/Bayesian-lumber-strength/R/longer-term")
library(rstan)
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




##--------------stan fitting----------------
standata <- list(N_R20 = length(R20_data),N_R100 = length(R100_data),N_R20R100 = length(R20R100_data),
                 X_R20 = R20_data, X_R100 = R100_data, X_R20R100 = R20R100_data,
                 l = l, thresh = 0.91)
init_dmg <- function() {
  list(mu = 45, sigma = 20, alpha =0.5)
}
options(mc.cores = parallel::detectCores()) 
# dmg_mod <- stan_model("longerterm.stan")
dmgfit <- rstan::sampling(dmg_mod, data = standata, init = init_dmg)

print(dmgfit)
traceplot(dmgfit)



standata <- list(N_R20 = length(R20_data),N_R100 = length(R100_data),N_R20R100 = length(R20R100_data),
                 X_R20 = R20_data, X_R100 = R100_data, X_R20R100 = R20R100_data,
                 l = l, thresh = 0.85)
init_dmg <- function() {
  list(mu = 45, sigma = 20, alpha =0.97)
}
options(mc.cores = parallel::detectCores()) 
dmg_mod <- stan_model("threegroups.stan")
dmgfit <- rstan::sampling(dmg_mod, data = standata, init = init_dmg)

print(dmgfit)
traceplot(dmgfit)
# 0.90 -1798.90
# 


summary(dmgfit)[[1]][4,1]
thresh_seq <- seq(from = .7, to = 0.99,by = .01)





logpost <- c()
for (i in 1:length(thresh_seq)){
  standata <- list(N_R20 = length(R20_data),N_R100 = length(R100_data),N_R20R100 = length(R20R100_data),
                   X_R20 = R20_data, X_R100 = R100_data, X_R20R100 = R20R100_data,
                   l = l, thresh = thresh_seq[i])
  init_dmg <- function() {
    list(mu = 45, sigma = 20, alpha =0.999)
    
  }
  dmgfit <- rstan::sampling(dmg_mod, data = standata, init = init_dmg)
  logpost[i] <- summary(dmgfit)[[1]][4,1]
  
  print(i)
}

saveRDS(logpost, file = 'logpost.rds')

logpost <- readRDS('logpost.rds')

plot(thresh_seq[c(10:30)],logpost[c(10:30)],type = 'l')


# thresh = 0.91,  -1505.40
# thresh = 0.81,  -1504.97 
cbind(thresh_seq, logpost)
which.max(logpost)
thresh_seq[which.max(logpost)]


thresh <-  .82
dmglik <- function(theta){
  # mu <- theta[1]
  # sigma <- theta[2]
  mu <- 47.91
  sigma <- 18.84
  alpha <- theta
  lik <- 0
  for (jj in 1:length(R20R100_data)){
    if(R20R100_data[jj]/alpha < l/thresh){
      lik <- lik - log(alpha) + dnorm(R20R100_data[jj]/alpha, mean =  mu,sd = sigma, log = TRUE)
    }else if(R20R100_data[jj]/alpha > l/thresh && R20R100_data[jj] > l/thresh){
      lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)
    }
    else if(R20R100_data[jj]/alpha > l/thresh && R20R100_data[jj] < l/thresh){
      lik <- lik-10000
    }
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  return(-1*lik)
}



alphaseq <- seq(from = 0.1,to = 1.2, by = 0.01)
likseq <- c()
for (iseq in 1:length(alphaseq) ){
  likseq[iseq] <- dmglik(alphaseq[iseq])
  
}

plot(alphaseq,likseq,type = 'l')







# thresh <- 0.81
thresh <-  0.81
cdfdmg <- function(alpha, ystar){


  mu <- 49.67
  sigma <- 20.20

 cdf <-  pnorm(min(ystar/alpha,l/thresh ),mean = mu, sd = sigma) + max(0,pnorm (ystar ,mean = mu, sd = sigma) - pnorm(l/thresh,mean = mu, sd = sigma))
  return(cdf)
}


xseq <- seq(0,100, by = .1)
yseq <- c()
for (jj in 1:length(xseq)){
yseq[jj] <- cdfdmg(1.03, xseq[jj])
}

plot(xseq, yseq , type = 'l',main = "alpha = 1.03, thresh = 0.81")



## thresh = 0.91
thresh <-  0.91
cdfdmg <- function(alpha, ystar){
  
  
  mu <- 49.71
  sigma <- 20.14
  
  cdf <-  pnorm(min(ystar/alpha,l/thresh ),mean = mu, sd = sigma) + max(0,pnorm (ystar ,mean = mu, sd = sigma) - pnorm(l/thresh,mean = mu, sd = sigma))
  return(cdf)
}


xseq <- seq(0,100, by = .1)
yseq <- c()
for (jj in 1:length(xseq)){
  yseq[jj] <- cdfdmg(.56, xseq[jj])
}

plot(xseq, yseq , type = 'l',main = "alpha = 0.56, thresh = 0.91")


## thresh = 0.9, lp__  -1505.41
## thresh = 0.8, lp__ -1504.86
## thresh = .7, lp__ -1505.36
## thresh = 0.82, lp__ -1505.05


# hist(bending_mar[,1])
# shapiro.test(bending_mar[1:138,1])
# shapiro.test(sqrt(bending_mar[1:138,1]))
# hist(log(bending_mar[1:138,1]))
# 
# 
# FITnorm <- fitdistrplus::fitdist(R100_data, "norm")
# 
# shapiro.test(R100_data)
# 
# plot(FIT)
FITnorm <- fitdistrplus::fitdist(R100_data, "norm")
FITweibull <- fitdistrplus::fitdist(R100_data, "weibull")

FITnorm$bic
FITweibull$bic
FITweibull$estimate








# mu = 47.91, sd = 18.84
R100_data <- rweibull(138, shape = 2.76, scale = 53.997)
l <- qweibull(.2, shape = 2.76, scale = 53.997)


thresh <- 0.9
alpha <- .5

samples <- rweibull(500, shape = 2.76, scale = 53.997)

R20_data <- samples[which(samples < l)]
R20R100_data <- samples[which(samples > l)]
for (ii in 1: length(R20R100_data)){
  if(l > thresh*R20R100_data[ii]){
    R20R100_data[ii] <- alpha * R20R100_data[ii]
  }
}


length(R20R100_data)
sum(l > thresh*R20R100_data)


##--------------stan fitting----------------
standata <- list(N_R20 = length(R20_data),N_R100 = length(R100_data),N_R20R100 = length(R20R100_data),
                 X_R20 = R20_data, X_R100 = R100_data, X_R20R100 = R20R100_data,
                 l = l, thresh = .9)
init_dmg <- function() {
  list(mu = 45, sigma = 20, alpha =0.5)
}
options(mc.cores = parallel::detectCores()) 
# dmg_mod <- stan_model("longerterm.stan")
dmgfit <- rstan::sampling(dmg_mod, data = standata, init = init_dmg,iter = 4000)

print(dmgfit)
traceplot(dmgfit)

## simulation





###-------------new model ----------------
setwd("~/Downloads/github/Bayesian-lumber-strength/R/longer-term")
library(rstan)
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

length(R20R100_data)
length(R20_data)
length(R100_data)



# R20R100_observed <- dataset[[3]][,1]
dmglik <- function(theta,thresh){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  est_y <-c()
  cond <- c()
  # mu <- 47.91
  # sigma <- 18.84
  # alpha <- theta
  # R100 group
  lik <- sum(dnorm(R100_data, mean = mu, sd = sigma, log = TRUE))
  # group 1
  lik <- lik + sum(dnorm(R20_data, mean = mu, sd = sigma, log = TRUE) - log(pnorm(l, mean  = mu, sd = sigma)))
  # group 2
  lik <- lik + 97*log(pnorm(l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma ))
  cond <-c()
  likindi <- c()
  for (jj in 1:length(R20R100_data)){
    if(R20R100_data[jj]  < l/thresh){
      # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
      lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      est_y[jj] <- alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh
      cond[jj] <- 1
      }
    else
    {
      lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)- log(1 - pnorm(l/thresh, mean = mu, sd = sigma))
      est_y[jj] <- R20R100_data[jj]
      cond[jj] <- 0
      }
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  #liktable <- data.frame(R20R100_data,est_y,cond)
  return(-1*lik)
  #return(liktable)
}
dmglik(theta = c(47.91,18.84,.8), thresh = .5)
optimout <- optim(c(47.91,18.84,.8),thresh = .9,
                  dmglik, 
                  lower = c(20,3,.001),
                  upper = c(80,30,10), 
                  method = "L-BFGS-B")
optimout




threshseq <- seq(from = .1, to = .99, by = .01)
neglik <- c()
for (jj in 1:length(threshseq) ){
  optimout <- optim(c(47.91,18.84,.8),thresh = threshseq[jj],
                    dmglik, 
                    lower = c(20,10,.001),
                    upper = c(80,30,10), 
                    method = "L-BFGS-B")
  neglik[jj] <- optimout$value
}

plot(threshseq, neglik, type = 'l')

threshseq[which.min(neglik)]



optimout <- optim(c(47.91,18.84,.8),thresh = .39,
                  dmglik, 
                  lower = c(20,10,.001),
                  upper = c(80,30,10), 
                  method = "L-BFGS-B")
optimout$par


dmglik <- function(theta,thresh){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  est_y <-c()
  cond <- c()
  # mu <- 47.91
  # sigma <- 18.84
  # alpha <- theta
  # R100 group
  lik <- sum(dnorm(R100_data, mean = mu, sd = sigma, log = TRUE))
  # group 1
  lik <- lik + sum(dnorm(R20_data, mean = mu, sd = sigma, log = TRUE) - log(pnorm(l, mean  = mu, sd = sigma)))
  # group 2
  lik <- lik + dataset[[2]]*log(pnorm(l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma ))
  cond <-c()
  likindi <- c()
  for (jj in 1:length(R20R100_data)){
    if(R20R100_data[jj]  < l/thresh){
      # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
      lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      est_y[jj] <- alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh
      cond[jj] <- 1
    }
    else
    {
      lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)- log(1 - pnorm(l/thresh, mean = mu, sd = sigma))
      est_y[jj] <- R20R100_data[jj]
      cond[jj] <- 0
    }
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  liktable <- data.frame(R20R100_data,est_y,cond)
  #return(-1*lik)
  return(liktable)
}

dmglik(c(52.819832 ,2.600892,1.192509), 
       thresh = .39)
