setwd("~/Downloads/github/Bayesian-lumber-strength/R/longer-term")
require(rstan)
bending <- read.csv("bending-pl.csv", header = T)
# convert psi to Mpa
#  1 thousand psi = 6.895 MPa
bending[,1] <- bending[,1]/1000*6.895

R100_realdata <- bending[bending[,2] == "R100",1]
R20_realdata <- bending[bending[,2] == "R20",1]
R20R100_realdata <- bending[bending[,2] == "R20R100",1]

N_R100 <- length(R100_realdata)
N_R20 <- length(R20_realdata)
N_R20R100 <- length(R20R100_realdata)


# mu = 47.91, sd = 18.84
R100_data <- rnorm(138, mean = 47.91, sd =  18.84)
l <- qnorm(.2, mean =  47.91, sd = 18.84)

thresh <- 0.6
alpha <- .8

samples <- rnorm(500, mean = 47.91, sd =  18.84)

R20_data <- samples[which(samples < l)]
R20R100_data <- samples[which(samples > l)]
for (ii in 1: length(R20R100_data)){
  if(l > thresh*R20R100_data[ii]){
    R20R100_data[ii] <- alpha * R20R100_data[ii] + l/thresh - alpha*l/thresh
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
dmg_mod <- stan_model("longerterm.stan")
dmgfit <- rstan::sampling(dmg_mod, data = standata, init = init_dmg)
print(dmgfit)
traceplot(dmgfit)

# alpha = 0.5
# thresh = 0.75, -2737.99  
# thresh = 0.8,  -2749.09
# thresh= 0.9，  -2771.01 




# true thresh .9 lp__ -2745.12  
# thresh = 0.8, lp__ -2758.72 
# thresh = 0.82, lp__ -2764.83
# thresh = 0.84, lp__  -2767.06
# thresh = 0.86, lp__  -2767.43
# thresh = 0.88, lp__  -2767.55  

# 



thresh <-  .6
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




optimdmg <- optimize(dmglik, c(0,20))
optimdmg$minimum
optimCheck::optim_proj(.75,dmglik, xrng = 1)


dmglik(.5)
dmglik(.7)

dmgvalue <- c()
xseq <- seq(from = 0, to = 1, by= 0.1)
for (jj in 1:length(xseq)){
  dmgvalue[jj] <- dmglik(xseq[jj])
}

plot(xseq, dmgvalue,type = 'l')





l <- qnorm(.2, mean =  47.91, sd = 18.84)

thresh <- 0.5
alpha <- .8

samples <- rnorm(100000, mean = 47.91, sd =  18.84)

R20_data <- samples[which(samples < l)]
R20R100_data <- samples[which(samples > l)]

R20R100_observed <- c()
damage_table <- c()
for (ii in 1: length(R20R100_data)){
  if(l > thresh*R20R100_data[ii]){
    R20R100_observed[ii] <- 1/alpha * R20R100_data[ii] + l/thresh - 1/alpha*l/thresh
    # 1 means 'damaged'
    damage_table[ii] <- 1
  }
  else{
    # 0 means 'nondamaged'
    R20R100_observed[ii] <- R20R100_data[ii]
    damage_table[ii] <- 0
  }
}
damage_data <- cbind(R20R100_observed,damage_table,R20R100_data)

sum(damage_data[,2] == 0)
sum(damage_data[,2] == 1)

damage_data <- damage_data[damage_data[,1] > l,]
R20R100_observed<- damage_data[,1]





# cbind(dmglik(.8), damage_data)
# sum(dmglik(.6)$likindi)
# sum(dmglik(.8)$likindi)

# R20R100_observed <- R20R100_observed[which(R20R100_observed>l)]
thresh <- 0.5
dmglik <- function(theta){
  # mu <- theta[1]
  # sigma <- theta[2]
  mu <- 47.91
  sigma <- 18.84
  alpha <- theta
  data <- c()
  est_y <-c()
  lik <- 0
  cond <-c()
  likindi <- c()
  for (jj in 1:length(R20R100_observed)){
    data[jj] <- R20R100_observed[jj]
    if(R20R100_observed[jj]  < l/thresh){
      # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
      lik <- lik + log(alpha) + dnorm(alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      likindi[jj] <- log(alpha) + dnorm(alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE)- log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      est_y[jj] <- alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh
      cond[jj] <- 1
      }
    else{
      lik <- lik + dnorm(R20R100_observed[jj], mean = mu, sd = sigma,log = TRUE)
      likindi[jj] <- dnorm(R20R100_observed[jj], mean = mu, sd = sigma,log = TRUE)
      est_y[jj] <- data[jj]
      cond[jj] <- 0
    }
    
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  liktable <- data.frame(data,likindi,est_y,cond)
  return(-1*lik)
  #return(liktable)
}

# alphaseq <- seq(from = 0.01,to = 1.2, by = 0.01)
# likseq <- c()
# for (iseq in 1:length(alphaseq) ){
#   likseq[iseq] <- dmglik(alphaseq[iseq])
#   
# }
# 
# plot(alphaseq, likseq,type = 'l')
optimize(dmglik, c(0,1.2))
# 
# 
# 
# dmglik(.1)
# dmglik(.12)
# negpdfdmg(.1)
# 
# 
# 
# 
# optimdmg <- optimize(dmglik, c(0,20))
# optimdmg$minimum
# optimCheck::optim_proj(.1,dmglik, xrng = 1)
# 
# # R20R100_data/alpha < l/thresh
# # 
# # R20R100_data[1]
# # log(1/.2*dnorm(R20R100_data[1]/.2, 45,12))
# # 
# # log(1/.3*dnorm(R20R100_data[1]/.3, 47,18))
# # log(1/.1*dnorm(R20R100_data[1]/.1, 47,18))
# 
# 
# # dmglik(c(47,18,.2)) 
# # dmglik(c(47,18,.1)) 
# # dmglik(c(47,18,1)) 
# 
# 
# 
# # dmglik(c(47,18,.8)) 
# # optimout$value
# 
# 
# optimout <- optim(.1, dmglik, method = 'Brent', lower = 0, upper = 19)
# optimout$convergence
# optimout$par
# 
# optimout$value
# 
# 
# 
# optimCheck::optim_proj(optimout$par,dmglik, xrng = 20)
# 
# # write the likelihood for all the groups
# 
# 
# 
# 
# cdfdmg <- function(alpha, ystar){
#   
#   
#   mu <- 47.91
#   sigma <- 18.84
#   
#  cdf <-  pnorm(min(ystar/alpha,l/thresh ),mean = mu, sd = sigma) + max(0,pnorm (ystar ,mean = mu, sd = sigma) - pnorm(l/thresh,mean = mu, sd = sigma))
#   return(cdf)
# }
# 
# cdfdmg(.1, R20R100_data[2])
# 
# xseq <- seq(0,100, by = .1)
# yseq <- c()
# for (jj in 1:length(xseq)){
# yseq[jj] <- cdfdmg(.1, xseq[jj])  
# }
# 
# plot(xseq, yseq , type = 'l')
# 
negpdfdmg <- function(alpha){
  pdf <- c()
  for (jj in 1:length(R20R100_data)){
    pdf[jj] <- (cdfdmg(alpha, R20R100_data[jj] + 0.0001) - cdfdmg(alpha, R20R100_data[jj]))/.0001

  }
  # remove 0
  id <- pdf == 0
  if(sum(id) != 0){
    return(10000)
  }
  else{
    return(-1*sum(log(pdf)))
  }

  # return(-1*sum(log(pdf)))
  # return(pdf)
}

# optimout <- optim(.1,negpdfdmg , method = 'Brent', lower = 0, upper = 19)

sum(negpdfdmg(.8) == 0)
negpdfdmg(.1)
# 
# 
# 
# 
# optimresult <- optimize(negpdfdmg, c(0,2))
# optimCheck::optim_proj(.8,negpdfdmg, xrng = 1)
# optimresult
# negpdfdmg(.2)
# 
# 
# negpdfdmg (.1,R20R100_data)
# 
# 
# 
# 
# negpdfdmg (.3,R20R100_data) - dmglik(.3)
# negpdfdmg (3,R20R100_data) - dmglik(3)
# 
# 
# optimfunction <- function(alpha){
#   negpdfdmg (alpha,R20R100_data)
# }
# 
# optimout <- optim(.1, optimfunction, method = 'Brent', lower = 0, upper = 19)
# 
# optimout$par
# 
# 
# 
# optimCheck::optim_proj(1.5,negpdfdmg, xrng = 20)



x <- rnorm(1000,0,1)
y <- .2*x[which(x>1)]
#y <- .2*x
likx <- function(alpha){
  lik <- 0
  for(i in 1:length(y)){
    lik <- lik + dnorm(y[i]/alpha,log = TRUE) - log(alpha) 
  }
  
  return(-1*lik)
}

alphaseq <- seq(from = .1, to =1, by = 0.01 )
likvalue <- c()
for(i in 1:length(alphaseq)){
  likvalue[i] <- likx(alphaseq[i])
}
plot(alphaseq, likvalue, type = 'l')
optimize(likx, c(0,1))
