# Group 1: Y < l
# Group 2: Y>l, Y* < l
# Group 3: Y*>l



## damaged y
y <- rnorm(1000,0,1)
# can only observe x > 1
# y = alpha * x, is the observed data
# how to estimate alpha
ystar <- .8*y[which(y>1)]
ystar <- ystar[which(ystar >1)]

# try to use mle to estimate alpha
liky <- function(alpha){
  lik <- 0
  for(i in 1:length(ystar)){
    lik <- lik + dnorm(ystar[i]/alpha,log = TRUE) - log(alpha) - log(1-pnorm(1/alpha))
  }
  
  # minimize the negative likelihood
  return(-1*lik)
}

alphaseq <- seq(from = .1, to =1, by = 0.01 )
likvalue <- c()
for(i in 1:length(alphaseq)){
  likvalue[i] <- liky(alphaseq[i])
}
plot(alphaseq, likvalue, type = 'l')
optimize(liky, c(0,1))



## experiment for group 2
  ## damaged y
  y <- rnorm(1000,0,1)
  # can only observe x > 1
  # y = alpha * x, is the observed data
  # how to estimate alpha
  g2 <- length(y[which(y>1.2)])
  (1-pnorm(1.2))*1000
  g2
# try to use mle to estimate alpha
liky <- function(alpha){
  (1 - pnorm(alpha) - g2/1000)^2
}

alphaseq <- seq(from = .1, to =5, by = 0.01 )
likvalue <- c()
for(i in 1:length(alphaseq)){
  likvalue[i] <- liky(alphaseq[i])
}
plot(alphaseq, likvalue, type = 'l')
optimize(liky, c(0,5))





## damaged y
y <- rnorm(10000,mean = 20, sd = 10)
# can only observe x > 1
# y = alpha * x, is the observed data
# how to estimate alpha
ystar <- y[which(y<25)]


# try to use mle to estimate alpha
liky <- function(theta){
 mu <- theta[1]
 sigma <- theta[2]
 lik <- sum(dnorm(ystar,mean = mu, sd = sigma, log = TRUE) - log(pnorm(25, mean = mu, sd = sigma)))
  
  # minimize the negative likelihood
  return(-1*lik)
}
#liky(c(20,10))
# alphaseq <- seq(from = .1, to =1, by = 0.01 )
# likvalue <- c()
# for(i in 1:length(alphaseq)){
#   likvalue[i] <- liky(alphaseq[i])
# }
# plot(alphaseq, likvalue, type = 'l')
optim(c(20,10),liky)$par


# undamaged y
y <- rnorm(10000,0,1)
# can only observe x > 1
# y = alpha * x, is the observed data
# how to estimate alpha
ystar <- 1/.8*y[which(y>1 & y<2)] - 0.5
g2g3 <- length(ystar)
ystar <- ystar[which(ystar >1)]
g3 <- length(ystar)
g2 <- g2g3 - g3

10000*(pnorm(1.2) - pnorm(1))
g2


# try to use mle to estimate alpha
liky <- function(alpha){
  lik <- 0
  
  for(i in 1:length(ystar)){
    lik <- lik + dnorm((ystar[i] + .5)*alpha,log = TRUE) + log(alpha) - log(pnorm(2.5*alpha)-pnorm(1.5*alpha))
  }
  lik <- lik -  (1 - pnorm(alpha) - g2/length(y))^2
  # minimize the negative likelihood
  return(-1*lik)
}

alphaseq <- seq(from = .1, to =1, by = 0.01 )
likvalue <- c()
for(i in 1:length(alphaseq)){
  likvalue[i] <- liky(alphaseq[i])
}
plot(alphaseq, likvalue, type = 'l')
optimize(liky, c(0,1))







l <- qnorm(.2, mean =  47.91, sd = 18.84)
thresh <- 0.5
alpha <- .8
samples <- rnorm(10000, mean = 47.91, sd =  18.84)
R20R100_data <- samples[which(samples > l & samples < l/thresh)]
R20R100_observed <- 1/alpha * R20R100_data + l/thresh - 1/alpha*l/thresh
R20R100_observed <- R20R100_observed[which(R20R100_observed>l)]
length(R20R100_observed)

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
    # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
    # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
    # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
    lik <- lik + log(alpha) + dnorm(alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))

    
    
  }

  return(-1*lik)
  #return(liktable)
}

alphaseq <- seq(from = 0.01,to = 1.2, by = 0.01)
likseq <- c()
for (iseq in 1:length(alphaseq) ){
  likseq[iseq] <- dmglik(alphaseq[iseq])

}

plot(alphaseq, likseq,type = 'l')
optimize(dmglik, c(0,1.2))






l <- qnorm(.2, mean =  47.91, sd = 18.84)
thresh <- 0.5
alpha <- .8
samples <- rnorm(1000, mean = 47.91, sd =  18.84)
R20R100_data <- samples[which(samples > l)]
R20R100_observed <- 1/alpha * R20R100_data + l/thresh - 1/alpha*l/thresh
R20R100_observed <- R20R100_observed[which(R20R100_observed>l)]


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
      # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
      lik <- lik + log(alpha) + dnorm(alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(1-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      likindi[jj] <- log(alpha) + dnorm(alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE)- log(1-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      est_y[jj] <- alpha*R20R100_observed[jj]+l/thresh-alpha*l/thresh
      cond[jj] <- 1

    
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  liktable <- data.frame(data,likindi,est_y,cond)
  return(-1*lik)
  #return(liktable)
}

alphaseq <- seq(from = 0.01,to = 1.2, by = 0.01)
likseq <- c()
for (iseq in 1:length(alphaseq) ){
  likseq[iseq] <- dmglik(alphaseq[iseq])
  
}

plot(alphaseq, likseq,type = 'l')
optimize(dmglik, c(0,1.2))





mu <- 47.91
sd <-  18.84
sigma <- sd
alpha <- .8


N <- 35000
samples <- rnorm(N, mean = mu, sd = sd)
l <- qnorm(.2, mean =  mu, sd = sd)

thresh <- 0.5
dataset <- list()


# Group 1: Y < l,216, which is dataset[[1]]
R20_data <- samples[which(samples < l)]
dataset[[1]] <- R20_data

# Group 2: Y>l, Y* < l
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


g2g3 <- cbind(R20R100_observed,damage_table,R20R100_data)


# group3 data
g3 <- g2g3[g2g3[,1] > l,]

# Group 2: Y>l, Y* < l, dataset[[2]] is the number of pieces in group 2
g2 <-  g2g3[g2g3[,1] < l,]
dataset[[2]] <- dim(g2)[1]
# dataset[[2]] <- N - length(dataset[[1]]) - dim(damage_data)[1] 
# dataset[[2]]/N
# check the probalibility of that interval
pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - pnorm(l, mean = mu, sd = sd)


dmglikgroup2 <- function(alpha){
  lik <- sum(dnorm(g2, mean = mu, sd = sigma, log = TRUE) - 
               log(pnorm(l, mean  = mu, sd = sigma)))
  return(-1*lik)
}
optimize(dmglikgroup2,c(0,1))


dmglikgroup2 <- function(alpha){
  pxga <- pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - pnorm(l, mean = mu, sd = sd)
  lik <-  - (pxga - dataset[[2]]/N)^2
  return(-1*lik)
}
# optimize(dmglikgroup2,c(0,1))

# dmglikgroup2 <- function(alpha){
#   pxga <- pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - pnorm(l, mean = mu, sd = sd)
#   lik <-  dbinom(dataset[[2]], size = N, prob = pxga, log = TRUE)
#   return(-1*lik)
# }
optimize(dmglikgroup2,c(0,1))



# Group 3: Y*>l
dataset[[3]] <- g3
  
  
R20R100_data <- dataset[[3]][,1]

# the broken piece in the group 3 
# id <- dataset[[3]][,2] == 1
# R20R100_observed <- dataset[[3]][id,1]



dmglikgroup1 <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  lik <- sum(dnorm(R20_data, mean = mu, sd = sigma, log = TRUE) - log(pnorm(l, mean  = mu, sd = sigma)))
  return(-1*lik)
}

optim(c(50,10), dmglikgroup1)$par

# R20R100_observed <- dataset[[3]][,1]
dmglik <- function(theta,thresh){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  # mu <- 47.91
  # sigma <- 18.84
  # alpha <- theta
  est_y <-c()
  cond <- c()
  # group 1
  lik <- sum(dnorm(R20_data, mean = mu, sd = sigma, log = TRUE) - log(pnorm(l, mean  = mu, sd = sigma)))
  # group 2
  # lik <- lik + dataset[[2]]*log(pnorm(l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma ))
  #lik <- lik + dataset[[2]]*log(pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma ))
  # p(l<x<alpha*l+l/thresh - alpha*l/thresh)
  pxga <- pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma)
  # lik <- lik + dbinom(dataset[[2]], size = N, prob = pxga, log = TRUE)
  lik <- lik + dataset[[2]]*log(pxga)
  likindi <- c()
  for (jj in 1:length(R20R100_data)){
    if(R20R100_data[jj]  < l/thresh){
    # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
    # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
    # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
      
    # conditional likelihood 
    # lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
    # lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE)
    lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE)- log(1-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      
    est_y[jj] <- alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh
    cond[jj] <- 1
    }
    else
    {
      
    # lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)- log(1 - pnorm(l/thresh, mean = mu, sd = sigma))
    lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE) - log(1 - pnorm(l, mean = mu, sd = sigma))
    
    est_y[jj] <- R20R100_data[jj]
    cond[jj] <- 0
    }
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  #liktable <- data.frame(R20R100_data,est_y,cond)
  return(-1*lik)
  #return(liktable)
}

# cbind(dmglik(c(47.91,18.84,.8),thresh = .5),dataset[[3]])

optimout <- optim(c(47.91,18.84,.8),thresh = .5,
                  dmglik, 
                  lower = c(20,10,.001),
                  upper = c(80,30,2), 
                  method = "L-BFGS-B")
optimout$par


# threshseqsim <- seq(from = .1, to = 1.2, by = .01)
# negliksim <- c()
# for (jj in 1:length(threshseqsim) ){
#   optimout <- optim(c(47.91,18.84,.8),thresh = threshseqsim[jj],
#                                    dmglik,
#                                    lower = c(20,10,.001),
#                                    upper = c(80,30,2),
#                                    method = "L-BFGS-B")
#   negliksim[jj] <- optimout$value
# }
# 
# plot(threshseqsim, negliksim, type = 'l')


# dmglik(theta =c(47.91,18.84,.8), thresh = .5)
# dmglik(theta = optimout$par, thresh = .5)


# 
# 
alphaseq <- seq(from = 0.1,to = 1.1, by = 0.01)
likseq <- c()
for (iseq in 1:length(alphaseq) ){
  likseq[iseq] <- dmglik(c(47.91,18.84,alphaseq[iseq]),thresh = .5)

}

plot(alphaseq, likseq,type = 'l')
# optimize(dmglik, c(0,1.2))




## specify the strength of pieces in group 2 


mu <- 47.91
sd <-  18.84
sigma <- 18.84
alpha <- .8


N <- 35000
samples <- rnorm(N, mean = mu, sd = sd)
l <- qnorm(.2, mean =  mu, sd = sd)

thresh <- 0.5
dataset <- list()


# Group 1: Y < l,216, which is dataset[[1]]
R20_data <- samples[which(samples < l)]
dataset[[1]] <- R20_data

# Group 2: Y>l, Y* < l
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


g2g3 <- cbind(R20R100_observed,damage_table,R20R100_data)


# group3 data
g3 <- g2g3[g2g3[,1] > l,]

# Group 2: Y>l, Y* < l, dataset[[2]] is the number of pieces in group 2
g2 <-  g2g3[g2g3[,1] < l,]
# dataset[[2]] <- dim(g2)[1]
dataset[[2]] <- g2[,1]  


max(alpha*dataset[[2]]+l/thresh-alpha*l/thresh)
alpha*l + l/thresh - alpha*l/thresh



# dataset[[2]] <- N - length(dataset[[1]]) - dim(damage_data)[1] 
# length(dataset[[2]])/N
# check the probalibility of that interval
pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - pnorm(l, mean = mu, sd = sd)
# dmglikgroup2 <- function(alpha){
#   pxga <- pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - pnorm(l, mean = mu, sd = sd)
#   lik <-  - (pxga - dataset[[2]]/N)^2
#   return(-1*lik)
# }
# optimize(dmglikgroup2,c(0,1))

dmglikgroup2 <- function(alpha){
  # pxga <- pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - pnorm(l, mean = mu, sd = sd)
  # lik <-  dbinom(dataset[[2]], size = N, prob = pxga, log = TRUE)
  # return(-1*lik)
  
  lik <- 0
  for (jj in 1:length(dataset[[2]])){
    lik <- lik + log(alpha) + 
      dnorm(alpha*dataset[[2]][jj]+l/thresh-alpha*l/thresh,
            mean =  mu,sd = sd, log = TRUE) - 
      log(pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sd) - 
            pnorm(l, mean = mu, sd = sd))
  }
  return(-1*lik)
}
optimize(dmglikgroup2,c(0,1))



max(alpha*dataset[[2]]+l/thresh-alpha*l/thresh)
alpha*l + l/thresh - alpha*l/thresh
optimize(dmglikgroup2,c(0,1))

dmglikgroup2(.95)

# Group 3: Y*>l
dataset[[3]] <- g3


R20R100_data <- dataset[[3]][,1]

# the broken piece in the group 3 
# id <- dataset[[3]][,2] == 1
# R20R100_observed <- dataset[[3]][id,1]



dmglikgroup1 <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  lik <- sum(dnorm(R20_data, mean = mu, sd = sigma, log = TRUE) - log(pnorm(l, mean  = mu, sd = sigma)))
  return(-1*lik)
}

optim(c(50,10), dmglikgroup1)$par

# R20R100_observed <- dataset[[3]][,1]
dmglik <- function(theta,thresh){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  # mu <- 47.91
  # sigma <- 18.84
  # alpha <- theta
  est_y <-c()
  cond <- c()
  # group 1
  lik <- sum(dnorm(R20_data, mean = mu, sd = sigma, log = TRUE) - log(pnorm(l, mean  = mu, sd = sigma)))
  # group 2
  # lik <- lik + dataset[[2]]*log(pnorm(l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma ))
  #lik <- lik + dataset[[2]]*log(pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma ))
  # p(l<x<alpha*l+l/thresh - alpha*l/thresh)
  
  
  # pxga <- pnorm(alpha*l+l/thresh - alpha*l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma)
  #  lik <- lik + dbinom(dataset[[2]], size = N, prob = pxga, log = TRUE)
  # for (jj in 1:length(dataset[[2]])){
  #   lik <- lik + log(alpha) + dnorm(alpha*dataset[[2]][jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE) - log(pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma) - pnorm(l, mean = mu, sd = sigma))
  # }
   
   
   
  # lik <- lik + dataset[[2]]*log(pxga)
  likindi <- c()
  for (jj in 1:length(R20R100_data)){
    if(R20R100_data[jj]  < l/thresh){
      # lik <- lik - log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # likindi[jj] <- -1*log(alpha) + dnorm((R20R100_observed[jj]+ alpha*l/thresh - l/thresh)/alpha, mean =  mu,sd = sigma, log = TRUE)
      # est_y[jj] <- (data[jj]+ alpha*l/thresh - l/thresh)/alpha
      
      # conditional likelihood 
       lik <- lik + 
         log(alpha) + 
         dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, 
               mean =  mu,sd = sigma, log = TRUE) - 
         log(pnorm(l/thresh, 
                   mean = mu, sd = sigma)-
               pnorm(alpha*l + l/thresh - alpha*l/thresh, 
                     mean = mu, sd = sigma))
      # lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE)
      # lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, mean =  mu,sd = sigma, log = TRUE)- log(1-pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma))
      
      # lik <- lik + log(alpha) + dnorm(alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh, 
      #                                 mean =  mu,sd = sigma, log = TRUE) +  
      #   log(pnorm(l/thresh, mean = mu, sd = sigma)-
      #         pnorm(alpha*l + l/thresh - alpha*l/thresh, mean = mu, sd = sigma)) - 
      #   log(pnorm(l/thresh, mean = mu, sd = sigma)-pnorm(l, mean = mu, sd = sigma))
      # 
      est_y[jj] <- alpha*R20R100_data[jj]+l/thresh-alpha*l/thresh
      cond[jj] <- 1
    }
    else
    {
      
       lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)- log(1 - pnorm(l/thresh, mean = mu, sd = sigma))
       # lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)
      
      est_y[jj] <- R20R100_data[jj]
      cond[jj] <- 0
    }
  }
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
  #liktable <- data.frame(R20R100_data,est_y,cond)
  return(-1*lik)
  #return(liktable)
}

# cbind(dmglik(c(47.91,18.84,.8),thresh = .5),dataset[[3]])

optimout <- optim(c(47.91,18.84,.8),thresh = .5,
                  dmglik, 
                  lower = c(20,10,.001),
                  upper = c(80,30,2), 
                  method = "L-BFGS-B")
optimout$par





