require(rstan)
require(MASS)
require(loo)
require(optimCheck)



mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.8


N <- 300

# alpha[6] alpha for R20, R40,R60, T20,T40,T60
# alpha <- c(.1,.2,.3,.4,.5,.6) 


# alpha <- c(.8,.95,.85,.8,.95,.85) 
alpha <- c(.8,.85,.9,.8,.85,.9) 

prop <- 0.7


##------proof loading-----####
R_pf <- c(qnorm(0.2,mu[1],sigma[1]), qnorm(0.4,mu[1],sigma[1]),
          qnorm(0.6,mu[1],sigma[1]))
T_pf <- c(qnorm(0.2,mu[2],sigma[2]),qnorm(0.4,mu[2],sigma[2]),
          qnorm(0.6,mu[2],sigma[2]))




###---------R20---------#####

sd_x <- sigma[1]
sd_y <- sigma[2]
Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
bvn1 <- mvrnorm(N, mu = c(mu[1],mu[2]), Sigma = Sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X","bvn1_Y")
samples <- bvn1

sample_size <- nrow(samples)
N <- sample_size
l <- R_pf[1]
PFY_ob <- matrix(0, nrow = N, ncol = 3)
for (i in 1:nrow(PFY_ob)){
  if(samples[i,1] <= l){
    PFY_ob[i,1] <- samples[i,1]
    PFY_ob[i,3] <- 1
  }
  else{
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*prop, alpha[1],1)
    PFY_ob[i,3] <- 0
  }
}

R20_data <- PFY_ob

##-----------R40------------#####

sd_x <- sigma[1]
sd_y <- sigma[2]
Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
bvn1 <- mvrnorm(N, mu = c(mu[1],mu[2]), Sigma = Sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X","bvn1_Y")
samples <- bvn1

sample_size <- nrow(samples)
N <- sample_size
l <- R_pf[2]
PFY_ob <- matrix(0, nrow = N, ncol = 3)
for (i in 1:nrow(PFY_ob)){
  if(samples[i,1] <= l){
    PFY_ob[i,1] <- samples[i,1]
    PFY_ob[i,3] <- 1
  }
  else{
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*prop, alpha[2],1)
    PFY_ob[i,3] <- 0
  }
}

R40_data <- PFY_ob


##--------------R60-------------####
sd_x <- sigma[1]
sd_y <- sigma[2]
Sigma<- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
bvn1 <- mvrnorm(N, mu = c(mu[1],mu[2]), Sigma = Sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X","bvn1_Y")
samples <- bvn1

sample_size <- nrow(samples)
N <- sample_size
l <- R_pf[3]
PFY_ob <- matrix(0, nrow = N, ncol = 3)
for (i in 1:nrow(PFY_ob)){
  if(samples[i,1] <= l){
    PFY_ob[i,1] <- samples[i,1]
    PFY_ob[i,3] <- 1
  }
  else{
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*prop, alpha[3],1)
    PFY_ob[i,3] <- 0
  }
}

R60_data <- PFY_ob


##---------T20-------####
sd_x <- sigma[2]
sd_y <- sigma[1]
Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
bvn1 <- mvrnorm(N, mu = c(mu[2],mu[1]), Sigma = Sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X","bvn1_Y")
samples <- bvn1

sample_size <- nrow(samples)
N <- sample_size
l <- T_pf[1]
PFY_ob <- matrix(0, nrow = N, ncol = 3)
for (i in 1:nrow(PFY_ob)){
  if(samples[i,1] <= l){
    PFY_ob[i,1] <- samples[i,1]
    PFY_ob[i,3] <- 1
  }
  else{
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*prop, alpha[4],1)
    PFY_ob[i,3] <- 0
  }
}

T20_data <- PFY_ob

##---------T40-------####
sd_x <- sigma[2]
sd_y <- sigma[1]
Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
bvn1 <- mvrnorm(N, mu = c(mu[2],mu[1]), Sigma = Sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X","bvn1_Y")
samples <- bvn1

sample_size <- nrow(samples)
N <- sample_size
l <- T_pf[2]
PFY_ob <- matrix(0, nrow = N, ncol = 3)
for (i in 1:nrow(PFY_ob)){
  if(samples[i,1] <= l){
    PFY_ob[i,1] <- samples[i,1]
    PFY_ob[i,3] <- 1
  }
  else{
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*prop, alpha[5],1)
    PFY_ob[i,3] <- 0
  }
}

T40_data <- PFY_ob

##------------T60--------------###

sd_x <- sigma[2]
sd_y <- sigma[1]
Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
bvn1 <- mvrnorm(N, mu = c(mu[2],mu[1]), Sigma = Sigma ) # from MASS package
colnames(bvn1) <- c("bvn1_X","bvn1_Y")
samples <- bvn1

sample_size <- nrow(samples)
N <- sample_size
l <- T_pf[3]
PFY_ob <- matrix(0, nrow = N, ncol = 3)
for (i in 1:nrow(PFY_ob)){
  if(samples[i,1] <= l){
    PFY_ob[i,1] <- samples[i,1]
    PFY_ob[i,3] <- 1
  }
  else{
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*prop, alpha[6],1)
    PFY_ob[i,3] <- 0
  }
}

T60_data <- PFY_ob

##-----T100-----######

R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])




options(mc.cores = parallel::detectCores())
dmg_mod <- stan_model("newmodel_test.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5,
       alpha = rep(1,6),thresh = 0.9 )
}
dmg_fit <- sampling(object = dmg_mod,
                    data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                N_x = length(T100_data),N_y = length(R100_data),
                                X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                t_x = R100_data,t_y = T100_data,
                                l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                    init = init_dmg)
print(dmg_fit)


dmg_fit <- readRDS('dmg_fit_thresh.rds')
print(dmg_fit)


# real data 
dmg_fit <- readRDS('real_fit_thresh.rds')
print(dmg_fit)



saveRDS(dmg_fit, file = 'dmg_fit_thresh.rds')

saveRDS(dmg_fit, file = 'real_fit_thresh.rds')





yseq <- seq(from = 1, to = 100, by =.1)
ystarseq <- ((1/0.2)*yseq + 32/.8 - 1/0.2*32/.8)*ifelse(.8*yseq<32 & yseq > 32,1,0) + yseq*ifelse(.8*yseq>=32 | yseq < 32,1,0)

plot(yseq,ystarseq, type = 'l', ylim =c(0,100))
# the dash line is no damaged plot 
abline(yseq, yseq, lty = 2)



yseq <- rnorm(100000,mean = 48, sd = 19)
plot(density(yseq))
plot(ecdf(yseq))


yseq <- rnorm(100000,mean = 48, sd = 19)
ystarseq <- ((1/0.5)*yseq + 32/.8 - 1/0.5*32/.8)*ifelse(.8*yseq<32 & yseq > 32,1,0) + yseq*ifelse(.8*yseq>=32 | yseq < 32,1,0)
plot(density(ystarseq))
abline(v = 32)
abline(v = 32/.8)
plot(ecdf(ystarseq))
abline(v = 32)
abline(v = 32/.8)

ystarcdf <- ecdf(ystarseq)
ystarcdf[2]
## cdf of ystar
x <- seq(from = 0, to = 100, by = .01)
xcdf <- pnorm(x, mean = 48, sd = 19)
id <- which(x >32 &x<32/.8)
xcdf[id] <- pnorm(x[id]*.5+32/.8 - 32/.8*.5, mean = 48,sd = 19)

plot(x,xcdf,type = 'l')

alpha <- .5
c <- .8
l <- 32
xcdf <- c()
for(ii in 1:length(x)){
  if(x[ii] < 1/alpha*l + l/c - 1/alpha*l/c){
    xcdf[ii] <- pnorm(x[ii], mean = 48, sd = 19)
  }
  else if(x[ii] > 1/alpha*l + l/c - 1/alpha*l/c & x[ii] <l){
    xcdf[ii] <- pnorm(x[ii], mean = 48, sd = 19) + 
      pnorm(alpha*x[ii] + l/c - alpha*l/c, mean = 48, sd = 19) - 
      pnorm(l, mean = 48, sd = 19)
  }
  else if(x[ii] > l & x[ii] < l/c){
    xcdf[ii] <- pnorm(alpha*x[ii] + l/c - alpha*l/c, mean = 48, sd = 19)  
  }
  else if(x[ii] > l/c){
    xcdf[ii] <- pnorm(x[ii], mean = 48, sd = 19)
  }
  
}



## pdf 
xpdf <- c()
for(ii in 1:length(x)){
  if(x[ii] < 1/alpha*l + l/c - 1/alpha*l/c){
    xpdf[ii] <- dnorm(x[ii], mean = 48, sd = 19)
  }
  else if(x[ii] > 1/alpha*l + l/c - 1/alpha*l/c & x[ii] <l){
    xpdf[ii] <- dnorm(x[ii], mean = 48, sd = 19) + 
      alpha*dnorm(alpha*x[ii] + l/c - alpha*l/c, mean = 48, sd = 19) 
  }
  else if(x[ii] > l & x[ii] < l/c){
    xpdf[ii] <- alpha*dnorm(alpha*x[ii] + l/c - alpha*l/c, mean = 48, sd = 19)  
  }
  else if(x[ii] > l/c){
    xpdf[ii] <- dnorm(x[ii], mean = 48, sd = 19)
  }
  
}
plot(x,xpdf,type = 'l')


dmglik <- function(x, alpha){
  mu <- 48
  sigma <- 19
  lik <- 0
  
  id1 <-  x < 1/alpha*l + l/c - 1/alpha*l/c
  lik <- sum(dnorm(x[id1], mean = mu, sd = sigma, log = TRUE))
  
  id2 <- x> 1/alpha*l + l/c - 1/alpha*l/c & x <l
  lik <- lik + sum(log(dnorm(x[id2], mean = mu, sd = sigma) + 
    alpha*dnorm(alpha*x[id2] + l/c - alpha*l/c, mean = mu, sd = sigma)))
  
  id3 <- x > l & x < l/c
  lik <- lik + sum(log(alpha*dnorm(alpha*x[id3] + l/c - alpha*l/c, mean = mu, sd = sigma)))
  
  id4 <- x > l/c
  lik <- lik + sum(dnorm(x[id4], mean = mu, sd = sigma, log = TRUE))
  
  return(-lik)
  
}


alphaseq <- seq(from = .01, to = .99, by = 0.01)
likvalue <- c()
for (jj in 1:length(alphaseq)){
  likvalue[jj] <- dmglik(ystarseq, alphaseq[jj])
}

plot(alphaseq, likvalue, type = 'l')



dmglik <- function(x, theta){
  mu <- theta[1] 
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  lik <- 0
  
  id1 <-  x < 1/alpha*l + l/c - 1/alpha*l/c
  lik <- sum(dnorm(x[id1], mean = mu, sd = sigma, log = TRUE))
  
  id2 <- x> 1/alpha*l + l/c - 1/alpha*l/c & x <l
  lik <- lik + sum(log(dnorm(x[id2], mean = mu, sd = sigma) + 
                         alpha*dnorm(alpha*x[id2] + l/c - alpha*l/c, mean = mu, sd = sigma)))
  
  id3 <- x > l & x < l/c
  lik <- lik + sum(log(alpha*dnorm(alpha*x[id3] + l/c - alpha*l/c, mean = mu, sd = sigma)))
  
  id4 <- x > l/c
  lik <- lik + sum(dnorm(x[id4], mean = mu, sd = sigma, log = TRUE))
  
  return(-lik)
  
}

optim(c(48,19,.5,.5),dmglik, x = ystarseq,
      lower = c(20,10,.001,.001),
      upper = c(80,30,2,2), 
      method = "L-BFGS-B")$par





meanreg_logel <- function(y, X, beta) {
  yXb <- y - X %*% beta
  G <- as.numeric(yXb) * X
  flexEL::logEL(G = G, supp_adj = TRUE)
}

# the true value of beta
beta_true <- c(1,2)
# the sample size n = 200
n <- 200
# simulate the covariate X = c(1,x)
X <- cbind(1, rnorm(n))
# simulate the noise term epsilon
eps <- rnorm(n)
#simulate the response data 
y <- c(X %*% beta_true + eps)



meanreg_logel(y,X, c(2,2))
neg_meanreg_logel <- function(beta) {
  -meanreg_logel(y,X,beta)
}
nlm(neg_meanreg_logel, c(0.75, 1.25))


?invgamma::rinvgamma

x <- seq(0, 1, length.out = 101)
# good choice 
dx <- dbeta(x, 2, 2)
plot(x,dx, type = 'l')


## another choice 
dx <- dbeta(x, 3, 3)
plot(x,dx, type = 'l')


plotInvGamma(shape, scale)
