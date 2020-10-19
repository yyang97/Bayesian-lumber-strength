##-----------set environment and check--------------------
setwd("F:/study/Research/Bayesian-lumber-strength/R/Single simulation/no damage")
#install.packages("jsonlite", type = "source")
require(rstan)
require(MASS)
require(loo)
devtools::find_rtools()

##-----mu sigma-----####
mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.8

##------proof loading-----####
R_pf <- c(qnorm(0.2,mu[1],sigma[1]), qnorm(0.4,mu[1],sigma[1]),
          qnorm(0.6,mu[1],sigma[1]))
T_pf <- c(qnorm(0.2,mu[2],sigma[2]),qnorm(0.4,mu[2],sigma[2]),
          qnorm(0.6,mu[2],sigma[2]))
N <- 100


##-----------real parameter-------------
real_par <- c(mu,sigma,rho)



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
    PFY_ob[i,2] <- samples[i,2]
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
    PFY_ob[i,2] <- samples[i,2]
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
    PFY_ob[i,2] <- samples[i,2]
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
    PFY_ob[i,2] <- samples[i,2]
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
    PFY_ob[i,2] <- samples[i,2]
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
    PFY_ob[i,2] <- samples[i,2]
    PFY_ob[i,3] <- 0
  }
}

T60_data <- PFY_ob

##-----T100-----######

R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])



##---------------stan fitting-----------####

dmg_mod <- stan_model("damage.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R20 = 1,
       alpha_R40 = 1,alpha_R60 = 1,alpha_T20 = 1,alpha_T40 =1,alpha_R60 = 1 )
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
                    control = list(adapt_delta = 0.8),init = init_dmg)




# result 
print(dmg_fit,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
                       'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))

# pairs
pairs(dmg_fit,pars = c('rho','alpha_R20','alpha_R40',
                       'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))

# estimation and check
estimation_dmg <-  summary(dmg_fit)[[1]][1:5,1]

rbind(estimation_dmg,real_par)
# 95% credible interval 
CI <-  as.numeric(real_par >summary(dmg_fit)[[1]][1:5,4])*
  as.numeric(real_par <summary(dmg_fit)[[1]][1:5,8])

# LOOIC
dmg_loo <- loo(dmg_fit)

# LOOIC 5917.2


##-----------------------without dmg----------------------
nondmg_mod <- stan_model("nondamage.stan")
init_nondmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5)
}
nondmg_fit <- sampling(object = nondmg_mod,
                       data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                   N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                   N_x = length(T100_data),N_y = length(R100_data),
                                   X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                   X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                   t_x = R100_data,t_y = T100_data,
                                   l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                   l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                       control = list(adapt_delta = 0.8),init = init_nondmg)

# result
print(nondmg_fit,pars = c('mu','sigma','rho'))

# estimation 
estimation_nondmg <- summary(nondmg_fit)[[1]][1:5,1]
rbind(estimation_nondmg,real_par)
# LOOIC

nondmg_loo <- loo(nondmg_fit)


# compare LOOIC
## The preferred model will be at the first row. 
loo_compare(dmg_loo,nondmg_loo)
