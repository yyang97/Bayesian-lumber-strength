##-----------set environment and check--------------------
#setwd("~/Downloads/github/Bayesian-lumber-strength/R/Simulation")
require(rstan)
require(MASS)
require(loo)



options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


stanresult <- list()
for (nsim in 1:100){
  
  print(nsim)
##-----mu sigma-----####
mu <- c(44.73,5.49)
sigma <- c(12.74,1.05)

rho <- 0.8

##------proof loading-----####
R_pf <- c(qnorm(0.2,mu[1],sigma[1]), qnorm(0.4,mu[1],sigma[1]),
          qnorm(0.6,mu[1],sigma[1]))
T_pf <- c(qnorm(0.2,mu[2],sigma[2]),qnorm(0.4,mu[2],sigma[2]),
          qnorm(0.6,mu[2],sigma[2]))
N <- 87


##--------------alpha-----#####

# alpha_R <- 10* c(0.67,1.83,1.29)
# alpha_T <- 10*c(85.23,104.09,65.95)

#alpha_R <- 10* c(0.67,1.83,1.29)
#alpha_T <- 10*c(85.23,104.09,65.95)
 alpha_R <- c(0,0,0)
 alpha_T <- c(0,0,0)
##-----------real parameter-------------
real_par <- c(rho,alpha_R,alpha_T)



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

# damage data
alpha <- alpha_R[1]
id_y <- which(PFY_ob[,3] == 0)
R20_data <- PFY_ob
R20_data[id_y,2] <- PFY_ob[id_y,2] - alpha/PFY_ob[id_y,2]


remove(PFY_ob)





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

# damage data
alpha <- alpha_R[2]
id_y <- which(PFY_ob[,3] == 0)
R40_data <- PFY_ob
R40_data[id_y,2] <- PFY_ob[id_y,2] - alpha/PFY_ob[id_y,2]


remove(PFY_ob)



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

# damage data
alpha <- alpha_R[3]
id_y <- which(PFY_ob[,3] == 0)
R60_data <- PFY_ob
R60_data[id_y,2] <- PFY_ob[id_y,2] - alpha/PFY_ob[id_y,2]


remove(PFY_ob)




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

# damage data
alpha <- alpha_T[1]
id_y <- which(PFY_ob[,3] == 0)
T20_data <- PFY_ob
T20_data[id_y,2] <- PFY_ob[id_y,2] - alpha/PFY_ob[id_y,2]


remove(PFY_ob)




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

# damage data
alpha <- alpha_T[2]
id_y <- which(PFY_ob[,3] == 0)
T40_data <- PFY_ob
T40_data[id_y,2] <- PFY_ob[id_y,2] - alpha/PFY_ob[id_y,2]


remove(PFY_ob)




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

# damage data
alpha <- alpha_T[3]
id_y <- which(PFY_ob[,3] == 0)
T60_data <- PFY_ob
T60_data[id_y,2] <- PFY_ob[id_y,2] - alpha/PFY_ob[id_y,2]


remove(PFY_ob)




##-----T100-----######

R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])



##---------------stan fitting-----------####

dmg_mod <- stan_model("damage_real.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R20 = 1,
       alpha_R40 = 1,alpha_R60 = 1,alpha_T20 = 1,alpha_T40 =1,alpha_R60 = 1 )
}

data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
            N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
            N_x = length(T100_data),N_y = length(R100_data),
            X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
            X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
            t_x = R100_data,t_y = T100_data,
            l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
            l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3])
dmg_fit <- sampling(object = dmg_mod,
                    data = data,
                    control = list(adapt_delta = 0.8),init = init_dmg)

stanresult[[nsim]] <- summary(dmg_fit)[[1]]

saveRDS(stanresult,file = "stanresult.rds")
}
#dmg_fit <- readRDS("dmg_fit.rds") 

# traceplot(dmg_fit,pars = c('rho','alpha_R20','alpha_R40',
#                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
# result
# print(dmg_fit,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
#                       'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))

# pairs
# pairs(dmg_fit,pars = c('rho','alpha_R20','alpha_R40',
#                        'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
# 
# loo_fit <- loo(dmg_fit)


#waic(dmg_fit)
# dmg_fit_onechain <- readRDS("onechain.rds")
# print(dmg_fit_onechain,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
#                        'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
# 
# 
# a <- summary(dmg_fit_onechain,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
#                                      'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
#saveRDS(dmg_fit, file = "onechain.rds")
#hist(dmg_alpha[[6]])
# # estimation and check
# estimation_dmg <-  summary(dmg_fit)[[1]][5:11,1]
# 
# rbind(estimation_dmg,real_par)
# # 95% credible interval 
# CI <-  as.numeric(real_par >summary(dmg_fit)[[1]][5:11,4])*
#   as.numeric(real_par <summary(dmg_fit)[[1]][5:11,8])
# 
# # LOOIC
# dmg_loo <- loo(dmg_fit)
# 
# # LOOIC 5917.2
# 
# 17066.1
##-----------------------without dmg----------------------
# nondmg_mod <- stan_model("nondamage.stan")
# init_nondmg <- function() {
#   list(mu = c(35,8), sigma = c(10,1), rho = .5)
# }
# nondmg_fit <- sampling(object = nondmg_mod,
#                        data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
#                                    N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
#                                    N_x = length(T100_data),N_y = length(R100_data),
#                                    X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
#                                    X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
#                                    t_x = R100_data,t_y = T100_data,
#                                    l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
#                                    l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
#                        control = list(adapt_delta = 0.8),init = init_nondmg)
# 
# # result
# print(nondmg_fit,pars = c('mu','sigma','rho'))
# pairs(nondmg_fit,pars = c('mu','sigma','rho'))
# # estimation 
# estimation_nondmg <- summary(nondmg_fit)[[1]][5,1]
# rbind(estimation_nondmg,real_par[1])
# # LOOIC,17068.9
# 
# nondmg_loo <- loo(nondmg_fit)
# 
# 
# # compare LOOIC
# ## The preferred model will be at the first row. 
# loo_compare(dmg_loo,nondmg_loo)
# 
# 
# ##----------------model with only R40----------------------------------
# R40dmg_mod <- stan_model("only_alphaR40.stan")
# init_R40dmg <- function() {
#   list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R40 = 1)
# }
# 
# 
# R40dmg_fit <- sampling(object = R40dmg_mod,
#                        data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
#                                    N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
#                                    N_x = length(T100_data),N_y = length(R100_data),
#                                    X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
#                                    X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
#                                    t_x = R100_data,t_y = T100_data,
#                                    l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
#                                    l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
#                        control = list(adapt_delta = 0.8),init = init_R40dmg)
# print(R40dmg_fit,pars = c('mu','sigma','rho','alpha_R40'))
# pairs(R40dmg_fit,pars = c('mu','sigma','rho','alpha_R40'))
# # LOOIC,17060.1
# 
# loo_R40dmg <- loo(R40dmg_fit)
# loo_compare(loo_R40dmg ,nondmg_loo)
# 
# 
# 
# 
# ##-------------------Model 4 let alpha_R40 be real----------------
# R40dmgreal_mod <- stan_model("only_alphaR40_real.stan")
# init_R40dmg <- function() {
#   list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R40 = 1)
# }
# 
# 
# R40dmgreal_fit <- sampling(object = R40dmgreal_mod,
#                            data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
#                                        N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
#                                        N_x = length(T100_data),N_y = length(R100_data),
#                                        X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
#                                        X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
#                                        t_x = R100_data,t_y = T100_data,
#                                        l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
#                                        l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
#                            control = list(adapt_delta = 0.8),init = init_R40dmg)
# print(R40dmgreal_fit,pars = c('mu','sigma','rho','alpha_R40'))
# pairs(R40dmgreal_fit,pars = c('mu','sigma','rho','alpha_R40'))
# 
# 
# # LOOIC, 17060.0
# loo_R40dmgrel <- loo(R40dmgreal_fit)
# loo_compare(loo_R40dmg, loo_R40dmgrel)
# 



# ##-------------------Model 5 let all alpha's be real----------------
# dmgreal_mod <- stan_model("damage_real.stan")
# 
# init_dmg <- function() {
#   list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R20 = 1,
#        alpha_R40 = 1,alpha_R60 = 1,alpha_T20 = 1,alpha_T40 =1,alpha_R60 = 1 )
# }
# 
# dmgreal_fit <- sampling(object = dmgreal_mod,
#                         data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
#                                     N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
#                                     N_x = length(T100_data),N_y = length(R100_data),
#                                     X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
#                                     X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
#                                     t_x = R100_data,t_y = T100_data,
#                                     l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
#                                     l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
#                         control = list(adapt_delta = 0.8),init = init_dmg)
# 
# print(dmgreal_fit,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
#                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
# 
# 
# pairs(dmgreal_fit,pars = c('alpha_R20','alpha_R40',
#                           'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
# traceplot(dmgreal_fit,pars = c('alpha_R20','alpha_R40',
#                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
# 
# 

stanresult <- readRDS("stanresult.rds")
Rhat <- matrix(NA,nrow = 12, ncol = 100)
for (ii in 1:dim(Rhat)[2]){
  Rhat[,ii] <- stanresult[[ii]][,10]
  
}

which.max(Rhat)

which(colMeans(Rhat) > 1.1)

Rhat[,98]
# run longer if unstable 
# 
# > which(colMeans(Rhat) > 1.1)
# [1] 39 51 62 63 73 88 98
# > Rhat[,39]
# [1] 1.002134 1.008719 1.031091 1.031661 1.612425 1.215041 1.341164 1.386088 1.162415
# [10] 1.211965 1.315456 1.002185
# > Rhat[,98]
# [1] 1.527733 1.123418 1.899453 1.154466 1.567547 1.110722 1.294807 1.196106 1.125926
# [10] 1.341846 1.173881 1.365159
