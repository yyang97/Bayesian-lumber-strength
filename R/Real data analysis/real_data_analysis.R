require(rstan)
require(MASS)
require(loo)
library(readxl)
library(dplyr)
library(moments)
###------------------data preprocessing--------------
#setwd("~/Downloads/github/Bayesian-lumber-strength/R/Real data analysis")

summary_all08122013 <- read_excel("summary_all08122013.xlsx")

summary_all08122013
summary_all08122013 <- select(summary_all08122013,"Group","Broken","MOR","UTS")

## T group 
T100_data <-summary_all08122013[which(summary_all08122013$Group == "T100"),]
T60_data <-summary_all08122013[which(summary_all08122013$Group == "T60"),]
T40_data <-summary_all08122013[which(summary_all08122013$Group == "T40"),]
T20_data <-summary_all08122013[which(summary_all08122013$Group == "T20"),]
T100_data$UTS <- as.numeric(T100_data$UTS)
T100_data <- T100_data$UTS
## G group 
R100_data <-summary_all08122013[which(summary_all08122013$Group == "R100"),]
R60_data <-summary_all08122013[which(summary_all08122013$Group == "R60"),]
R40_data <-summary_all08122013[which(summary_all08122013$Group == "R40"),]
R20_data <-summary_all08122013[which(summary_all08122013$Group == "R20"),]
R100_data <-as.numeric(R100_data$MOR)
R100_data 
##-----------------T substitute NA to 0 -----------

id <- T60_data$Broken == 1
T60_data$MOR[id] <- 0
T60_data$UTS[id] <- as.numeric(T60_data$UTS[id])
T60_data$UTS[!id] <- 0
T60_data$UTS <- as.numeric(T60_data$UTS)
T60_data$MOR <- as.numeric(T60_data$MOR)


id <- T40_data$Broken == 1
T40_data$MOR[id] <- 0
T40_data$UTS[id] <- as.numeric(T40_data$UTS[id])
T40_data$UTS[!id] <- 0
T40_data$UTS <- as.numeric(T40_data$UTS)
T40_data$MOR <- as.numeric(T40_data$MOR)


id <- T20_data$Broken == 1
T20_data$MOR[id] <- 0
T20_data$UTS[id] <- as.numeric(T20_data$UTS[id])
T20_data$UTS[!id] <- 0
T20_data$UTS <- as.numeric(T20_data$UTS)
T20_data$MOR <- as.numeric(T20_data$MOR)



##-----------------R substitute NA to 0 ------

id <- R60_data$Broken == 0
R60_data$MOR[id] <- 0
R60_data$UTS[id] <- as.numeric(R60_data$UTS[id])
R60_data$UTS[!id] <- 0
R60_data$UTS <- as.numeric(R60_data$UTS)
R60_data$MOR <- as.numeric(R60_data$MOR)



id <- R40_data$Broken == 0
R40_data$MOR[id] <- 0
R40_data$UTS[id] <- as.numeric(R40_data$UTS[id])
R40_data$UTS[!id] <- 0
R40_data$UTS <- as.numeric(R40_data$UTS)
R40_data$MOR <- as.numeric(R40_data$MOR)

id <- R20_data$Broken == 0
R20_data$MOR[id] <- 0
R20_data$UTS[id] <- as.numeric(R20_data$UTS[id])
R20_data$UTS[!id] <- 0
R20_data$UTS <- as.numeric(R20_data$UTS)
R20_data$MOR <- as.numeric(R20_data$MOR)




##---------------Convert psi to Mpa------------ 

# 1 thousand psi = 6.895 MPa

c = 6.895

R20_data$MOR <- c*R20_data$MOR
R20_data$UTS <- c*R20_data$UTS


R40_data$MOR <- c*R40_data$MOR
R40_data$UTS <- c*R40_data$UTS

R60_data$MOR <- c*R60_data$MOR
R60_data$UTS <- c*R60_data$UTS

R100_data <- c*R100_data

T20_data$MOR <- c*T20_data$MOR
T20_data$UTS <- c*T20_data$UTS

T40_data$MOR <- c*T40_data$MOR
T40_data$UTS <- c*T40_data$UTS

T60_data$MOR <- c*T60_data$MOR
T60_data$UTS <- c*T60_data$UTS


T100_data <- c*T100_data

##-------check normal fitting---------
# 
# library(fitdistrplus)
# 
# # check R100_data
# FIT <- fitdist(R100_data, "norm")    ## note: it is "norm" not "normal"
# plot(FIT)    ## use method `plot.fitdist`
# FIT$estimate # mean = 45.679 sd = 12.900
# FIT$bic # bic = 1394.022
# # good normal fitting 
# 
# 
# # check T100_data
# 
# shapiro.test(T100_data) # p = 0.0001841
# shapiro.test(log(T100_data)) # p-value = 0.002327
# shapiro.test(sqrt(T100_data)) # p-value = 0.3988
# 


##------proof loading-------------- 
R_pf <- c* c(4.956690733, 6.110714122, 7.092435407)
T_pf <- sqrt(c* c(2.962390379, 3.986497991, 4.916102264))


##--------convert UTS to sqrt(UTS)

R20_data$UTS <- sqrt(R20_data$UTS)
R40_data$UTS <- sqrt(R40_data$UTS)
R60_data$UTS <- sqrt(R60_data$UTS)
R100_data <- R100_data


T20_data$UTS <- sqrt(T20_data$UTS)
T40_data$UTS <- sqrt(T40_data$UTS)
T60_data$UTS <- sqrt(T60_data$UTS)
T100_data <- sqrt(T100_data)



##-----check outlier and delete it-------

# R20
id <- which((R20_data$MOR < R_pf[1]) == F)
R20_data <- R20_data[-id,]

# R40
id <- which((R40_data$MOR < R_pf[2]) == F)

# R60
id <- which((R60_data$MOR < R_pf[3]) == F)
R60_data <- R60_data[-id,]

# T20
id <- which((T20_data$UTS < T_pf[1]) == F)
T20_data <- T20_data[-id,]

# T40
id <- which((T40_data$UTS < T_pf[2]) == F)
T40_data <- T40_data[-id,]

# T60
id <- which((T60_data$UTS < T_pf[3]) == F)
T60_data <- T60_data[-id,]


##---- Convert the data so stan can use it -----
R20_data <- cbind(R20_data$MOR,R20_data$UTS,R20_data$Broken)
R40_data <- cbind(R40_data$MOR,R40_data$UTS,R40_data$Broken)
R60_data <- cbind(R60_data$MOR,R60_data$UTS,R60_data$Broken)
T20_data <- cbind(T20_data$UTS,T20_data$MOR,T20_data$Broken)
T40_data <- cbind(T40_data$UTS,T40_data$MOR,T40_data$Broken)
T60_data <- cbind(T60_data$UTS,T60_data$MOR,T60_data$Broken)



##------data preprocessing completed -----


##------ Stan for all alpha-------------
dmg_mod <- stan_model("damage.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R20 = 1,
       alpha_R40 = 1,alpha_R60 = 1,alpha_T20 = 1,alpha_T40 =1,alpha_R60 = 1 )
}
set.seed(2020)
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

print(dmg_fit,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))


#summary(dmg_fit,pars = 'rho')$summary[,c('mean','2.5%','97.5%')]


pairs(dmg_fit,pars = c('rho','alpha_R20','alpha_R40',
                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))

dmg_alpha <- extract(dmg_fit,pars = c('alpha_R20','alpha_R40',
                                     'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))
dmg_wood <- extract(dmg_fit,pars = c('mu','sigma','rho'))
sapply(dmg_wood, moments::skewness)
sapply(dmg_alpha, moments::skewness)


loo_dmg <- loo(dmg_fit)
##-------without alpha---------

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
print(nondmg_fit,pars = c('mu','sigma','rho'))
pairs(nondmg_fit,pars = c('mu','sigma','rho'))

# LOOIC

loo_nondamage <- loo(nondmg_fit)


## LOOIC comparison between damage model and nondamage model
## The preferred model will be at the first row. 
loo_compare(loo_dmg, loo_nondamage)

##----------------model with only R40----------------------------------
R40dmg_mod <- stan_model("only_alphaR40.stan")
init_R40dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R40 = 1)
}


R40dmg_fit <- sampling(object = R40dmg_mod,
                       data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                   N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                   N_x = length(T100_data),N_y = length(R100_data),
                                   X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                   X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                   t_x = R100_data,t_y = T100_data,
                                   l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                   l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                       control = list(adapt_delta = 0.8),init = init_R40dmg)
print(R40dmg_fit,pars = c('mu','sigma','rho','alpha_R40'))
pairs(R40dmg_fit,pars = c('mu','sigma','rho','alpha_R40'))
# LOOIC

loo_R40dmg <- loo(R40dmg_fit)
extract<- extract(R40dmg_fit)$'alpha_R40'
##
loo_compare(loo_R40dmg, loo_nondamage)






##-------------------Model 4 let alpha_R40 be real----------------
R40dmgreal_mod <- stan_model("only_alphaR40_real.stan")
init_R40dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R40 = 1)
}


R40dmgreal_fit <- sampling(object = R40dmgreal_mod,
                       data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                   N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                   N_x = length(T100_data),N_y = length(R100_data),
                                   X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                   X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                   t_x = R100_data,t_y = T100_data,
                                   l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                   l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                       control = list(adapt_delta = 0.8),init = init_R40dmg)
print(R40dmgreal_fit,pars = c('mu','sigma','rho','alpha_R40'))
pairs(R40dmgreal_fit,pars = c('mu','sigma','rho','alpha_R40'))


R40alpha <- extract(R40dmgreal_fit, pars = c('alpha_R40'))
moments::skewness(R40alpha[[1]])


loo_R40dmgrel <- loo(R40dmgreal_fit)
loo_compare(loo_R40dmg, loo_R40dmgrel)

##-------------------Model 5 let all alpha's be real----------------
dmgreal_mod <- stan_model("damage_real.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R20 = 1,
       alpha_R40 = 1,alpha_R60 = 1,alpha_T20 = 1,alpha_T40 =1,alpha_R60 = 1 )
}

dmgreal_fit <- sampling(object = dmgreal_mod,
                    data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                N_x = length(T100_data),N_y = length(R100_data),
                                X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                t_x = R100_data,t_y = T100_data,
                                l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                    control = list(adapt_delta = 0.8),init = init_dmg)

print(dmgreal_fit,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
                       'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))


pairs(dmgreal_fit,pars = c('rho','alpha_R20','alpha_R40',
                       'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))

# 4976
loo_dmgrel <- loo(dmgreal_fit)
##-------------model 6 real alpha--------------
R40dmg_mod <- stan_model("only_alphaR40_real.stan")
init_R40dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R40 = 1)
}
R40dmg_fit <- sampling(object = R40dmg_mod,
                       data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                   N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                   N_x = length(T100_data),N_y = length(R100_data),
                                   X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                   X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                   t_x = R100_data,t_y = T100_data,
                                   l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                   l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                       control = list(adapt_delta = 0.8),init = init_R40dmg)

 
##----------------model with only R40----------------------------------
R40T40dmg_mod <- stan_model("alphaR40T40.stan")
init_R40T40dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha_R40 = 1,alpha_T40 = 1)
}


R40T40dmg_fit <- sampling(object = R40T40dmg_mod,
                       data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                   N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                   N_x = length(T100_data),N_y = length(R100_data),
                                   X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                   X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                   t_x = R100_data,t_y = T100_data,
                                   l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                   l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                       control = list(adapt_delta = 0.8),init = init_R40T40dmg)
print(R40T40dmg_fit,pars = c('mu','sigma','rho','alpha_R40','alpha_T40'))
pairs(R40T40dmg_fit,pars = c('mu','sigma','rho','alpha_R40','alpha_T40'))
# LOOIC, 5057.9

loo_R40T40dmg_fit <- loo(R40T40dmg_fit)

##---------------Posterior predictive checks------------------------####


N = 87


extract_mu <- extract(nondmg_fit)$'mu'[3001:4000,]
extract_sigma <- extract(nondmg_fit)$'sigma'[3001:4000,]
extract_rho <- extract(nondmg_fit)$'rho'[3001:4000]


# R20

t_10_R20 <- rep(0,1000)
t_50_R20 <- rep(0,1000)
t_90_R20 <- rep(0,1000)


for(j in 1:1000){
  
  mu_x <- extract_mu[j,1]
  sd_x <- extract_sigma[j,1]
  mu_y <- extract_mu[j,2]
  sd_y <- extract_sigma[j,2]
  rho <- extract_rho[j]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp <- 0.2
  l <- quantile(samples[1:(pp*N),1],0.2)[[1]]
  PFY_ob_rep <- matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] <- samples[i,1]
      PFY_ob_rep[i,3] <- 1
    }
    else{
      PFY_ob_rep[i,2] <- samples[i,2]
      PFY_ob_rep[i,3] <- 0
    }
  }
  PFY_y_rep <- samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep <- which(PFY_ob_rep[,3] == 0)
  
  
  
  t_10_R20[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.1)
  t_90_R20[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.9)
  t_50_R20[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.5)
  
}


id_y <- which(R20_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10_R20)
abline(v = quantile(R20_data[id_y,2],0.1),col = 'red')
hist(t_50_R20)
abline(v = quantile(R20_data[id_y,2],0.5),col = 'red')
hist(t_90_R20)
abline(v = quantile(R20_data[id_y,2],0.9),col = 'red')



p_10_R20 <-  mean(quantile(R20_data[id_y,2],0.1)<t_10_R20)
p_50_R20 <- mean(quantile(R20_data[id_y,2],0.5)<t_50_R20)
p_90_R20 <-  mean(quantile(R20_data[id_y,2],0.9)<t_90_R20)



# R40

t_10_R40 <- rep(0,1000)
t_50_R40 <- rep(0,1000)
t_90_R40 <- rep(0,1000)

for(j in 1:1000){
  
  mu_x <- extract_mu[j,1]
  sd_x <- extract_sigma[j,1]
  mu_y <- extract_mu[j,2]
  sd_y <- extract_sigma[j,2]
  rho <- extract_rho[j]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp <- 0.4
  l <- quantile(samples[1:(pp*N),1],0.4)[[1]]
  PFY_ob_rep <- matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] <- samples[i,1]
      PFY_ob_rep[i,3] <- 1
    }
    else{
      PFY_ob_rep[i,2] <- samples[i,2]
      PFY_ob_rep[i,3] <- 0
    }
  }
  PFY_y_rep <- samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep <- which(PFY_ob_rep[,3] == 0)
  
  
  
  t_10_R40[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.1)
  t_90_R40[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.9)
  t_50_R40[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.5)
  
}



par(mfrow = c(2,2))
id_y <- which(R40_data[,3] == 0)
hist(t_10_R40)
abline(v = quantile(R40_data[id_y,2],0.1),col = 'red')
hist(t_50_R40)
abline(v = quantile(R40_data[id_y,2],0.5),col = 'red')
hist(t_90_R40)
abline(v = quantile(R40_data[id_y,2],0.9),col = 'red')



p_10_R40 <-  mean(quantile(R40_data[id_y,2],0.1)<t_10_R40)
p_50_R40 <- mean(quantile(R40_data[id_y,2],0.5)<t_50_R40)
p_90_R40 <-  mean(quantile(R40_data[id_y,2],0.9)<t_90_R40)






# R60

t_10_R60 <- rep(0,1000)
t_50_R60 <- rep(0,1000)
t_90_R60 <- rep(0,1000)

for(j in 1:1000){
  
  mu_x <- extract_mu[j,1]
  sd_x <- extract_sigma[j,1]
  mu_y <- extract_mu[j,2]
  sd_y <- extract_sigma[j,2]
  rho <- extract_rho[j]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp <- 0.6
  l <- quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep <- matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] <- samples[i,1]
      PFY_ob_rep[i,3] <- 1
    }
    else{
      PFY_ob_rep[i,2] <- samples[i,2]
      PFY_ob_rep[i,3] <- 0
    }
  }
  PFY_y_rep <- samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep <- which(PFY_ob_rep[,3] == 0)
  
  
  
  t_10_R60[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.1)
  t_90_R60[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.9)
  t_50_R60[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.5)
  
}


par(mfrow = c(2,2))

id_y <- which(R60_data[,3] == 0)
hist(t_10_R60)
abline(v = quantile(R60_data[id_y,2],0.1),col = 'red')
hist(t_50_R60)
abline(v = quantile(R60_data[id_y,2],0.5),col = 'red')
hist(t_90_R60)
abline(v = quantile(R60_data[id_y,2],0.9),col = 'red')



p_10_R60 <-  mean(quantile(R60_data[id_y,2],0.1)<t_10_R60)
p_50_R60 <- mean(quantile(R60_data[id_y,2],0.5)<t_50_R60)
p_90_R60 <-  mean(quantile(R60_data[id_y,2],0.9)<t_90_R60)




# T20

t_10_T20 <- rep(0,1000)
t_50_T20 <- rep(0,1000)
t_90_T20 <- rep(0,1000)

for(j in 1:1000){
  
  mu_x <- extract_mu[j,2]
  sd_x <- extract_sigma[j,2]
  mu_y <- extract_mu[j,1]
  sd_y <- extract_sigma[j,1]
  rho <- extract_rho[j]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp <- 0.2
  l <- quantile(samples[1:(pp*N),1],0.2)[[1]]
  PFY_ob_rep <- matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] <- samples[i,1]
      PFY_ob_rep[i,3] <- 1
    }
    else{
      PFY_ob_rep[i,2] <- samples[i,2]
      PFY_ob_rep[i,3] <- 0
    }
  }
  PFY_y_rep <- samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  
  
  
  t_10_T20[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.1)
  t_90_T20[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.9)
  t_50_T20[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.5)
  
}


id_y <- which(T20_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10_T20)
abline(v = quantile(T20_data[id_y,2],0.1),col = 'red')
hist(t_50_T20)
abline(v = quantile(T20_data[id_y,2],0.5),col = 'red')
hist(t_90_T20)
abline(v = quantile(T20_data[id_y,2],0.9),col = 'red')



p_10_T20 <-  mean(quantile(T20_data[id_y,2],0.1)<t_10_T20)
p_50_T20 <- mean(quantile(T20_data[id_y,2],0.5)<t_50_T20)
p_90_T20 <-  mean(quantile(T20_data[id_y,2],0.9)<t_90_T20)


# T40

t_10_T40 <- rep(0,1000)
t_50_T40 <- rep(0,1000)
t_90_T40 <- rep(0,1000)


for(j in 1:1000){
  
  mu_x <- extract_mu[j,2]
  sd_x <- extract_sigma[j,2]
  mu_y <- extract_mu[j,1]
  sd_y <- extract_sigma[j,1]
  rho <- extract_rho[j]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp <- 0.4
  l <- quantile(samples[1:(pp*N),1],0.4)[[1]]
  PFY_ob_rep <- matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] <- samples[i,1]
      PFY_ob_rep[i,3] <- 1
    }
    else{
      PFY_ob_rep[i,2] <- samples[i,2]
      PFY_ob_rep[i,3] <- 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  
  
  
  t_10_T40[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.1)
  t_90_T40[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.9)
  t_50_T40[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.5)
  
}


id_y <- which(T40_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10_T40)
abline(v = quantile(T40_data[id_y,2],0.1),col = 'red')
hist(t_50_T40)
abline(v = quantile(T40_data[id_y,2],0.5),col = 'red')
hist(t_90_T40)
abline(v = quantile(T40_data[id_y,2],0.9),col = 'red')



p_10_T40 <-  mean(quantile(T40_data[id_y,2],0.1)<t_10_T40)
p_50_T40 <- mean(quantile(T40_data[id_y,2],0.5)<t_50_T40)
p_90_T40 <-  mean(quantile(T40_data[id_y,2],0.9)<t_90_T40)



# T60

t_10_T60 <- rep(0,1000)
t_50_T60 <- rep(0,1000)
t_90_T60 <- rep(0,1000)


for(j in 1:1000){
  
  mu_x <- extract_mu[j,2]
  sd_x <- extract_sigma[j,2]
  mu_y <- extract_mu[j,1]
  sd_y <- extract_sigma[j,1]
  rho <- extract_rho[j]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp <- 0.6
  l <- quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep <- matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] <- samples[i,1]
      PFY_ob_rep[i,3] <- 1
    }
    else{
      PFY_ob_rep[i,2] <- samples[i,2]
      PFY_ob_rep[i,3] <- 0
    }
  }
  PFY_y_rep <- samples[(pp*N+1):N,2]
  
  # damage data
  id_y_rep <- which(PFY_ob_rep[,3] == 0)
  

  t_10_T60[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.1)
  t_90_T60[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.9)
  t_50_T60[j] <- quantile(PFY_ob_rep[id_y_rep,2],0.5)
  
}


id_y <- which(T60_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10_T60)
abline(v = quantile(T60_data[id_y,2],0.1),col = 'red')
hist(t_50_T60)
abline(v = quantile(T60_data[id_y,2],0.5),col = 'red')
hist(t_90_T60)
abline(v = quantile(T60_data[id_y,2],0.9),col = 'red')



p_10_T60 <-  mean(quantile(T60_data[id_y,2],0.1)<t_10_T60)
p_50_T60 <- mean(quantile(T60_data[id_y,2],0.5)<t_50_T60)
p_90_T60 <-  mean(quantile(T60_data[id_y,2],0.9)<t_90_T60)



##--------------plot the PPC for R group------------


par(mfrow = c(3,3))
id_y <- which(R20_data[,3] == 0)
hist(t_10_R20,main = "PPC for R20 (1)",xlab = "T(y) = quantile(y,0.1)")
abline(v = quantile(R20_data[id_y,2],0.1),col = 'red')
hist(t_50_R20,main = "PPC for R20 (2)",xlab = "T(y) = quantile(y,0.5)")
abline(v = quantile(R20_data[id_y,2],0.5),col = 'red')
hist(t_90_R20,main = "PPC for R20 (3)",xlab = "T(y) = quantile(y,0.9)")
abline(v = quantile(R20_data[id_y,2],0.9),col = 'red')

id_y <- which(R40_data[,3] == 0)
hist(t_10_R40,main = "PPC for R40 (1)",xlab = "T(y) = quantile(y,0.1)")
abline(v = quantile(R40_data[id_y,2],0.1),col = 'red')
hist(t_50_R40,main = "PPC for R40 (2)",xlab = "T(y) = quantile(y,0.5)")
abline(v = quantile(R40_data[id_y,2],0.5),col = 'red')
hist(t_90_R40,main = "PPC for R40 (3)",xlab = "T(y) = quantile(y,0.9)")
abline(v = quantile(R40_data[id_y,2],0.9),col = 'red')

id_y <- which(R60_data[,3] == 0)
hist(t_10_R60,main = "PPC for R60 (1)",xlab = "T(y) = quantile(y,0.1)")
abline(v = quantile(R60_data[id_y,2],0.1),col = 'red')
hist(t_50_R60,main = "PPC for R60 (2)",xlab = "T(y) = quantile(y,0.5)")
abline(v = quantile(R60_data[id_y,2],0.5),col = 'red')
hist(t_90_R60,main = "PPC for R60 (3)",xlab = "T(y) = quantile(y,0.9)")
abline(v = quantile(R60_data[id_y,2],0.9),col = 'red')



##---------------plot the ppc for T group--------------

par(mfrow = c(3,3))
id_y <- which(T20_data[,3] == 0)
hist(t_10_T20,main = "PPC for T20 (1)",xlab = "T(y) = quantile(y,0.1)")
abline(v = quantile(T20_data[id_y,2],0.1),col = 'red')
hist(t_50_T20,main = "PPC for T20 (2)",xlab = "T(y) = quantile(y,0.5)")
abline(v = quantile(T20_data[id_y,2],0.5),col = 'red')
hist(t_90_T20,main = "PPC for T20 (3)",xlab = "T(y) = quantile(y,0.9)")
abline(v = quantile(T20_data[id_y,2],0.9),col = 'red')

id_y <- which(T40_data[,3] == 0)
hist(t_10_T40,main = "PPC for T40 (1)",xlab = "T(y) = quantile(y,0.1)")
abline(v = quantile(T40_data[id_y,2],0.1),col = 'red')
hist(t_50_T40,main = "PPC for T40 (2)",xlab = "T(y) = quantile(y,0.5)")
abline(v = quantile(T40_data[id_y,2],0.5),col = 'red')
hist(t_90_T40,main = "PPC for T40 (3)",xlab = "T(y) = quantile(y,0.9)")
abline(v = quantile(T40_data[id_y,2],0.9),col = 'red')

id_y <- which(T60_data[,3] == 0)
hist(t_10_T60,main = "PPC for T60 (1)",xlab = "T(y) = quantile(y,0.1)")
abline(v = quantile(T60_data[id_y,2],0.1),col = 'red')
hist(t_50_T60,main = "PPC for T60 (2)",xlab = "T(y) = quantile(y,0.5)")
abline(v = quantile(T60_data[id_y,2],0.5),col = 'red')
hist(t_90_T60,main = "PPC for T60 (3)",xlab = "T(y) = quantile(y,0.9)")
abline(v = quantile(T60_data[id_y,2],0.9),col = 'red')

