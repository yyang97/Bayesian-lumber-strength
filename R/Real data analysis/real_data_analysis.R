setwd("F:/study/Research/Bayesian-lumber-strength/Data")

require(rstan)
require(MASS)
install.packages("jsonlite", type = "source")
devtools::find_rtools()
#load("real_data_analysis.RData")

##-------------convert log to sqrt-------------###
# load("all_data_one_model.RData")
# id <- which(R20_data[,3] == 0)
# R20_data[id,2] <- sqrt(exp(R20_data[id,2]))
# 
# id <- which(R40_data[,3] == 0)
# R40_data[id,2] <- sqrt(exp(R40_data[id,2]))
# 
# id <- which(R60_data[,3] == 0)
# R60_data[id,2] <- sqrt(exp(R60_data[id,2]))
# 
# 
# id <- which(T20_data[,3] == 1)
# T20_data[id,1] <- sqrt(exp(T20_data[id,1]))
# 
# 
# id <- which(T40_data[,3] == 1)
# T40_data[id,1] <- sqrt(exp(T40_data[id,1]))
# 
# 
# id <- which(T60_data[,3] == 1)
# T60_data[id,1] <- sqrt(exp(T60_data[id,1]))
# 
# T100_data <- sqrt(exp(T100_data))

###------------------data preprocessing--------------
library(readxl)
library(dplyr)
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
T60_data$MOR[id] = 0
T60_data$UTS[id] = as.numeric(T60_data$UTS[id])
T60_data$UTS[!id] = 0
T60_data$UTS <- as.numeric(T60_data$UTS)
T60_data$MOR <- as.numeric(T60_data$MOR)


id <- T40_data$Broken == 1
T40_data$MOR[id] = 0
T40_data$UTS[id] = as.numeric(T40_data$UTS[id])
T40_data$UTS[!id] = 0
T40_data$UTS <- as.numeric(T40_data$UTS)
T40_data$MOR <- as.numeric(T40_data$MOR)


id <- T20_data$Broken == 1
T20_data$MOR[id] = 0
T20_data$UTS[id] = as.numeric(T20_data$UTS[id])
T20_data$UTS[!id] = 0
T20_data$UTS <- as.numeric(T20_data$UTS)
T20_data$MOR <- as.numeric(T20_data$MOR)



##-----------------R substitute NA to 0 ------

id <- R60_data$Broken == 0
R60_data$MOR[id] = 0
R60_data$UTS[id] = as.numeric(R60_data$UTS[id])
R60_data$UTS[!id] = 0
R60_data$UTS <- as.numeric(R60_data$UTS)
R60_data$MOR <- as.numeric(R60_data$MOR)



id <- R40_data$Broken == 0
R40_data$MOR[id] = 0
R40_data$UTS[id] = as.numeric(R40_data$UTS[id])
R40_data$UTS[!id] = 0
R40_data$UTS <- as.numeric(R40_data$UTS)
R40_data$MOR <- as.numeric(R40_data$MOR)

id <- R20_data$Broken == 0
R20_data$MOR[id] = 0
R20_data$UTS[id] = as.numeric(R20_data$UTS[id])
R20_data$UTS[!id] = 0
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

library(fitdistrplus)

# check R100_data
FIT <- fitdist(R100_data, "norm")    ## note: it is "norm" not "normal"
plot(FIT)    ## use method `plot.fitdist`
FIT$estimate # mean = 45.679 sd = 12.900
FIT$bic # bic = 1394.022
# good normal fitting 


# check T100_data

shapiro.test(T100_data) # p = 0.0001841
shapiro.test(log(T100_data)) # p-value = 0.002327
shapiro.test(sqrt(T100_data)) # p-value = 0.3988



##------proof loading-------------- 
R_pf =c* c(4.956690733, 6.110714122, 7.092435407)
T_pf =sqrt(c* c(2.962390379, 3.986497991, 4.916102264))


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

##------ Stan for all alpha-------------
setwd("F:/study/Research/Bayesian-lumber-strength/R/Real data analysis")
dmg_mod_waic <- stan_model("test.stan")

dmg_fit_waic <- sampling(object = dmg_mod_waic,
                         data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                     N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                     N_x = length(T100_data),N_y = length(R100_data),
                                     X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                     X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                     t_x = R100_data,t_y = T100_data,
                                     l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                     l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                         control = list(adapt_delta = 0.8))

print(dmg_fit_waic,pars = c('mu','sigma','rho','alpha_R20','alpha_R40',
                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))


pairs(dmg_fit_waic,pars = c('rho','alpha_R20','alpha_R40',
                            'alpha_R60','alpha_T20','alpha_T40','alpha_T60'))

pairs((extract(dmg_fit_waic)))

##-------without alpha---------






##----------------- Stan for R20------------

R20_data <- cbind(R20_data$MOR,R20_data$UTS,R20_data$Broken)
R20_data <- R20_data[-48,]
#R20_data[,1] < R_pf[1]


dmg_mod <- stan_model("damage.stan")
dmg_fit_R20 <- sampling(object = dmg_mod,
                    data = list(N_SPLD = nrow(R20_data),N_x = length(R100_data),
                                N_y = length(T100_data),
                                X = R20_data,t_x = R100_data,t_y = T100_data,
                                l=R_pf[1]),
                    control = list(adapt_delta = 0.8))

print(dmg_fit_R20)
theta <- c('mu','sigma','rho','alpha')
pairs(dmg_fit_20,pars = theta)



cor(extract(dmg_fit_R20)$'rho',extract(dmg_fit_R20)$'alpha')


for(j in 1:1000){
  N = 87
  mu_x <- extract(dmg_fit_R20)$'mu'[j+3000,1]
  sd_x <- extract(dmg_fit_R20)$'sigma'[j+3000,1]
  mu_y <- extract(dmg_fit_R20)$'mu'[j+3000,2]
  sd_y <- extract(dmg_fit_R20)$'sigma'[j+3000,2]
  rho <- extract(dmg_fit_R20)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R20)$'alpha'[j+3000]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  quantile(PFY_dmg_rep[id_y_rep,2],.5)
  
  t_min[j] <- min(PFY_dmg_rep[id_y_rep,2])
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_max[j] <- max(PFY_dmg_rep[id_y_rep,2])
  t_sd[j] <- sd(PFY_dmg_rep[id_y_rep,2])
  
}

id_y <- which(R20_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_min)
abline(v = min(R20_data[id_y,2]),col = 'red')
hist(t_mean)
abline(v = mean(R20_data[id_y,2]),col = 'red')
hist(t_max)
abline(v = max(R20_data[id_y,2]),col = 'red')
hist(t_sd)
abline(v = sd(R20_data[id_y,2]),col = 'red')




##----------------- Stan for R40------------

R40_data <- cbind(R40_data$MOR,R40_data$UTS,R40_data$Broken)
R20_data <- R20_data[-48,]
#R40_data[,1] < R_pf[2]


dmg_mod <- stan_model("damage.stan")
dmg_fit_R40 <- sampling(object = dmg_mod,
                    data = list(N_SPLD = nrow(R40_data),N_x = length(R100_data),
                                N_y = length(T100_data),
                                X = R40_data,t_x = R100_data,t_y = T100_data,
                                l=R_pf[2]),
                    control = list(adapt_delta = 0.8))

print(dmg_fit_R40)
theta <- c('mu','sigma','rho','alpha')
pairs(dmg_fit_R40,pars = theta)
pairs(dmg_fit_R40,pars = c('rho','alpha'))




for(j in 1:1000){
  N = 87
  mu_x <- extract(dmg_fit_R40)$'mu'[j+3000,1]
  sd_x <- extract(dmg_fit_R40)$'sigma'[j+3000,1]
  mu_y <- extract(dmg_fit_R40)$'mu'[j+3000,2]
  sd_y <- extract(dmg_fit_R40)$'sigma'[j+3000,2]
  rho <- extract(dmg_fit_R40)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R40)$'alpha'[j+3000]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  quantile(PFY_dmg_rep[id_y_rep,2],.5)
  
  t_min[j] <- min(PFY_dmg_rep[id_y_rep,2])
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_max[j] <- max(PFY_dmg_rep[id_y_rep,2])
  t_sd[j] <- sd(PFY_dmg_rep[id_y_rep,2])
  
}

id_y <- which(R40_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_min)
abline(v = min(R40_data[id_y,2]),col = 'red')
hist(t_mean)
abline(v = mean(R40_data[id_y,2]),col = 'red')
hist(t_max)
abline(v = max(R40_data[id_y,2]),col = 'red')
hist(t_sd)
abline(v = sd(R40_data[id_y,2]),col = 'red')
##----------------- Stan for R40------------

R60_data <- cbind(R60_data$MOR,R60_data$UTS,R60_data$Broken)
R60_data <- R60_data[-which((R60_data[,1] < R_pf[3]) == F),]



dmg_mod <- stan_model("damage.stan")
dmg_fit_R60 <- sampling(object = dmg_mod,
                    data = list(N_SPLD = nrow(R60_data),N_x = length(R100_data),
                                N_y = length(T100_data),
                                X = R60_data,t_x = R100_data,t_y = T100_data,
                                l=R_pf[3]),
                    control = list(adapt_delta = 0.8))

print(dmg_fit_R60)
theta <- c('mu','sigma','rho','alpha')
pairs(dmg_fit_R60,pars = theta)

pairs(dmg_fit_R60,pars = c('rho','alpha'))

for(j in 1:1000){
  N = 87
  mu_x <- extract(dmg_fit_R60)$'mu'[j+3000,1]
  sd_x <- extract(dmg_fit_R60)$'sigma'[j+3000,1]
  mu_y <- extract(dmg_fit_R60)$'mu'[j+3000,2]
  sd_y <- extract(dmg_fit_R60)$'sigma'[j+3000,2]
  rho <- extract(dmg_fit_R60)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R60)$'alpha'[j+3000]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  quantile(PFY_dmg_rep[id_y_rep,2],.5)
  
  t_min[j] <- min(PFY_dmg_rep[id_y_rep,2])
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_max[j] <- max(PFY_dmg_rep[id_y_rep,2])
  t_sd[j] <- sd(PFY_dmg_rep[id_y_rep,2])
  
}

id_y <- which(R60_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_min)
abline(v = min(R60_data[id_y,2]),col = 'red')
hist(t_mean)
abline(v = mean(R60_data[id_y,2]),col = 'red')
hist(t_max)
abline(v = max(R60_data[id_y,2]),col = 'red')
hist(t_sd)
abline(v = sd(R60_data[id_y,2]),col = 'red')
##--------------------Stan T20---------------
# T20
T20_data <- cbind(T20_data$UTS,T20_data$MOR,T20_data$Broken)
T20_data <- T20_data[-74,] 
#T20_data[,1] < sqrt(T_pf[1])


wood_mod <- stan_model("wood.stan")
wood_fit <- sampling(object = wood_mod,
                     data = list(N_SPLD = nrow(T20_data),N_x = length(T100_data),
                                 N_y = length(R100_data),
                                 X = T20_data,t_x = T100_data,t_y = R100_data,
                                 l=T_pf[1]),
                     control = list(adapt_delta = 0.8))

print(wood_fit)
pairs(wood_fit)




for(j in 1:1000){
  N = 87
  mu_x <- extract(wood_fit)$'mu'[j+3000,1]
  sd_x <- extract(wood_fit)$'sigma'[j+3000,1]
  mu_y <- extract(wood_fit)$'mu'[j+3000,2]
  sd_y <- extract(wood_fit)$'sigma'[j+3000,2]
  rho <- extract(wood_fit)$'rho'[j+3000]
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)

  t_min[j] <- min(PFY_ob_rep[id_y_rep,2])
  t_mean[j] <- mean(PFY_ob_rep[id_y_rep,2])
  t_max[j] <- max(PFY_ob_rep[id_y_rep,2])
  t_sd[j] <- sd(PFY_ob_rep[id_y_rep,2])
  
}

id_y <- which(T20_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_min)
abline(v = min(T20_data[id_y,2]),col = 'red')
hist(t_mean)
abline(v = mean(T20_data[id_y,2]),col = 'red')
hist(t_max)
abline(v = max(T20_data[id_y,2]),col = 'red')
hist(t_sd)
abline(v = sd(T20_data[id_y,2]),col = 'red')




dmg_fit_T20 <- sampling(object = dmg_mod,
                     data = list(N_SPLD = nrow(T20_data),N_x = length(T100_data),
                                 N_y = length(R100_data),
                                 X = T20_data,t_x = T100_data,t_y = R100_data,
                                 l=T_pf[1]),
                     control = list(adapt_delta = 0.99))

print(dmg_fit_T20)


# T 40
T40_data <- cbind(T40_data$UTS,T40_data$MOR,T40_data$Broken)
T40_data <- T40_data[-81,] 
#T40_data[,1] < (T_pf[2])

dmg_fit_T40 <- sampling(object = dmg_mod,
                     data = list(N_SPLD = nrow(T40_data),N_x = length(T100_data),
                                 N_y = length(R100_data),
                                 X = T40_data,t_x = T100_data,t_y = R100_data,
                                 l=T_pf[2]),
                     control = list(adapt_delta = 0.99))

print(dmg_fit_T40)




# T 60
T60_data <- cbind(T60_data$UTS,T60_data$MOR,T60_data$Broken)
T60_data <- T60_data[-(20:21),] 
# T60_data[,1] < T_pf[3]

dmg_fit_T60 <- sampling(object = dmg_mod,
                        data = list(N_SPLD = nrow(T60_data),N_x = length(T100_data),
                                    N_y = length(R100_data),
                                    X = T60_data,t_x = T100_data,t_y = R100_data,
                                    l=T_pf[3]),
                        control = list(adapt_delta = 0.99))

print(dmg_fit_T60)





##--------------------R20 R40 R60--------------------############

dmg_mod_R <- stan_model("damage_combined.stan")

dmg_fit_R <- sampling(object = dmg_mod_R,
                        data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                    N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                    N_x = length(T100_data),N_y = length(R100_data),
                                    X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                    X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                    t_x = R100_data,t_y = T100_data,
                                    l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                    l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                        control = list(adapt_delta = 0.8))
print(dmg_fit_R)



alpha = c('alpha_R20','alpha_R40','alpha_R60','alpha_T20','alpha_T40','alpha_T60')
traceplot(dmg_fit_R,pars = alpha)
pairs(dmg_fit_R,pars =c ('rho',alpha))






##-----------------no constraint--------######


load("all_data_one_model.RData")

print(dmg_fit_R)


# R20

t_10 <- rep(1000,0)
t_mean <- rep(1000,0)
t_90 <- rep(1000,0)
t_50 <- rep(1000,0)

N = 87

for(j in 1:1000){
  
  mu_x <- extract(dmg_fit_R)$'mu'[j+3000,1]
  sd_x <- extract(dmg_fit_R)$'sigma'[j+3000,1]
  mu_y <- extract(dmg_fit_R)$'mu'[j+3000,2]
  sd_y <- extract(dmg_fit_R)$'sigma'[j+3000,2]
  rho <- extract(dmg_fit_R)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R)$'alpha_R20'[j+3000]
  # mu_x = 6.48
  # sd_x = 1.85
  # mu_y = 1.43
  # sd_y = 0.4
  # rho = 0.2
  # alpha = 0.3
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  
  t_10[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.1)
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_90[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.9)
  t_50[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.5)
  
}


id_y <- which(R20_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10)
abline(v = quantile(R20_data[id_y,2],0.1),col = 'red')
hist(t_mean)
abline(v = mean(R20_data[id_y,2]),col = 'red')
hist(t_90)
abline(v = quantile(R20_data[id_y,2],0.9),col = 'red')
hist(t_50)
abline(v = quantile(R20_data[id_y,2],0.5),col = 'red')


p_10 <-  mean(quantile(R20_data[id_y,2],0.1)<t_10)
p_mean <-  mean(mean(R20_data[id_y,2])<t_mean)
p_90 <-  mean(quantile(R20_data[id_y,2],0.9)<t_90)
p_50 <- mean(quantile(R20_data[id_y,2],0.5)<t_50)

p_10
p_mean
p_90
p_50

##----R40------

for(j in 1:1000){
  
  mu_x <- extract(dmg_fit_R)$'mu'[j+3000,1]
  sd_x <- extract(dmg_fit_R)$'sigma'[j+3000,1]
  mu_y <- extract(dmg_fit_R)$'mu'[j+3000,2]
  sd_y <- extract(dmg_fit_R)$'sigma'[j+3000,2]
  rho <- extract(dmg_fit_R)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R)$'alpha_R40'[j+3000]
  # mu_x = 6.48
  # sd_x = 1.85
  # mu_y = 1.43
  # sd_y = 0.4
  # rho = 0.2
  # alpha = 0.3
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  
  t_10[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.1)
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_90[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.9)
  t_50[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.5)
  
}


id_y <- which(R40_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10)
abline(v = quantile(R40_data[id_y,2],0.1),col = 'red')
hist(t_mean)
abline(v = mean(R40_data[id_y,2]),col = 'red')
hist(t_90)
abline(v = quantile(R40_data[id_y,2],0.9),col = 'red')
hist(t_50)
abline(v = quantile(R40_data[id_y,2],0.5),col = 'red')


p_10 <-  mean(quantile(R40_data[id_y,2],0.1)<t_10)
p_mean <-  mean(mean(R40_data[id_y,2])<t_mean)
p_90 <-  mean(quantile(R40_data[id_y,2],0.9)<t_90)
p_50 <- mean(quantile(R40_data[id_y,2],0.5)<t_50)

p_10
p_mean
p_90
p_50

#######-------R60----------######

for(j in 1:1000){
  
  mu_x <- extract(dmg_fit_R)$'mu'[j+3000,1]
  sd_x <- extract(dmg_fit_R)$'sigma'[j+3000,1]
  mu_y <- extract(dmg_fit_R)$'mu'[j+3000,2]
  sd_y <- extract(dmg_fit_R)$'sigma'[j+3000,2]
  rho <- extract(dmg_fit_R)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R)$'alpha_R60'[j+3000]
  # mu_x = 6.48
  # sd_x = 1.85
  # mu_y = 1.43
  # sd_y = 0.4
  # rho = 0.2
  # alpha = 0.3
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  
  t_10[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.1)
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_90[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.9)
  t_50[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.5)
  
}


id_y <- which(R60_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10)
abline(v = quantile(R60_data[id_y,2],0.1),col = 'red')
hist(t_mean)
abline(v = mean(R60_data[id_y,2]),col = 'red')
hist(t_90)
abline(v = quantile(R60_data[id_y,2],0.9),col = 'red')
hist(t_50)
abline(v = quantile(R60_data[id_y,2],0.5),col = 'red')


p_10 <-  mean(quantile(R60_data[id_y,2],0.1)<t_10)
p_mean <-  mean(mean(R60_data[id_y,2])<t_mean)
p_90 <-  mean(quantile(R60_data[id_y,2],0.9)<t_90)
p_50 <- mean(quantile(R60_data[id_y,2],0.5)<t_50)


p_10
p_mean
p_90
p_50
#####-------T20------######


for(j in 1:1000){
  
  mu_x <- extract(dmg_fit_R)$'mu'[j+3000,2]
  sd_x <- extract(dmg_fit_R)$'sigma'[j+3000,2]
  mu_y <- extract(dmg_fit_R)$'mu'[j+3000,1]
  sd_y <- extract(dmg_fit_R)$'sigma'[j+3000,1]
  rho <- extract(dmg_fit_R)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R)$'alpha_T20'[j+3000]
  # mu_x = 6.48
  # sd_x = 1.85
  # mu_y = 1.43
  # sd_y = 0.4
  # rho = 0.2
  # alpha = 0.3
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  
  t_10[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.1)
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_90[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.9)
  t_50[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.5)
  
}


id_y <- which(T20_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10)
abline(v = quantile(T20_data[id_y,2],0.1),col = 'red')
hist(t_mean)
abline(v = mean(T20_data[id_y,2]),col = 'red')
hist(t_90)
abline(v = quantile(T20_data[id_y,2],0.9),col = 'red')
hist(t_50)
abline(v = quantile(T20_data[id_y,2],0.5),col = 'red')


p_10 <-  mean(quantile(T20_data[id_y,2],0.1)<t_10)
p_mean <-  mean(mean(T20_data[id_y,2])<t_mean)
p_90 <-  mean(quantile(T20_data[id_y,2],0.9)<t_90)
p_50 <- mean(quantile(T20_data[id_y,2],0.5)<t_50)


p_10
p_mean
p_90
p_50

##----------T40-------######




for(j in 1:1000){
  
  mu_x <- extract(dmg_fit_R)$'mu'[j+3000,2]
  sd_x <- extract(dmg_fit_R)$'sigma'[j+3000,2]
  mu_y <- extract(dmg_fit_R)$'mu'[j+3000,1]
  sd_y <- extract(dmg_fit_R)$'sigma'[j+3000,1]
  rho <- extract(dmg_fit_R)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R)$'alpha_T40'[j+3000]
  # mu_x = 6.48
  # sd_x = 1.85
  # mu_y = 1.43
  # sd_y = 0.4
  # rho = 0.2
  # alpha = 0.3
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  
  t_10[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.1)
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_90[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.9)
  t_50[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.5)
  
}


id_y <- which(T40_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10)
abline(v = quantile(T40_data[id_y,2],0.1),col = 'red')
hist(t_mean)
abline(v = mean(T40_data[id_y,2]),col = 'red')
hist(t_90)
abline(v = quantile(T40_data[id_y,2],0.9),col = 'red')
hist(t_50)
abline(v = quantile(T40_data[id_y,2],0.5),col = 'red')


p_10 <-  mean(quantile(T40_data[id_y,2],0.1)<t_10)
p_mean <-  mean(mean(T40_data[id_y,2])<t_mean)
p_90 <-  mean(quantile(T40_data[id_y,2],0.9)<t_90)
p_50 <- mean(quantile(T40_data[id_y,2],0.5)<t_50)


p_10
p_mean
p_90
p_50

###-------T60---------######



for(j in 1:1000){
  
  mu_x <- extract(dmg_fit_R)$'mu'[j+3000,2]
  sd_x <- extract(dmg_fit_R)$'sigma'[j+3000,2]
  mu_y <- extract(dmg_fit_R)$'mu'[j+3000,1]
  sd_y <- extract(dmg_fit_R)$'sigma'[j+3000,1]
  rho <- extract(dmg_fit_R)$'rho'[j+3000]
  alpha <- extract(dmg_fit_R)$'alpha_T60'[j+3000]
  # mu_x = 6.48
  # sd_x = 1.85
  # mu_y = 1.43
  # sd_y = 0.4
  # rho = 0.2
  # alpha = 0.3
  mu <- c(mu_x,mu_y)
  sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
  bvn1 <- mvrnorm(N, mu = mu, Sigma = sigma ) # from MASS package
  colnames(bvn1) <- c("bvn1_X","bvn1_Y")
  samples <- bvn1
  
  pp = 0.6
  l = quantile(samples[1:(pp*N),1],0.6)[[1]]
  PFY_ob_rep = matrix(0, nrow = pp*N, ncol = 3)
  for (i in 1:nrow(PFY_ob_rep)){
    if(samples[i,1] <= l){
      PFY_ob_rep[i,1] = samples[i,1]
      PFY_ob_rep[i,3] = 1
    }
    else{
      PFY_ob_rep[i,2] = samples[i,2]
      PFY_ob_rep[i,3] = 0
    }
  }
  PFY_y_rep = samples[(pp*N+1):N,2]
  
  # damage data
  
  id_y_rep = which(PFY_ob_rep[,3] == 0)
  PFY_dmg_rep = PFY_ob_rep
  PFY_dmg_rep[id_y_rep,2] = PFY_ob_rep[id_y_rep,2] - alpha/PFY_ob_rep[id_y_rep,2]
  
  remove(PFY_ob_rep)
  
  
  t_10[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.1)
  t_mean[j] <- mean(PFY_dmg_rep[id_y_rep,2])
  t_90[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.9)
  t_50[j] <- quantile(PFY_dmg_rep[id_y_rep,2],0.5)
  
}


id_y <- which(T60_data[,3] == 0)
par(mfrow = c(2,2))
hist(t_10)
abline(v = quantile(T60_data[id_y,2],0.1),col = 'red')
hist(t_mean)
abline(v = mean(T60_data[id_y,2]),col = 'red')
hist(t_90)
abline(v = quantile(T60_data[id_y,2],0.9),col = 'red')
hist(t_50)
abline(v = quantile(T60_data[id_y,2],0.5),col = 'red')


p_10 <-  mean(quantile(T60_data[id_y,2],0.1)<t_10)
p_mean <-  mean(mean(T60_data[id_y,2])<t_mean)
p_90 <-  mean(quantile(T60_data[id_y,2],0.9)<t_90)
p_50 <- mean(quantile(T60_data[id_y,2],0.5)<t_50)


p_10
p_mean
p_90
p_50





##-----------------normal testing for UTS data-------------

T100_data <- exp(T100_data)



fit <- fitdistr(T100_data, "normal")
para <- fit$estimate
hist(T100_data, prob = TRUE)
curve(dnorm(T100_data, para[1], para[2]), col = 2, add = TRUE)



library(fitdistrplus)

FIT <- fitdist(T100_data, "norm")    ## note: it is "norm" not "normal"
class(FIT)
# [1] "fitdist"

plot(FIT)    ## use method `plot.fitdist`
FIT$estimate
FIT$bic
