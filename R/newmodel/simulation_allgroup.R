setwd("~/Downloads/github/Bayesian-lumber-strength/R/newmodel")
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

# prop is thresh 
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



# 
# int_function_R <- function(ystar, mu,sigma, rho,l){
#   
#   a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
#   
#   dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)
#   
#   
# }
# 
# 
# PFY_lik_R <- function(mu, sigma, rho,alpha,l,data){
#   
#   (sum(dnorm(data[data[,3] == 1,1],
#              mean = mu[1], 
#              sd = sigma[1], 
#              log = T)
#   ) + 
#     sum(log(1/alpha*(int_function_R(data[data[,3] == 0,2]/alpha,mu,sigma,rho,l)-
#                        int_function_R(data[data[,3] == 0,2]/alpha,mu,sigma,rho,1/prop*l)) +  int_function_R(data[data[,3] == 0,2],mu,sigma,rho,1/prop*l)
#     ))
#   
#   
#   )
#   
#   
#   
# }
# 
# 
# 
# int_function_T <- function(ystar, mu,sigma, rho,l){
#   
#   a_l <- (l - mu[2]-rho*(sigma[2]/sigma[1])*(ystar-mu[1]))/(sigma[2]*sqrt(1-rho^2));
#   
#   dnorm(ystar, mean = mu[1], sd = sigma[1],log = F)*pnorm(a_l,lower.tail = F, log.p = F)
#   
#   
# }
# 
# 
# PFY_lik_T <- function(mu, sigma, rho,alpha,l,data){
#   
#   (sum(dnorm(data[data[,3] == 1,1],
#              mean = mu[2], 
#              sd = sigma[2], 
#              log = T)
#   ) + 
#     sum(log(1/alpha*(int_function_T(data[data[,3] == 0,2]/alpha,mu,sigma,rho,l)-
#                        int_function_T(data[data[,3] == 0,2]/alpha,mu,sigma,rho,1/prop*l)) +  int_function_T(data[data[,3] == 0,2],mu,sigma,rho,1/prop*l)
#     ))
#   
#   
#   )
#   
#   
#   
# }
# # sum(dnorm(data[data[,3] == 1,1], mean = mu[1], sd = sigma[1], log = T))
# 
# nlogpost <- function(theta){
#   mu <- theta[1:2]
#   sigma <- theta[3:4]
#   rho <- theta[5]
#   alpha <- theta[6:11]
#   lik <- PFY_lik_R(mu,sigma,rho,alpha[1],R_pf[1],R20_data) + 
#     PFY_lik_R(mu,sigma,rho,alpha[2],R_pf[2],R40_data) +
#     PFY_lik_R(mu,sigma,rho,alpha[3],R_pf[3],R60_data) + 
#     PFY_lik_T(mu,sigma,rho,alpha[4],T_pf[1],T20_data) +
#     PFY_lik_T(mu,sigma,rho,alpha[5],T_pf[2],T40_data) +
#     PFY_lik_T(mu,sigma,rho,alpha[6],T_pf[3],T60_data) +
#     sum(dnorm(T100_data,
#           mean = mu[2],
#           sd = sigma[2], log = T))+
#     sum(dnorm(R100_data,
#               mean = mu[1],
#               sd = sigma[1], log = T))
#   return(-1*lik)
# }

# nlogpost(c(mu,sigma,rho,alpha))




# optimresult <- optim(c(mu,sigma,rho,alpha),nlogpost,method = 'L-BFGS-B',
#                    lower = rep(0,11), upper = c(rep(Inf,4), 1,rep(Inf,6)))
# optimresult <- optim(c(mu,sigma,rho,alpha[1:3]),nlogpost)
# optimresult$par


# optim_proj(fun = nlogpost, # objective function
#            xsol = optimresult$par, # candidate optimum
#            maximize = TRUE, # look for max
#            xnames = c("mu_x","mu_y","sig_x","sig_y","rho","alpha")) # pretty names on plot




options(mc.cores = parallel::detectCores())
dmg_mod <- stan_model("newmodel.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5,
       alpha = rep(1,6) )
}
dmg_fit <- sampling(object = dmg_mod,
                    data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                N_x = length(T100_data),N_y = length(R100_data),
                                X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                t_x = R100_data,t_y = T100_data,
                                l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3] ,thresh = .7),
                    init = init_dmg)


dmg_fit <- readRDS("dmgfit.rds")
print(dmg_fit)

traceplot(dmg_fit, pars = 'alpha')
traceplot(dmg_fit, pars = c('mu','sigma','rho'))







options(mc.cores = parallel::detectCores())
dmg_mod <- stan_model("newmodel_test.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5,
       alpha = rep(1,6),thresh = 1,
       eta = 1, mu_eta = 0, eta_eta2 = .25)
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
# saveRDS(dmg_fit, file = 'dmgfit_ig.rds')
#dmg_fit <- readRDS("dmg_fit_thresh.rds")
print(dmg_fit)

result <- round(c(summary(dmg_fit)[[1]][16,1],
                  summary(dmg_fit)[[1]][15,c(1,4,6,8)]),3)
saveRDS(result, file = 'thresh_ig1.rds')

# prop .75,  -8462.49 
# prop .9,   -8463.54 
# prop .5,   -8463.92


# dmg_fit <- readRDS("dmg_fit.rds")

# 1.
# In the simulation, how the different thresholds affect the results?
# if alpha = rep(1,6), different thresholds won't affect the result (estimation and lp__)
# if alpha is not close to 1, inaccurate thresholds can cause bias estimation of alpha 
# if different thresholds cannot change the lp__, it is a good evidence that there is no damage  
# but I have not yet tested how close alpha to 1 it can happen 
# should different groups have different thresholds?



# 2.
# Do alpha's go larger or smaller when it moves from T20 to T60, try alpha = 0.9 or 0.8, 0.95
# tried alpha <- c(.8,.95,.85,.8,.95,.85)  good 
# also alpha <- c(.8,.85,.9,.8,.85,.9) 

# 3.
# may need some explanations why the threshold cannot affect the results in the real data (may alpha is too close to 1, can use simulation to check)
# look at the model 
# Y_i* = Y_i * alpha * I(l > x_i/2) +  Y_i * I(l < x_i/2),

# 4. 
# supposing we know the threshold, and check whether the estimation is unbiased
# Yes, it is unbiased. 


# 5.
# if we don't know the threshold, and how to proceed ? 
# (may pick the highest lp__ that the threshold provides)
# make a plot lp__ vs prop on the same dataset 
# when alpha <- c(.8,.85,.9,.8,.85,.9), true prop = 0.75



# 6. move the outliers to 100 groups 


# prop 0.5, -8763.81
# prop 0.55,  -8761.56   
# prop 0.6, -8753.34
# prop .65,  -8722.71
# prop .7,  -8689.28
# prop .75,  -8681.00
# prop .8,  -8699.05


plot(seq(.5,.8,.05), c(-8763.81, -8761.56,-8753.34, -8722.71, -8689.28, -8681.00,-8699.05),type = 'l', xlab = 'threshold', ylab = 'lp')



# real data
# prop .2, -2481.83
# prop .3,  -2481.82
# prop .4, -2481.86
# prop .5,  -2482.00
# prop .6, -2482.05
# prop .7, -2482.00
# prop .8, -2481.98,-2483.25 
# prop .85,-2481.52  , -2482.97 
# prop .9, -2481.07 ,  -2482.42
# prop .95, ,         -2483.12

# 6. move the outliers to 100 groups 
# 7. to find a better threshold between .85, .9
