##-----------set environment and check--------------------
#setwd("~/Downloads/github/Bayesian-lumber-strength/R/Simulation")
require(rstan)
require(MASS)
require(loo)
options(mc.cores = parallel::detectCores())

load('allSTAN.rda')



nsim <- 10
loo_value <- vector("list",nsim)
rho_value <- vector("list",nsim)
for (jj in 1:nsim){
  loo_fit <- vector("list",64)
  show(jj)
  

    ##-----mu sigma-----####
  mu <- c(44.73,5.49)
  sigma <- c(12.74,1.05)
  
  rho <- 0.8
  
  ##------proof loading-----####
  R_pf <- c(qnorm(0.2,mu[1],sigma[1]), qnorm(0.4,mu[1],sigma[1]),
            qnorm(0.6,mu[1],sigma[1]))
  T_pf <- c(qnorm(0.2,mu[2],sigma[2]),qnorm(0.4,mu[2],sigma[2]),
            qnorm(0.6,mu[2],sigma[2]))
  N <- 87*2
  
  
  ##--------------alpha-----#####
  
 alpha_R <- 2* c(0,1.83,0)
 alpha_T <- 2*c(0,104.09,65.95)

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
  
  allcombs <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
  rho_est <- matrix(NA,nrow = 64, ncol = 3)
  for (ii in 1:64) {
    show(paste('model', ii))
     
     ci <- allcombs[ii,]
     # initial value 
     init_value <- function(){
       init <- list (mu = c(35,8), sigma = c(10,1), rho = .5)
       if (ci[1] == 1){init$alpha_R20 <- 1}
       if (ci[2] == 1){init$alpha_R40 <- 1}
       if (ci[3] == 1){init$alpha_R60 <- 1}
       if (ci[4] == 1){init$alpha_T20 <- 1}
       if (ci[5] == 1){init$alpha_T40 <- 1}
       if (ci[6] == 1){init$alpha_T60 <- 1}
       init
     }
     
    dmg_fit <- sampling(object = stanmods[[ii]],
                        data = list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                                    N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                                    N_x = length(T100_data),N_y = length(R100_data),
                                    X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                                    X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                                    t_x = R100_data,t_y = T100_data,
                                    l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                                    l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3]),
                        control = list(adapt_delta = 0.8),
                        refresh = 0,init = init_value)
    
    
    loo_fit[[ii]] <- loo(dmg_fit)
    rho_est[ii,] <- summary(dmg_fit,pars = 'rho')$summary[,c('mean','2.5%','97.5%')]

   }
  loocomp <- loo_compare(loo_fit[[1]],loo_fit[[2]],loo_fit[[3]],loo_fit[[4]],loo_fit[[5]],loo_fit[[6]],loo_fit[[7]],loo_fit[[8]],loo_fit[[9]],loo_fit[[10]],loo_fit[[11]],loo_fit[[12]],loo_fit[[13]],loo_fit[[14]],loo_fit[[15]],loo_fit[[16]],loo_fit[[17]],loo_fit[[18]],loo_fit[[19]],loo_fit[[20]],loo_fit[[21]],loo_fit[[22]],loo_fit[[23]],loo_fit[[24]],loo_fit[[25]],loo_fit[[26]],loo_fit[[27]],loo_fit[[28]],loo_fit[[29]],loo_fit[[30]],loo_fit[[31]],loo_fit[[32]],loo_fit[[33]],loo_fit[[34]],loo_fit[[35]],loo_fit[[36]],loo_fit[[37]],loo_fit[[38]],loo_fit[[39]],loo_fit[[40]],loo_fit[[41]],loo_fit[[42]],loo_fit[[43]],loo_fit[[44]],loo_fit[[45]],loo_fit[[46]],loo_fit[[47]],loo_fit[[48]],loo_fit[[49]],loo_fit[[50]],loo_fit[[51]],loo_fit[[52]],loo_fit[[53]],loo_fit[[54]],loo_fit[[55]],loo_fit[[56]],loo_fit[[57]],loo_fit[[58]],loo_fit[[59]],loo_fit[[60]],loo_fit[[61]],loo_fit[[62]],loo_fit[[63]],loo_fit[[64]])

  loo_value[[jj]] <- loocomp
  rho_value[[jj]] <- rho_est
  
  saveRDS(loo_value,file = 'loo_fit_partial_2times.rds')
  saveRDS(rho_value,file = 'rho_2times.rds')
  
}








