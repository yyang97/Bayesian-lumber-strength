##-----------set environment and check--------------------
#setwd("~/Downloads/github/Bayesian-lumber-strength/R/Simulation")
require(rstan)
require(MASS)
require(loo)

load('allSTANcomp.rda')
rstan_options(auto_write = TRUE)


nsim <- 100
loo_chosen <- c()
loo_value <- matrix(NA, nrow = nsim, ncol = 64)
for (jj in 1:nsim){
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
  
  # alpha_R <- 10* c(0.67,1.83,1.29)
  # alpha_T <- 10*c(85.23,104.09,65.95)
  
  alpha_R <- 10* c(0.67,1.83,1.29)
  alpha_T <- 10*c(85.23,104.09,65.95)

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
    
    
    loo_fit <- loo(dmg_fit)
    loo_value[jj,ii] <- loo_fit[[1]][3,1]
   }
  
  

  
  save(loo_value,file = 'loo_chosen_full.rda')
}







