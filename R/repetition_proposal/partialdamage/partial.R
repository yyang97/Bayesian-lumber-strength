require(rstan)
require(MASS)
require(loo)
nsims <- 1000
# with test, it has an additional term 1/p(X>l)
# sim_mean <- matrix(NA, ncol = 12, nrow = nsims)
# sim_upper <- matrix(NA, ncol = 12, nrow = nsims)
# sim_sd <- matrix(NA, ncol = 12, nrow = nsims)
# sim_lower <- matrix(NA, ncol = 12, nrow = nsims)

sim_mean_test <- matrix(NA, ncol = 12, nrow = nsims)
sim_upper_test <- matrix(NA, ncol = 12, nrow = nsims)
sim_sd_test <- matrix(NA, ncol = 12, nrow = nsims)
sim_lower_test <- matrix(NA, ncol = 12, nrow = nsims)



for(isim in 1:nsims){
  # show(isim)
  mu <- c(45,5.5)
  sigma <- c(13,1)
  
  rho <- 0.7
  
  
  N <- 87
  
  # alpha[6] alpha for R20, R40,R60, T20,T40,T60
  # alpha <- c(.1,.2,.3,.4,.5,.6) 
  
  # alpha <- c(.9,.8,.7,.9,.8,.7) 
  # alpha <- c(1,1,1,1,1,1) 
   alpha <- c(.7,1,.9,.7,1,.9) 
  
  # prop is thresh 
  thresh <- 0.7
  prop <- thresh
  
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
  
 
  dmgtest_mod <- stan_model("newmodel_thresh.stan")
  options(mc.cores = parallel::detectCores())
  
  init_dmg <- function() {
    list(mu = c(45,5), sigma = c(13,1), rho = .7, alpha = rep(1,6),thresh = .7 )
  }
  
  data_stan <- list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                    N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                    N_x = length(T100_data),N_y = length(R100_data),
                    X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                    X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                    t_x = R100_data,t_y = T100_data,
                    l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                    l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3])
  
  realtest_fit <- sampling(object = dmgtest_mod,
                           data = data_stan,init = init_dmg, iter = 5000)
  # print(realtest_fit, par = c('mu', 'sigma','rho','alpha','thresh'))
  sim_mean_test[isim,] <- summary(realtest_fit, par = c('mu', 'sigma','rho','alpha','thresh'))$summary[,1]
  sim_sd_test[isim,] <- summary(realtest_fit, par = c('mu', 'sigma','rho','alpha','thresh'))$summary[,3]
  sim_lower_test[isim,] <- summary(realtest_fit, par = c('mu', 'sigma','rho','alpha','thresh'))$summary[,4]
  sim_upper_test[isim,] <- summary(realtest_fit, par = c('mu', 'sigma','rho','alpha','thresh'))$summary[,8]
  simtest_res <- list(mean = sim_mean_test, sd = sim_sd_test, upper = sim_upper_test, lower = sim_lower_test)
  saveRDS(simtest_res, file = "sim_partial.rds")
  
  
  
 
  
  
    show(isim)
}

