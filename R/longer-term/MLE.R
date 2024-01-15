require(MASS)

require(optimCheck)



mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.8


N <- 300

# alpha[6] alpha for R20, R40,R60, T20,T40,T60
# alpha <- c(.1,.2,.3,.4,.5,.6) 


# alpha <- c(.8,.95,.85,.8,.95,.85) 
alpha <- c(.8,.85,.9,.8,.85,.9) 

# eta is thresh 
eta <- 0.7


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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta, alpha[1],1)
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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta, alpha[2],1)
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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta, alpha[3],1)
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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta, alpha[4],1)
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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta, alpha[5],1)
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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta, alpha[6],1)
    PFY_ob[i,3] <- 0
  }
}

T60_data <- PFY_ob

##-----T100-----######

R100_data <- rnorm(2*N,mean = mu[1],sd = sigma[1])
T100_data <- rnorm(2*N,mean = mu[2],sd = sigma[2])




int_function_R <- function(ystar, mu,sigma, rho,l){

  a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));

  dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)


}


PFY_lik_R <- function(mu, sigma, rho,alpha,l,data){

  (sum(dnorm(data[data[,3] == 1,1],
             mean = mu[1],
             sd = sigma[1],
             log = T)
  ) +
    sum(log(1/alpha*(int_function_R(data[data[,3] == 0,2]/alpha,mu,sigma,rho,l)-
                       int_function_R(data[data[,3] == 0,2]/alpha,mu,sigma,rho,1/eta*l)) +  int_function_R(data[data[,3] == 0,2],mu,sigma,rho,1/eta*l)
    ))


  )



}



int_function_T <- function(ystar, mu,sigma, rho,l){

  a_l <- (l - mu[2]-rho*(sigma[2]/sigma[1])*(ystar-mu[1]))/(sigma[2]*sqrt(1-rho^2));

  dnorm(ystar, mean = mu[1], sd = sigma[1],log = F)*pnorm(a_l,lower.tail = F, log.p = F)


}


PFY_lik_T <- function(mu, sigma, rho,alpha,l,data){

  (sum(dnorm(data[data[,3] == 1,1],
             mean = mu[2],
             sd = sigma[2],
             log = T)
  ) +
    sum(log(1/alpha*(int_function_T(data[data[,3] == 0,2]/alpha,mu,sigma,rho,l)-
                       int_function_T(data[data[,3] == 0,2]/alpha,mu,sigma,rho,1/eta*l)) +  int_function_T(data[data[,3] == 0,2],mu,sigma,rho,1/eta*l)
    ))


  )



}
# sum(dnorm(data[data[,3] == 1,1], mean = mu[1], sd = sigma[1], log = T))

nlogpost <- function(theta){
  mu <- theta[1:2]
  sigma <- theta[3:4]
  rho <- theta[5]
  alpha <- theta[6:11]
  lik <- PFY_lik_R(mu,sigma,rho,alpha[1],R_pf[1],R20_data) +
    PFY_lik_R(mu,sigma,rho,alpha[2],R_pf[2],R40_data) +
    PFY_lik_R(mu,sigma,rho,alpha[3],R_pf[3],R60_data) +
    PFY_lik_T(mu,sigma,rho,alpha[4],T_pf[1],T20_data) +
    PFY_lik_T(mu,sigma,rho,alpha[5],T_pf[2],T40_data) +
    PFY_lik_T(mu,sigma,rho,alpha[6],T_pf[3],T60_data) +
    sum(dnorm(T100_data,
          mean = mu[2],
          sd = sigma[2], log = T))+
    sum(dnorm(R100_data,
              mean = mu[1],
              sd = sigma[1], log = T))
  return(-1*lik)
}

nlogpost(c(mu,sigma,rho,alpha))




optimresult <- optim(c(mu,sigma,rho,alpha),nlogpost,method = 'L-BFGS-B',
                   lower = rep(0,11), upper = c(rep(Inf,4), 1,rep(Inf,6)))
#optimresult <- optim(c(mu,sigma,rho,alpha[1:3]),nlogpost)
optimresult$par


optim_proj(fun = nlogpost, # objective function
           xsol = optimresult$par, # candidate optimum
           maximize = TRUE, # look for max
           xnames = c("mu_x","mu_y","sig_x","sig_y","rho","alpha")) # pretty names on plot



optimresult$par[5]

pdf_condy <- function(mu, sigma,rho,l){
  a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
  
  dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)/(1-pnorm((l - mu[1])/sigma[1]))
  
}


