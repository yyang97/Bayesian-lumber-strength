# 
# 
# sd_x <- sigma[1]
# sd_y <- sigma[2]
# Sigma <- matrix(c(sd_x^2, sd_x*sd_y*rho, sd_x*sd_y*rho, sd_y^2),2)
# biv <- mvrnorm(10000, mu = c(mu[1],mu[2]), Sigma = Sigma )
# 
# 
# id <- biv[,1] > R_pf[1]
# mean(biv[id,2])
# var(biv[id,2])
# 
# sigma[2]
# mu[2]
# R 20


alpha <- 0.95

mu <- c(45.13,5.51)
sigma <- c(12.92,1.08)

rho <- 0.8
eta <- 0.7


N <- 300


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
    PFY_ob[i,2] <- samples[i,2]*ifelse(l > samples[i,1]*eta,alpha,1)
    PFY_ob[i,3] <- 0
  }
}

R20_data <- PFY_ob

ystar <- R20_data[R20_data[,3] == 0,2]

mean(ystar)

pdf_condy <- function(ystar,mu, sigma,rho,l){
  a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
  
  dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)/(1-pnorm((l - mu[1])/sigma[1]))
  
}


expect_condy <- function(ystar,mu, sigma,rho,l){
  a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
  
  ystar*dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)/(1-pnorm((l - mu[1])/sigma[1]))
  
}
#' r is moment number number
#' r = 1 means E(y)
#' r = 2 means E(y^2)
moment_condy <- function(ystar,mu, sigma,rho,l, r){
  a_l <- (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
  
  ystar^(r)*dnorm(ystar, mean = mu[2], sd = sigma[2],log = F)*pnorm(a_l,lower.tail = F, log.p = F)/(1-pnorm((l - mu[1])/sigma[1]))
  
}


# expectation 
condy_mean <- stats::integrate(moment_condy,lower = -Inf, upper = Inf, mu = mu, sigma = sigma,
                 rho = rho, l = R_pf[1], r  = 1)$value
# var 
condy_var <- stats::integrate(moment_condy,lower = -Inf, upper = Inf, mu = mu, sigma = sigma,
                              rho = rho, l = R_pf[1], r  = 2)$value - condy_mean^2
# condy_mean
# condy_var

pval <- pnorm((mean(ystar) - condy_mean)/sqrt(condy_var/length(ystar)))

pval 
pval < 0.05

