realdata <- read.csv("bending-pl-4groups.csv",header = FALSE)
fac = factor(realdata[,2])
#nlevels(fac)
R20_4Y_R100 = realdata[which(realdada[,2] == "R20_4Y_R100"),1]
R20_4Y_R100 = R20_4Y_R100/1000*6.895

R20_4Y = realdata[which(realdada[,2] == "R20_4Y"),1]
R20_4Y = R20_4Y/1000*6.895
l = 4500/1000*6.895




negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(R20_4Y,mean = mu, sd = sigma, log = TRUE))
  lik2 <- 41 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(R20_4Y_R100,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik2-lik3)
  #return (-lik3)
}

theta0 <- c(mu,sigma,alpha,c)

optimout <- optim(theta0,negdmglik_model)
optimout$par
theta_est <- optimout$par





mu_hat <- 45.77
sigma_hat <- 17.94
c_hat <- 0.712


mu_truc <- truncnorm::etruncnorm(a = l/c_hat, b = Inf, mean = mu_hat, sd = sigma_hat)
var_truc <- truncnorm::vtruncnorm(a = l/c_hat, b = Inf, mean = mu_hat, sd = sigma_hat)


sample_mean <- mu_truc
sample_var <- 1/length(R20_4Y_R100 )*var_truc
test_stat <- (mean(R20_4Y_R100 ) - sample_mean)/sqrt(sample_var)
pnorm(test_stat)



s = 10

negdmglik_model <- function(alpha){
  mu <- 45.77
  sigma <- 17.94
  c <- 0.712
  
  lik1 <- sum(dnorm(R20_4Y,mean = mu, sd = sigma, log = TRUE))
  lik2 <- 41 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(R20_4Y_R100,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik2-lik3)
  #return (-lik3)
}
# 0.712
negdmglik_model(.71)

optimize(negdmglik_model,c(.5,.9))

alpha_seq <- seq(from = 0.5,to =0.9, by = 0.001)
obj <- sapply(alpha_seq, negdmglik_model)
plot(alpha_seq,obj, type = 'l')
