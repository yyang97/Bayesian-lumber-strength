source("dmg_func.R")

set.seed(8888)
# the orignal sample size 
N <- 300
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
# alpha <- 0.9
# c <- 0.8
# alpha <- 0.9
alpha <- 0.85
c <- 0.68
s <- 10

logit <- function(x){
  return(log(x/(1-x)))
}

expit <- function(x){
  1/(1+exp(-x))
}
## data 

y <- rnorm(N, mean = mu, sd = sigma )
y <- y[y>0]
y_obs <- list()
# group 1, y^* <y < l
y_obs$group1 <- y[which(y < l)]
# group 2, y^*<l<y
group23 <- y[which(y > l)]
group23_star <- dmg_model(group23,alpha,l,c,s)
y_obs$group2 <- length(group23[which(group23 > l &group23_star < l)])
y_obs$group3 <- group23_star[which(group23_star >l)]

negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  # alpha <- expit(theta[3])
  # c <- expit(theta[4])
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <-dmglik_jax(y_obs$group3,alpha,l,c,s,mu,sigma)
  return(-lik1-lik2-lik3)
  #return (-lik3)
}

#theta0 <- c(mu,sigma,logit(alpha),logit(c))
theta0 <- c(mu,sigma,alpha,c)
# optimCheck::optim_proj(theta0,
#                        negdmglik_model,xrng = .5,
#                        xnames = c("mu","sigma","alpha","c"))


# optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
#                   lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimout <- optim(theta0,negdmglik_model)

optimout$par
optimCheck::optim_proj(optimout$par,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))

# fisher<- numDeriv::hessian(negdmglik_model,optimout$par,
#                            method.args=list(r=6))
# theta_se <- sqrt(diag(chol2inv(chol(fisher))))
# theta_se 
# cbind(theta0,optimout$par,optimout$par - 1.96*theta_se,optimout$par + 1.96*theta_se )

# hypothesis test
# theta_hat <- optimout$par
# mu_hat <- theta_hat[1]
# sigma_hat <- theta_hat[2]
# alpha_hat <- expit(theta_hat[3])
# c_hat <- expit(theta_hat[4])
# 
# 
# optimCheck::optim_proj(c(mu_hat,sigma_hat,alpha_hat,c_hat),
#                        negdmglik_model,xrng = .5,
#                        xnames = c("mu","sigma","alpha","c"))

# test a ==1 
mu_truc <- truncnorm::etruncnorm(a = l, b = Inf, mean = mu_hat, sd = sigma_hat)
var_truc <- truncnorm::vtruncnorm(a = l, b = Inf, mean = mu_hat, sd = sigma_hat)


sample_mean <- mu_truc
sample_var <- 1/length(y_obs$group3)*var_truc
test_stat <- (mean(y_obs$group3) - sample_mean)/sqrt(sample_var)
test_stat
pnorm(test_stat)


# test a<= c
mu_truc <- truncnorm::etruncnorm(a = l/c_hat, b = Inf, mean = mu_hat, sd = sigma_hat)
var_truc <- truncnorm::vtruncnorm(a = l/c_hat, b = Inf, mean = mu_hat, sd = sigma_hat)


sample_mean <- mu_truc
sample_var <- 1/length(y_obs$group3)*var_truc
test_stat <- (mean(y_obs$group3) - sample_mean)/sqrt(sample_var)
pnorm(test_stat)
#pnorm(mean(y_obs$group3), mean = sample_mean, sd = sqrt(sample_var))



# 
# 
# N <- 300
# mu <- 48
# sigma <- 19
# #y <- rnorm(N,mean = mu, sd = sigma)
# l <-  32
# # alpha <- 0.6
# c <- 0.68
# s <- 10
# 
# 
# alpha_seq <- seq(from = 0.1, to = 0.99, by = 0.01)
# res <- matrix(NA, nrow = length(alpha_seq), ncol =5)
# colnames(res) <- c("alpha","c","p_val","alpha_hat","c_hat")
# for (ii in 46:length(alpha_seq)){
#   # data generating
#   print(ii)
#   alpha <- alpha_seq[ii]
#   set.seed(8888)
#   y <- rnorm(N, mean = mu, sd = sigma )
#   y <- y[y>0]
#   y_obs <- list()
#   # group 1, y^* <y < l
#   y_obs$group1 <- y[which(y < l)]
#   # group 2, y^*<l<y
#   group23 <- y[which(y > l)]
#   group23_star <- dmg_model(group23,alpha,l,c,s)
#   y_obs$group2 <- length(group23[which(group23 > l &group23_star < l)])
#   y_obs$group3 <- group23_star[which(group23_star >l)]
#   
#   # fitting the model
#   theta0 <- c(mu,sigma,logit(alpha),logit(c))
#   optimout <- optim(theta0,negdmglik_model)
#   # ,method = "L-BFGS-B",
#   #                   lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))
#   
#   # parameter estimation
#   theta_hat <- optimout$par
#   mu_hat <- theta_hat[1]
#   sigma_hat <- theta_hat[2]
#   alpha_hat <- expit(theta_hat[3])
#   c_hat <- expit(theta_hat[4])
#   
#   # hypothesis testing
#   mu_truc <- truncnorm::etruncnorm(a = l/c_hat, b = Inf, mean = mu_hat, sd = sigma_hat)
#   var_truc <- truncnorm::vtruncnorm(a = l/c_hat, b = Inf, mean = mu_hat, sd = sigma_hat)
#   
#   
#   sample_mean <- mu_truc
#   sample_var <- 1/length(y_obs$group3)*var_truc
#   # p_val is the probability of accepting H_0 when H_0 is true
#   p_val <- pnorm(mean(y_obs$group3), mean = sample_mean, sd = sqrt(sample_var))
#   res[ii,1] <- alpha
#   res[ii,2] <- c
#   res[ii,3] <- p_val
#   res[ii,4] <- alpha_hat
#   res[ii,5] <- c_hat
#   print(res[ii,])
# }
# p_val_matrix <- saveRDS(res, file = "p_val_matrix.rds")