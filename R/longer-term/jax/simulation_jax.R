source("dmg_func.R")
## plot 
set.seed(8888)
# the orignal sample size 
N <- 300
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
alpha <- 0.24
c <- 0.65
s <- 1

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

y_obs_g1 <- y_obs$group1
y_obs_g2 <- y_obs$group2
y_obs_g3 <- y_obs$group3




negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik2-lik3)
  #return (-lik3)
}

theta0 <- c(mu,sigma,alpha,c)
optimCheck::optim_proj(theta0,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))


optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
                  lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimCheck::optim_proj(optimout$par,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))

## using jax 
negdmglik_jax <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- dmglik_jax(y_obs$group3,alpha,l,c,s,mu,sigma)
  return(-lik1-lik2-lik3)
  #return (-lik3)
}

negdmglik_jax(theta0)
optimout <- optim(theta0,negdmglik_jax,method = "L-BFGS-B",
                  lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimCheck::optim_proj(optimout$par,
                       negdmglik_jax,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))



## py 
negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik2-lik3)
  #return (-lik3)
}


optimout <- optim(theta0,
                  negdmglik_np,
                  y_obs_g1 = y_obs_g1,
                  y_obs_g2 = y_obs_g2,
                  y_obs_g3 = y_obs_g2,l = l,s = s)

optimCheck::optim_proj(optimout$par,
                       negdmglik_py,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))

##--------python function test------

# test the sigmoid function 
sigmoid(2,1)
smooth_ind(2,1)

# test the dmgmodel 
dmg_model(y = 37,alpha = .45,l = 32, c = .65,s = 1)
dmgmodel_py(y = 37,alpha = .45,l = 32, c = .65,s = 1)


# test the yinverse 
dmg_inverse(ystar = 37,alpha = .45,l = 32, c = .65,s = 1)
dmginverse_py(ystar = 37,alpha = .45,l = 32, c = .65,s = 1)

# test the gradient 
dmg_inverse_grad(ystar = 37,alpha = .45,l = 32, c = .65,s = 1)
dmginvgrad_py(ystar = 37,alpha = .45,l = 32, c = .65,s = 1)


# dmglik 
dmglik(ystar = 37,alpha = .45,l = 32, c = .65,s = 1, mu = 48, sigma = 19)
dmglik_py(ystar = 37,alpha = .45,l = 32, c = .65,s = 1, mu = 48, sigma = 19)
dmglik_vmap(y_group = y_obs$group3,alpha = .45,l = 32, c = .65,s = 1, mu = 48, sigma = 19)
dmglik_np(y_group = y_obs$group3,alpha = .45,l = 32, c = .65,s = 1, mu = 48, sigma = 19)


sum(sapply(y_obs$group3,dmglik,alpha = .45,l = 32, c = .65,s = 1, mu = 48, sigma = 19))
#sum(sapply(y_obs$group3,dmglik_py,alpha,l,c,s,mu,sigma))



## neglik
y_obs
lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
lik2 <- y_obs$group2 * log(
  pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
    pnorm(l, mean = mu, sd = sigma)
)
lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))

-lik1-lik2-lik3

negdmglik_np(theta0,y_obs_g1,y_obs_g2,y_obs_g3,l,s)
