source("dmg_func.R")
##------------------test the smooth function-----
N <-300
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
# alpha <- 0.9
# c <- 0.8
alpha <- 0.45
c <- 0.65
s <- 1


x <- seq(from = -10, to= 10, length.out = 500)
ind_y <- ifelse(x>0,1,0)
plot(x,ind_y,type = "l")
# higher s, closer to the indicator function 
lines(x,sapply(x,smooth_ind, s = 1),type = "l",col = "red")
lines(x,sapply(x,smooth_ind, s = 10),type = "l",col = "blue")
lines(x,sapply(x,smooth_ind, s = 100),type = "l",col = "green")
legend("topleft", legend=c("Indicator", "s = 1", "s = 10", "s = 100"),
       col=c("black","red", "blue","green"), lty=1, cex=0.8)



# y and ystar indicator -------
y <- seq(from  = 0.01, to = 100, length = 200)
ystar <- ifelse(c*y>l,y,alpha*y)

plot(y,ystar, type = "l")


#------------y and ystar plot--------------
alpha <- 0.25
c <- 0.65
y <- seq(from  = 0.01, to = 100, length = 200)
ystar <- dmg_model(y,alpha,l,c,s =1)
plot(y,ystar,type = 'l',lty = 1)
lines(y,dmg_model(y,alpha = 0.45,l,c,s =1),lty = 2)
legend("topleft", legend = c("ystar alpha = 0.45","ystar alpha = 0.25"),lty=1:2)
#abline(v = l)
abline(v = l/c)



#-----------ystar and y plot-------
ystar <- seq(from  = 0.1, to = 100, length = 200)
y <- sapply(ystar,dmg_inverse,alpha = .9,l = l,c = c ,s = s)
plot(ystar,y,type = 'l',lty = 1)
lines(ystar,ystar,lty = 2)
legend("topleft", legend = c("y","ystar"),lty=1:2)
#abline(v = l)
abline(v = l/c)



## test the gradient
set.seed(8888)
N <- 1000
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
alpha <- 0.45
c <- 0.65
s <- 1

dmg_inverse_grad(40,alpha,l,c,s)
# 
# eps <- 1e-6
# (dmg_inverse(40+eps,alpha,l,c,s)-dmg_inverse(40-eps,alpha,l,c,s))/(2*eps)
# 
# 
# 
# dmg_inverse_numgrad <- function(ystar,alpha,l,c,s,eps){
#   abs((dmg_inverse(ystar+eps,alpha,l,c,s) -dmg_inverse(ystar-eps,alpha,l,c,s))/(2*eps))
# }  
# eps <- 1e-6
# dmg_inverse_grad(40,alpha,l,c,s)
# 
# dmg_inverse_numgrad(40,alpha,l,c,s,eps )

N <- 300
## suppose all observations are damaged 
set.seed(8888)
y <- rnorm(N, mean = mu, sd = sigma )
y <- y[y>0]
ystar <- dmg_model(y,alpha,l,c,s)

negg3 <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
  alpha <- theta[3]
  c <- theta[4]
  lik3 <- sum(sapply(ystar,dmglik,alpha,l,c,s,mu,sigma))
  return (-lik3)
}

theta0 <- c(47.9985218, 18.9804222,(0.4548490),  0.6515348)

negg3(theta0)

optimCheck::optim_proj(theta0,
                       negg3,
                       xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))

## test the simulaiton 





## plot 
set.seed(8888)
# the orignal sample size 
N <- 1000
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
# alpha <- 0.24
alpha <- 0.45
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
# increasing order now 
#y_obs$group3 <- y_obs$group3[order(y_obs$group3)]


theta0 <- c(mu,sigma,alpha,c)
optimCheck::optim_proj(theta0,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))
# 
# 
optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
                  lower = rep(0.1,4),upper = c(Inf,Inf,1,Inf))

optimCheck::optim_proj(optimout$par,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","alpha","c"))

negdmglik_alpha<- function(theta){

  alpha <- theta

  
  lik1 <- sum(dnorm(y_obs$group1,mean = mu, sd = sigma, log = TRUE))
  lik2 <- y_obs$group2 * log(
    pnorm(dmg_inverse(l,alpha,l,c,s), mean = mu, sd = sigma) -
      pnorm(l, mean = mu, sd = sigma)
  )
  #lik2 <- 0
  lik3 <- sum(sapply(y_obs$group3,dmglik,alpha,l,c,s,mu,sigma))
  return(-lik1-lik2-lik3)
}

alpha_seq <- seq(from = 0.1, to = 0.5, length = 200)
dmglik_value <- sapply(alpha_seq,negdmglik_alpha)
# dmg_num <- sapply(alpha_seq,function(alpha){log(dmg_inverse_numgrad(y_obs$group3[50],alpha,l,c,s,eps = 0.5))
# })
plot(alpha_seq,dmglik_value,type = "l")


optimize(negdmglik_alpha,c(.1,.2))
optimize(negdmglik_alpha,c(.15,.2))


alpha_seq <- seq(from = 0.15, to = 0.2, length = 100)
dmglik_value <- sapply(alpha_seq,negdmglik_alpha)
plot(alpha_seq,dmglik_value,type = "l")

optimize(negdmglik_alpha,c(.175,.182))

alpha_seq <- seq(from = 0.175, to = 0.182, length = 100)
dmglik_value <- sapply(alpha_seq,negdmglik_alpha)
plot(alpha_seq,dmglik_value,type = "l")

optimize(negdmglik_alpha,c(.179,.180))

which.min(dmglik_value)
alpha_seq[69]


negdmglik_alpha(0.1798) -negdmglik_alpha(0.18)


log(
  pnorm(dmg_inverse(l,0.1798,l,c,s), mean = mu, sd = sigma) -
    pnorm(l, mean = mu, sd = sigma)
)

log(
  pnorm(dmg_inverse(l,0.18,l,c,s), mean = mu, sd = sigma) -
    pnorm(l, mean = mu, sd = sigma)
)

# 0.44 diff
sum(sapply(y_obs$group3,dmglik,0.1798,l,c,s,mu,sigma)[-50])-sum(sapply(y_obs$group3,dmglik,0.180,l,c,s,mu,sigma)[-50])

sapply(y_obs$group3,dmglik,0.1798,l,c,s,mu,sigma)-sapply(y_obs$group3,dmglik,0.180,l,c,s,mu,sigma)


order(sapply(y_obs$group3,dmglik,0.1798,l,c,s,mu,sigma)-
  sapply(y_obs$group3,dmglik,0.180,l,c,s,mu,sigma),decreasing = TRUE )   


(sapply(y_obs$group3,dmglik,0.1798,l,c,s,mu,sigma)-
  sapply(y_obs$group3,dmglik,0.180,l,c,s,mu,sigma))[10]

y_obs$group3[50]

dmglik(36.09944,0.1798,l,c,s,mu,sigma)
dmglik(36.09944,0.18,l,c,s,mu,sigma)

log(dmg_inverse_grad(36.09944,0.1798,l,c,s))
log(dmg_inverse_grad(36.09944,0.18,l,c,s))



alpha_seq <- seq(from = 0.01, to = 0.99, length = 200)

sapply(alpha_seq, function(alpha)(log(dmg_inverse_grad(36.09944,alpha,l,c,s))))
plot(alpha_seq ,
     sapply(alpha_seq, function(alpha)(log(dmg_inverse_grad(36.09944,alpha,l,c,s)))),
     type = "l",ylab = "neglik"
)


alpha_seq <- seq(from = 0.1, to = 0.3, length = 200)

sapply(alpha_seq, function(alpha)(log(dmg_inverse_grad(36.09944,alpha,l,c,s))))
plot(alpha_seq ,
     sapply(alpha_seq, function(alpha)(log(dmg_inverse_grad(36.09944,alpha,l,c,s)))),
     type = "l",ylab = "neglik"
)


# the inverse function 
dmg_inverse <- function(ystar,alpha,l,c,s){
  uniroot((function (x) dmg_model(x,alpha,l,c,s) - ystar), lower = 0, upper = 1000)$root
}

plot(alpha_seq,
sapply(alpha_seq, dmg_inverse,ystar =36.0994, l,c,s),
type = "l", xlab = "alpha", ylab = "y")



dmg_inverse(ystar,alpha,l,c,s)
dnorm(y, mean = mu, sd = sigma, log = TRUE)


# test the first part
alpha_seq <- seq(from = 0.15, to = 0.25, length = 100)
plot(alpha_seq, 
     sapply(alpha_seq, function(alpha){dnorm(dmg_inverse(y_obs$group3[50],alpha,l,c,s), mean = mu, sd = sigma, log = TRUE)
}),
type = "l",xlab = "alpha", ylab = "normal part")


log(dmg_inverse_grad(dmg_inverse(y_obs$group3[50],alpha,l,c,s),alpha,l,c,s))



# test the gradient part
alpha_seq <- seq(from = 0.15, to = 0.25, length = 100)
plot(alpha_seq, 
     sapply(alpha_seq, function(alpha){log(dmg_inverse_grad(y_obs$group3[50],alpha,l,c,s))
}),
     type = "l",xlab = "alpha", ylab = "normal part")

plot(alpha_seq,
     sapply(alpha_seq, function(alpha){log(dmg_inverse_numgrad(y_obs$group3[50],alpha,l,c,s,eps = 0.5))
     }),
     type = "l",xlab = "alpha", ylab = "normal part")


gradient_test_func <- function(alpha){
  log(dmg_inverse_grad(y_obs$group3[50],alpha,l,c,s))
  
}




optimize(gradient_test_func,c(.18,..19), maximum = TRUE)

alpha1 <- 0.1964646
alpha2 <-  0.1944444
gradient_test_func(  0.1944444)




#-----------ystar and y plot-------
ystar <- seq(from  = 36, to = 37, length = 200)
y <- sapply(ystar,dmg_inverse,alpha = alpha1,l = l,c = c ,s = s)
plot(ystar,y,type = 'l',lty = 1,ylim = c(21,55))
lines(ystar,ystar,lty = 2)
legend("topleft", legend = c("y","ystar"),lty=1:2)
l/c

dmg_inverse_numgrad <- function(ystar,alpha,l,c,s,eps){
  abs((dmg_inverse(ystar+eps,alpha,l,c,s) -dmg_inverse(ystar-eps,alpha,l,c,s))/(2*eps))
}  

eps <- 1e-2
log((dmg_inverse(y_obs$group3[50]+eps,alpha2,l,c,s) - dmg_inverse(y_obs$group3[50]-eps,alpha2,l,c,s))/(2*eps))
