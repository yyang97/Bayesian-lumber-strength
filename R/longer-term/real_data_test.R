
# # real data 
# #setwd("~/Downloads/github/Bayesian-lumber-strength/R/longer-term")
bending <- read.csv("bending-pl.csv", header = T)
# convert psi to Mpa
#  1 thousand psi = 6.895 MPa
bending[,1] <- bending[,1]/1000*6.895
l <- 4500/1000*6.895
# bending_mar <- bending[1:195,]
# bending_dmg <- bending[195:341,]

R100_data <- bending[bending[,2] == "R100",1]
R20_data <- bending[bending[,2] == "R20",1]
R20R100_data <- bending[bending[,2] == "R20R100",1]



mu <-  47.326717
sgima <- 18.393412
c <- 0.68320185

mu_truc <- truncnorm::etruncnorm(a = l/c, b = Inf, mean = mu, sd = sigma)
var_truc <- truncnorm::vtruncnorm(a = l/c, b = Inf, mean = mu, sd = sigma)

sample_mean <- mu_truc
sample_var <- 1/length(R20R100_data)*var_truc



#pnorm(mean(R20R100_data), mean = sample_mean, sd = sqrt(sample_var))

negdmglik_model <- function(theta){
  mu <- theta[1]
  sigma <- theta[2]
 # alpha <- theta[3]
  c <- theta[3]
  
  lik1 <- sum(log(truncnorm::dtruncnorm(R20_data,
                                    a=-Inf,
                                    b=l, 
                                    mean = mu, sd = sigma)
))
  lik2 <- sum(dnorm(R100_data,mean = mu, sd = sigma, log = TRUE))
  
  #lik2 <- 0
  lik3 <- sum(log(truncnorm::dtruncnorm(R20R100_data,
                                        a=l/c,
                                        b=Inf, 
                                        mean = mu, sd = sigma)
  ))
  return(-lik1-lik2-lik3)
}

theta0 <- c(48,19,.9)
negdmglik_model(theta0 )
optimCheck::optim_proj(theta0,
                       negdmglik_model,xrng = .5,
                       xnames = c("mu","sigma","c"))

optimout <- optim(theta0,negdmglik_model,method = "L-BFGS-B",
                  lower = rep(0.1,3),upper = c(Inf,Inf,1))


log(truncnorm::dtruncnorm(R20_data,
                          a=-Inf,
                          b=l, 
                          mean = 48, sd = 19))
    

log(truncnorm::dtruncnorm(R20R100_data,
                      a=l/.85,
                      b=Inf, 
                      mean = 48, sd = 19))

    