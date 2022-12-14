require(rstan)
require(MASS)
require(loo)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#load('allSTAN.rda')


load("real_data.RData")

loo_fit <- list()

allcombs <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
for (ii in 1:64) {
  show(paste('model', ii))
  
  ci <- allcombs[ii,]
  # initial value 
  init_value <- function(){
    init <- list (mu = c(35,8), sigma = c(10,1), rho = .8)
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
  
}
loocomp <- loo_compare(loo_fit[[1]],loo_fit[[2]],loo_fit[[3]],loo_fit[[4]],loo_fit[[5]],loo_fit[[6]],loo_fit[[7]],loo_fit[[8]],loo_fit[[9]],loo_fit[[10]],loo_fit[[11]],loo_fit[[12]],loo_fit[[13]],loo_fit[[14]],loo_fit[[15]],loo_fit[[16]],loo_fit[[17]],loo_fit[[18]],loo_fit[[19]],loo_fit[[20]],loo_fit[[21]],loo_fit[[22]],loo_fit[[23]],loo_fit[[24]],loo_fit[[25]],loo_fit[[26]],loo_fit[[27]],loo_fit[[28]],loo_fit[[29]],loo_fit[[30]],loo_fit[[31]],loo_fit[[32]],loo_fit[[33]],loo_fit[[34]],loo_fit[[35]],loo_fit[[36]],loo_fit[[37]],loo_fit[[38]],loo_fit[[39]],loo_fit[[40]],loo_fit[[41]],loo_fit[[42]],loo_fit[[43]],loo_fit[[44]],loo_fit[[45]],loo_fit[[46]],loo_fit[[47]],loo_fit[[48]],loo_fit[[49]],loo_fit[[50]],loo_fit[[51]],loo_fit[[52]],loo_fit[[53]],loo_fit[[54]],loo_fit[[55]],loo_fit[[56]],loo_fit[[57]],loo_fit[[58]],loo_fit[[59]],loo_fit[[60]],loo_fit[[61]],loo_fit[[62]],loo_fit[[63]],loo_fit[[64]])

loo_value <- loocomp


saveRDS(loo_value,file = 'real_data_result.rds')


realcomp <- readRDS("real_data_result.rds")
#
#
realcomp
allcombs[3,]

allcombs[11,]


allcombs[3,]
# rho ~ U(0,1)
setwd("~/Downloads/github/Bayesian-lumber-strength/R")
realcomp <- readRDS("real_data_result.rds")
realcomp

# rho ~ U(0.7,0.85)
setwd("~/Downloads/github/Bayesian-lumber-strength/R/Real data analysis/rho07085")
realcomp_070085 <- readRDS("real_data_result.rds")
realcomp_070085

# rho ~ U(0,799, 0.801)
setwd("~/Downloads/github/Bayesian-lumber-strength/R/Real data analysis/rho07990801")
realcomp_0799081 <- readRDS("real_data_result.rds")
realcomp_0799081

# record the rho estimation of the three scenarios
# model 1 and model 3 

# Y* = a + bY






find_ci <- function(x, digits = 1) {
  round(c(lower = unname(x[1] - 2*x[2]),
          upper = unname(x[1] + 2*x[2])),
        digits = digits)
}
# check why 2 refers to 95\% interval in the paper 

t(apply(realcomp[-1, c("elpd_diff", "se_diff")], 1, find_ci))
t(apply(realcomp_070085[-1, c("elpd_diff", "se_diff")], 1, find_ci))
t(apply(realcomp_0799081[-1, c("elpd_diff", "se_diff")], 1, find_ci))

# mY* = mY - alpha/mY. alpha depends on the scale of Y 
# Scaling y will also scale alpha
# Y* = f(alpha, Y)
# fix mu, sigma, rho
# Y_i* = Y_i*I(L < 1/2*x_i) + alpha*Y_i*I(L > 1/2*x_i)
# p(x,y) = p(y|x)p(x)
# P(Y_i*, x_i > l) = p(y_i* | x_i > l)*p(x_i > l)


## A4 Q2


