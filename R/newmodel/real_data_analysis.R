require(rstan)
require(MASS)
require(loo)
library(readxl)
library(dplyr)
library(moments)
###------------------data preprocessing--------------
# setwd("~/Downloads/github/Bayesian-lumber-strength/R/Real data analysis")

summary_all08122013 <- read_excel("summary_all08122013.xlsx")

# summary_all08122013
summary_all08122013 <- select(summary_all08122013,"Group","Broken","MOR","UTS")

## T group 
T100_data <-summary_all08122013[which(summary_all08122013$Group == "T100"),]
T60_data <-summary_all08122013[which(summary_all08122013$Group == "T60"),]
T40_data <-summary_all08122013[which(summary_all08122013$Group == "T40"),]
T20_data <-summary_all08122013[which(summary_all08122013$Group == "T20"),]
T100_data$UTS <- as.numeric(T100_data$UTS)
T100_data <- T100_data$UTS
## G group 
R100_data <-summary_all08122013[which(summary_all08122013$Group == "R100"),]
R60_data <-summary_all08122013[which(summary_all08122013$Group == "R60"),]
R40_data <-summary_all08122013[which(summary_all08122013$Group == "R40"),]
R20_data <-summary_all08122013[which(summary_all08122013$Group == "R20"),]
R100_data <-as.numeric(R100_data$MOR)
R100_data 
##-----------------T substitute NA to 0 -----------

id <- T60_data$Broken == 1
T60_data$MOR[id] <- 0
T60_data$UTS[id] <- as.numeric(T60_data$UTS[id])
T60_data$UTS[!id] <- 0
T60_data$UTS <- as.numeric(T60_data$UTS)
T60_data$MOR <- as.numeric(T60_data$MOR)


id <- T40_data$Broken == 1
T40_data$MOR[id] <- 0
T40_data$UTS[id] <- as.numeric(T40_data$UTS[id])
T40_data$UTS[!id] <- 0
T40_data$UTS <- as.numeric(T40_data$UTS)
T40_data$MOR <- as.numeric(T40_data$MOR)


id <- T20_data$Broken == 1
T20_data$MOR[id] <- 0
T20_data$UTS[id] <- as.numeric(T20_data$UTS[id])
T20_data$UTS[!id] <- 0
T20_data$UTS <- as.numeric(T20_data$UTS)
T20_data$MOR <- as.numeric(T20_data$MOR)



##-----------------R substitute NA to 0 ------

id <- R60_data$Broken == 0
R60_data$MOR[id] <- 0
R60_data$UTS[id] <- as.numeric(R60_data$UTS[id])
R60_data$UTS[!id] <- 0
R60_data$UTS <- as.numeric(R60_data$UTS)
R60_data$MOR <- as.numeric(R60_data$MOR)



id <- R40_data$Broken == 0
R40_data$MOR[id] <- 0
R40_data$UTS[id] <- as.numeric(R40_data$UTS[id])
R40_data$UTS[!id] <- 0
R40_data$UTS <- as.numeric(R40_data$UTS)
R40_data$MOR <- as.numeric(R40_data$MOR)

id <- R20_data$Broken == 0
R20_data$MOR[id] <- 0
R20_data$UTS[id] <- as.numeric(R20_data$UTS[id])
R20_data$UTS[!id] <- 0
R20_data$UTS <- as.numeric(R20_data$UTS)
R20_data$MOR <- as.numeric(R20_data$MOR)




##---------------Convert psi to Mpa------------ 

# 1 thousand psi = 6.895 MPa

c = 6.895

R20_data$MOR <- c*R20_data$MOR
R20_data$UTS <- c*R20_data$UTS


R40_data$MOR <- c*R40_data$MOR
R40_data$UTS <- c*R40_data$UTS

R60_data$MOR <- c*R60_data$MOR
R60_data$UTS <- c*R60_data$UTS

R100_data <- c*R100_data

T20_data$MOR <- c*T20_data$MOR
T20_data$UTS <- c*T20_data$UTS

T40_data$MOR <- c*T40_data$MOR
T40_data$UTS <- c*T40_data$UTS

T60_data$MOR <- c*T60_data$MOR
T60_data$UTS <- c*T60_data$UTS


T100_data <- c*T100_data

##-------check normal fitting---------
# 
# library(fitdistrplus)
# 
# # check R100_data
# FIT <- fitdist(R100_data, "norm")    ## note: it is "norm" not "normal"
# plot(FIT)    ## use method `plot.fitdist`
# FIT$estimate # mean = 45.679 sd = 12.900
# FIT$bic # bic = 1394.022
# # good normal fitting 
# 
# 
# # check T100_data
# 
# shapiro.test(T100_data) # p = 0.0001841
# shapiro.test(log(T100_data)) # p-value = 0.002327
# shapiro.test(sqrt(T100_data)) # p-value = 0.3988
# 


##------proof loading-------------- 
R_pf <- c* c(4.956690733, 6.110714122, 7.092435407)
T_pf <- sqrt(c* c(2.962390379, 3.986497991, 4.916102264))


##--------convert UTS to sqrt(UTS)

R20_data$UTS <- sqrt(R20_data$UTS)
R40_data$UTS <- sqrt(R40_data$UTS)
R60_data$UTS <- sqrt(R60_data$UTS)
R100_data <- R100_data


T20_data$UTS <- sqrt(T20_data$UTS)
T40_data$UTS <- sqrt(T40_data$UTS)
T60_data$UTS <- sqrt(T60_data$UTS)
T100_data <- sqrt(T100_data)



##-----check outlier and delete it-------

# R20
id <- which((R20_data$MOR < R_pf[1]) == F)
R100_data <- c(R100_data,R20_data[id,]$MOR)
R20_data <- R20_data[-id,]

# R40
id <- which((R40_data$MOR < R_pf[2]) == F)


# R60
id <- which((R60_data$MOR < R_pf[3]) == F)
R100_data <- c(R100_data,R60_data[id,]$MOR)
R60_data <- R60_data[-id,]

# T20
id <- which((T20_data$UTS < T_pf[1]) == F)
T100_data <- c(T100_data,T20_data[id,]$UTS)

T20_data <- T20_data[-id,]

# T40
id <- which((T40_data$UTS < T_pf[2]) == F)
T100_data <- c(T100_data,T40_data[id,]$UTS)
T40_data <- T40_data[-id,]

# T60
id <- which((T60_data$UTS < T_pf[3]) == F)
T100_data <- c(T100_data,T60_data[id,]$UTS)
T60_data <- T60_data[-id,]


##---- Convert the data so stan can use it -----
R20_data <- cbind(R20_data$MOR,R20_data$UTS,R20_data$Broken)
R40_data <- cbind(R40_data$MOR,R40_data$UTS,R40_data$Broken)
R60_data <- cbind(R60_data$MOR,R60_data$UTS,R60_data$Broken)
T20_data <- cbind(T20_data$UTS,T20_data$MOR,T20_data$Broken)
T40_data <- cbind(T40_data$UTS,T40_data$MOR,T40_data$Broken)
T60_data <- cbind(T60_data$UTS,T60_data$MOR,T60_data$Broken)



##------data preprocessing completed -----

## single simulation 
dmg_mod <- stan_model("newmodel.stan")
options(mc.cores = parallel::detectCores())

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha = rep(1,6) )
}

data_stan <- list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                  N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                  N_x = length(T100_data),N_y = length(R100_data),
                  X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                  X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                  t_x = R100_data,t_y = T100_data,
                  l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                  l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3],prop = .89)
set.seed(2020)
real_fit <- sampling(object = dmg_mod,
                     data = data_stan,init = init_dmg)

print(real_fit)

-2482.44
##------ lp vs thresh -------------
dmg_mod <- stan_model("newmodel.stan")

init_dmg <- function() {
  list(mu = c(35,8), sigma = c(10,1), rho = .5, alpha = rep(1,6) )
}
set.seed(2020)


prop_seq <- seq(.7,.99,.01)
lpthresh <- matrix(NA, nrow = length(prop_seq),ncol = 2)
lpthresh[,1] <- prop_seq
for (ii in 1:length(prop_seq)){
  
  show (prop_seq[ii])
data_stan <- list(N_R20 = nrow(R20_data),N_R40 = nrow(R40_data),N_R60 = nrow(R60_data),
                  N_T20 = nrow(T20_data),N_T40 = nrow(T40_data),N_T60 = nrow(T60_data),
                  N_x = length(T100_data),N_y = length(R100_data),
                  X_R20 = R20_data,X_R40 = R40_data,X_R60 = R60_data,
                  X_T20 = T20_data,X_T40 = T40_data,X_T60 = T60_data,
                  t_x = R100_data,t_y = T100_data,
                  l_R20=R_pf[1],l_R40=R_pf[2],l_R60=R_pf[3],
                  l_T20=T_pf[1],l_T40=T_pf[2],l_T60=T_pf[3],prop = prop_seq[ii])

options(mc.cores = parallel::detectCores())
real_fit <- sampling(object = dmg_mod,
                         data = data_stan,init = init_dmg,refresh = 0, iter = 5000)
lpthresh[ii,2] <-  summary(real_fit )[1]$summary[12,1]
}

# saveRDS(lpthresh, file = 'lpthresh.rds')
lpthresh <- readRDS("lpthresh.rds")


realfit <- readRDS("realfit.rds")
print(realfit)

traceplot(realfit, pars = 'alpha')
traceplot(realfit, pars = c('mu','sigma','rho'))



plot(lpthresh[,1], lpthresh[,2] , type = 'l', xlab = 'threshold', ylab = 'lp')


lpthresh[which.max(lpthresh[,2]),]
lpthresh[6,]

real_fit <- readRDS("real_fit.rds")


simlp <- readRDS("simlp.rds")
colnames(simlp) <- c("threshold", "lp")
simlp[5,]
simlp[6,2] <- -8712.58

plot(simlp[,1], simlp[,2],type = 'l', xlab = 'threshold', ylab = 'logpost')


colnames(lpthresh) <- c("threshold", "lp")



dmg_fit <- readRDS("real_fit_thresh.rds")
print(dmg_fit)
