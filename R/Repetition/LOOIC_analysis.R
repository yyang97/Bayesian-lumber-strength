setwd("~/Downloads/github/Bayesian-lumber-strength/R/Repetition")
##---  Analysis part ---##
allcombs <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
colnames(allcombs) <- c('alpha_R20','alpha_R40','alpha_R60',
                        'alpha_T20','alpha_T40','alpha_T60')
allcombs

find_ci <- function(x, digits = 1) {
  round(c(lower = unname(x[1] - 2*x[2]),
          upper = unname(x[1] + 2*x[2])),
        digits = digits)
}
## --- no dmg reptition
load('loo_value.rda')
dim(loo_value)
loo_value[,1]

# real data size: 87
# repetition size: 2*87
# chosen model by LOOIC
apply(loo_value[1:55,],1,which.min)
# percent of choosing correct model 
mean(apply(loo_value[1:55,],1,which.min) == 1)
# check some wrong results
# repetition 3: close result
loo_value[3,9]
loo_value[3,1]
# repetition 4: close result
loo_value[4,33]
loo_value[4,1]

## partial 
# alpha_R <- 10* c(0,1.83,0)
# alpha_T <- 10*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model 
load('loo_chosen_partial.rda')
# the correct model
allcombs[51,]
# check how many repetitions are completed
loo_value[,1]
# the models LOOIC chosen
apply(loo_value[1:88,],1,which.min)
# the percent of LOOIC choosing correct model
mean(apply(loo_value[1:88,],1,which.min) == 51)
# check some wrong results
# check repetition 3: large differences 
loo_value[3,59]
loo_value[3,51]

allcombs[59,]

# check repetition 2: large differences
loo_value[2,56]
loo_value[2,51]

allcombs[56,]
##-------only R40 ---------
# alpha_R <- 1* c(0,1.3,0)
# alpha_T <- 0*c(0,104.09,65.95)
# correct model: c(0,1,0,0,0,0)
load('loo_chosen_onlyR40.rda')
# correct model: the 3rd model 
allcombs[3,]
# check how many repetitions are completed
loo_value[,1]

# check the model that LOOIC chose
apply(loo_value,1,which.min)
# the percent of LOOIC choosing correct model
mean(apply(loo_value,1,which.min) == 3)
# repetition 6: small difference
loo_value[6,19]
loo_value[6,3]
# repetition 2: small difference
loo_value[2,11]
loo_value[2,3]


# model with full damage 
load('loo_chosen_full.rda')
# alpha_R <- 10* c(0.67,1.83,1.29)
# alpha_T <- 10*c(85.23,104.09,65.95)
# the correct model: model 64 c(1,1,1,1,1,1)
allcombs[64,]

apply(loo_value[1:4,],1,which.min)
hist(loo_value[4,])
apply(loo_value,1,which.min)
loo_compare()

## 

## partial 
# alpha_R <- 10* c(0,1.83,0)
# alpha_T <- 10*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model 



# try 64 models on real data
# try 2* 5* for partial effects, try 2-5 repetitions 


##-------only R40 ---------
# alpha_R <- 1* c(0,1.3,0)
# alpha_T <- 0*c(0,104.09,65.95)
# correct model: the 3rd model 
load("loo_fit_onlyR40.rda")
loo_value[[1]]
require(loo)
length(loo_value)

loo_list <- paste0("loo_value[[6]][[",1:64,collapse = "]],")
loo_list <- paste0(loo_list,"]]")
loo_list



#loo_value[[1]][[11]]$estimates[3,1]
#loo_value[[1]][[31]]$estimates[3,1]



loocomp <- loo_compare(loo_value[[6]][[1]],loo_value[[6]][[2]],loo_value[[6]][[3]],loo_value[[6]][[4]],loo_value[[6]][[5]],loo_value[[6]][[6]],loo_value[[6]][[7]],loo_value[[6]][[8]],loo_value[[6]][[9]],loo_value[[6]][[10]],loo_value[[6]][[11]],loo_value[[6]][[12]],loo_value[[6]][[13]],loo_value[[6]][[14]],loo_value[[6]][[15]],loo_value[[6]][[16]],loo_value[[6]][[17]],loo_value[[6]][[18]],loo_value[[6]][[19]],loo_value[[6]][[20]],loo_value[[6]][[21]],loo_value[[6]][[22]],loo_value[[6]][[23]],loo_value[[6]][[24]],loo_value[[6]][[25]],loo_value[[6]][[26]],loo_value[[6]][[27]],loo_value[[6]][[28]],loo_value[[6]][[29]],loo_value[[6]][[30]],loo_value[[6]][[31]],loo_value[[6]][[32]],loo_value[[6]][[33]],loo_value[[6]][[34]],loo_value[[6]][[35]],loo_value[[6]][[36]],loo_value[[6]][[37]],loo_value[[6]][[38]],loo_value[[6]][[39]],loo_value[[6]][[40]],loo_value[[6]][[41]],loo_value[[6]][[42]],loo_value[[6]][[43]],loo_value[[6]][[44]],loo_value[[6]][[45]],loo_value[[6]][[46]],loo_value[[6]][[47]],loo_value[[6]][[48]],loo_value[[6]][[49]],loo_value[[6]][[50]],loo_value[[6]][[51]],loo_value[[6]][[52]],loo_value[[6]][[53]],loo_value[[6]][[54]],loo_value[[6]][[55]],loo_value[[6]][[56]],loo_value[[6]][[57]],loo_value[[6]][[58]],loo_value[[6]][[59]],loo_value[[6]][[60]],loo_value[[6]][[61]],loo_value[[6]][[62]],loo_value[[6]][[63]],loo_value[[6]][[64]])
loocomp
rownames(loocomp)[1]
find_ci <- function(x, digits = 1) {
   round(c(lower = unname(x[1] - 2*x[2]),
             upper = unname(x[1] + 2*x[2])),
           digits = digits)
   }
t(apply(loocomp[-1, c("elpd_diff", "se_diff")], 1, find_ci))






# loo_list <- paste0("loo_fit[[",1:64,collapse = "]],")
# loo_list <- paste0(loo_list,"]]")
# loocomp <- loo_compare(loo_fit[[1]],loo_fit[[2]],loo_fit[[3]],loo_fit[[4]],loo_fit[[5]],loo_fit[[6]],loo_fit[[7]],loo_fit[[8]],loo_fit[[9]],loo_fit[[10]],loo_fit[[11]],loo_fit[[12]],loo_fit[[13]],loo_fit[[14]],loo_fit[[15]],loo_fit[[16]],loo_fit[[17]],loo_fit[[18]],loo_fit[[19]],loo_fit[[20]],loo_fit[[21]],loo_fit[[22]],loo_fit[[23]],loo_fit[[24]],loo_fit[[25]],loo_fit[[26]],loo_fit[[27]],loo_fit[[28]],loo_fit[[29]],loo_fit[[30]],loo_fit[[31]],loo_fit[[32]],loo_fit[[33]],loo_fit[[34]],loo_fit[[35]],loo_fit[[36]],loo_fit[[37]],loo_fit[[38]],loo_fit[[39]],loo_fit[[40]],loo_fit[[41]],loo_fit[[42]],loo_fit[[43]],loo_fit[[44]],loo_fit[[45]],loo_fit[[46]],loo_fit[[47]],loo_fit[[48]],loo_fit[[49]],loo_fit[[50]],loo_fit[[51]],loo_fit[[52]],loo_fit[[53]],loo_fit[[54]],loo_fit[[55]],loo_fit[[56]],loo_fit[[57]],loo_fit[[58]],loo_fit[[59]],loo_fit[[60]],loo_fit[[61]],loo_fit[[62]],loo_fit[[63]],loo_fit[[64]])



## partial
# alpha_R <- 10* c(0,1.83,0)
# alpha_T <- 10*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model 

loo_partial <- readRDS("loo_fit_partial.rds")
loo_partial <- loo_partial[1:10]

sapply(loo_partial, function(x){rownames(x)[1]})

allcombs[59,]
allcombs[56,]
allcombs[51:64,]
loo_partial[[2]]
51 52 55 56 59 60 63 64

t(apply(loo_partial[[2]][-1, c("elpd_diff", "se_diff")], 1, find_ci))
t(apply(loo_partial[[4]][-1, c("elpd_diff", "se_diff")], 1, find_ci))
t(apply(loo_partial[[6]][-1, c("elpd_diff", "se_diff")], 1, find_ci))
t(apply(loo_partial[[10]][-1, c("elpd_diff", "se_diff")], 1, find_ci))


## partial
# alpha_R <- 2* c(0,1.83,0)
# alpha_T <- 2*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model

partial_2times <- readRDS("loo_fit_partial_2times.rds")
partial_2times[[4]]

t(apply(partial_2times[[1]][-1, c("elpd_diff", "se_diff")], 1, find_ci))




## partial
# alpha_R <- 5* c(0,1.83,0)
# alpha_T <- 5*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model

partial_5times <- readRDS("loo_fit_partial_5times.rds")
partial_5times[[3]]

t(apply(partial_5times[[2]][-1, c("elpd_diff", "se_diff")], 1, find_ci))






## partial
# alpha_R <- 3* c(0,1.83,0)
# alpha_T <- 3*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model

partial_3times <- readRDS("loo_fit_partial_3times.rds")
partial_3times[[3]]

t(apply(partial_3times[[3]][-1, c("elpd_diff", "se_diff")], 1, find_ci))

print(partial_3times[[3]],simplify = F)
# partial with rho 

partial_3times <- readRDS("loo_fit_partial_3times.rds")
partial_3times[[3]]

t(apply(partial_3times[[3]][-1, c("elpd_diff", "se_diff")], 1, find_ci))



rho_3times <- readRDS("rho_3times.rds")
modelname <- c()
for(ii in 1:64){
  modelname[ii] <- paste0("model",ii)
}
for(ii in 1:length(rho_3times)){
rownames(rho_3times[[ii]]) <- modelname
}         
         

t(apply(partial_3times[[3]][-1, c("elpd_diff", "se_diff")], 1, find_ci))
rho_3times[[3]][rownames(partial_3times[3][[1]]),]

partial_3times[[1]]







## partial
# alpha_R <- 5* c(0,1.83,0)
# alpha_T <- 5*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model

# partial with rho 

partial_5times <- readRDS("loo_fit_partial_5times.rds")
print(partial_5times[[6]],simplify = F)

t(apply(partial_5times[[6]][-1, c("elpd_diff", "se_diff")], 1, find_ci))



rho_5times <- readRDS("rho_5times.rds")
modelname <- c()
for(ii in 1:64){
  modelname[ii] <- paste0("model",ii)
}
for(ii in 1:length(rho_5times)){
  rownames(rho_5times[[ii]]) <- modelname
}         


t(apply(partial_5times[[6]][-1, c("elpd_diff", "se_diff")], 1, find_ci))

class(partial_5times[[6]])
rho_5times[[6]][rownames(partial_5times[6][[1]]),]

attributes(partial_5times[[6]])




## partial
# alpha_R <- 10* c(0,1.83,0)
# alpha_T <- 10*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model
# rho = 0.8
# partial with rho 

partial_10times <- readRDS("loo_fit_partial_10times.rds")
print(partial_10times[[1]],simplify = F)

t(apply(partial_10times[[1]][-1, c("elpd_diff", "se_diff")], 1, find_ci))



rho_10times <- readRDS("rho_10times.rds")
modelname <- c()
for(ii in 1:64){
  modelname[ii] <- paste0("model",ii)
}
for(ii in 1:length(rho_10times)){
  rownames(rho_10times[[ii]]) <- modelname
}         


t(apply(partial_10times[[2]][-1, c("elpd_diff", "se_diff")], 1, find_ci))

partial_10times[[2]]
rho_10times[[2]][rownames(partial_10times[2][[1]]),]

colnames(partial_10times[[1]])


cbind(partial_10times[[2]][,7],rho_10times[[2]][rownames(partial_10times[2][[1]]),])


matrix(partial_10times[[1]])
## partial
# alpha_R <- 2* c(0,1.83,0)
# alpha_T <- 2*c(0,104.09,65.95)
# true model: c(0,1,0,0,1,1)
# the 51th model

# partial with rho 

partial_2times <- readRDS("loo_fit_partial_2times.rds")
print(partial_2times[[1]],simplify = F)

t(apply(partial_2times[[1]][-1, c("elpd_diff", "se_diff")], 1, find_ci))



rho_2times <- readRDS("rho_2times.rds")
modelname <- c()
for(ii in 1:64){
  modelname[ii] <- paste0("model",ii)
}
for(ii in 1:length(rho_2times)){
  rownames(rho_2times[[ii]]) <- modelname
}         


t(apply(partial_2times[[2]][-1, c("elpd_diff", "se_diff")], 1, find_ci))

partial_2times[[2]]
rho_2times[[2]][rownames(partial_2times[2][[1]]),]


cbind(partial_2times[[2]][,7],rho_2times[[2]][rownames(partial_2times[2][[1]]),])

# given a dataset 
# task 1 : whether there is damage effect or not
# task 2: if there is, what is the estimation of the damage effect
# task 3: what is the estimation of rho



# high damage scenario
# 8 models similar looic
# do they also have similar rho?

# low damage scenario
# credible intervals of rho hat


