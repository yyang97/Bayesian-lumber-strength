#source("dmg_func.R")
## plot 
#set.seed(8800)
# the orignal sample size 
N <- 300
mu <- 48
sigma <- 19
#y <- rnorm(N,mean = mu, sd = sigma)
l <-  32
alpha <- 0.85
c <- 0.69
s <- 10



y <- rnorm(N, mean = mu, sd = sigma )
g23 <- y[y>l]

g23_star <- ifelse(l > c*g23,alpha*g23, g23)
#g23_star <- dmg_model(g23,alpha,l,c,s)
# if damaged, then set 1
dmg_cond <- ifelse(g23 > g23_star, 1,0)
g3_id <- g23_star > l

#cbind(g23_star,g23_smooth,dmg_ind)

# g3
g23_data <- cbind(g23,g23_star,dmg_cond)
g3_data <- g23_data[g3_id,]
# g3_data
# 
# l/c
mean(g23)
mean(g3_data[,1])
mean(g3_data[,2])

mu_truc <- truncnorm::etruncnorm(a = l/c, b = Inf, mean = mu, sd = sigma)
var_truc <- truncnorm::vtruncnorm(a = l/c, b = Inf, mean = mu, sd = sigma)

sample_mean <- mu_truc
sample_var <- 1/length(g3_data[,2])*var_truc
pnorm(mean(g3_data[,2]), mean = sample_mean, sd = sqrt(sample_var))



# test alpha ==1
mu_truc <- truncnorm::etruncnorm(a = l, b = Inf, mean = mu, sd = sigma)
var_truc <- truncnorm::vtruncnorm(a = l, b = Inf, mean = mu, sd = sigma)

sample_mean <- mu_truc
sample_var <- 1/length(g3_data[,2])*var_truc
pnorm(mean(g3_data[,2]), mean = sample_mean, sd = sqrt(sample_var))



#mean(g3_data[,2])
