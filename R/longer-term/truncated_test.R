# x ~ N(2,3^2)
set.seed(8888)
N <- 1000
x <- rnorm(N,mean = 2,sd = 3)
x_obs <- x[x >1 ]
x_obs

# test_statistics
test <- (mean(x_obs) -2*(1/sqrt(length(x_obs))))/3

test

mean(x_obs)

library(truncnorm)
truncnorm::ptruncnorm(test, a=(1), b=Inf, mean = 0, sd = 1)


mean_xobs <- etruncnorm(a=1, b=Inf, mean=2, sd=3)
sd_xobs <-sqrt(vtruncnorm(a=1, b=Inf, mean=2, sd=3))


sd_xobs*mean(x_obs) + mean_xobs*sqrt(length(x_obs))
