setwd("~/Downloads/github/Bayesian-lumber-strength/R/longer-term")
library(rstan)
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

mean(R100_data)


thresh <-  .71
dmglik <- function(theta){
  # mu <- theta[1]
  # sigma <- theta[2]
  mu <- 47.91
  sigma <- 18.84
  alpha <- theta
  data <- c()
  likindi <- c()
  est_y <-c()
  lik <- 0
  cond <-c()
  for (jj in 1:length(R20R100_data)){
    data[jj] <- R20R100_data[jj]
    if(R20R100_data[jj]/alpha < l/thresh){
      lik <- lik - log(alpha) + dnorm(R20R100_data[jj]/alpha, mean =  mu,sd = sigma, log = TRUE)
      likindi[jj] <-  - log(alpha) + dnorm(R20R100_data[jj]/alpha, mean =  mu,sd = sigma, log = TRUE)
      est_y[jj] <- data[jj]/alpha
      cond[jj] <- "damaged"
    }else if(R20R100_data[jj]/alpha > l/thresh && R20R100_data[jj] > l/thresh){
      lik <- lik + dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)
      likindi[jj] <- dnorm(R20R100_data[jj], mean = mu, sd = sigma,log = TRUE)
      est_y[jj] <- data[jj]
      cond[jj] <- 'nondamaged'
    }
    else if(R20R100_data[jj]/alpha > l/thresh && R20R100_data[jj] < l/thresh){
      lik <- lik-10000
      likindi[jj] <- -10000
      est_y[jj] <- "incompatible"
      cond[jj] <- "incompatible"
      
    }
  }
  # p(l > alpha*Y) has 97 pieces, also l > cy,  c is threshold
  # l > alpha*Y > cy, alpha >c?
  # Given 97 broken pieces, Y < l/alpha, Y < l/c,so l/c > l/alpha and alpha >c.
  # alpha = 0.97, thresh = 0.8, if we have a 0.97*Y <4500 in group 2, what should be Y
  # all the Y should be Y < 4639.175, of course 
  # of course .8*Y < l, Y < 5625
  # and we know Y > l = 4500, 
  # So 4500 < Y < 4639.175 
  # add all of the 3 conditions 
  
  # P(l < y < min(l/c,l/alpha)), group 2
  
   lik <- lik + log(pnorm(min(l/alpha,l/thresh),mean = mu, sd = sigma) - 
                   pnorm(l,mean = mu, sd = sigma))*97
   liktable <- data.frame(data,likindi,est_y,cond)
  # lik <-  lik + sum(dnorm(c(R100_data,R20_data), mean = mu, sd = sigma, log = TRUE))
   return(-1*lik)
  # return(list(liktable = liktable, likvalue = -1*lik))
}
optimize(dmglik, c(0.9,.999))
dmglik(.99)
dmglik(.9646)


dmglik(.98)$liktable
dmglik(.98)$likvalue
pnorm(min(31/.5, 31/.9),mean = 47.91, sd = 18.84) -  pnorm(31,mean = 47.91, sd = 18.84)
# alpha > c
# If alpha < c, then alpha*y < c*y, it means all the pieces are broken after the proof loading. 
# There should be no pieces in the group 3.



alphaseq <- seq(from = 0.1,to = 1, by = 0.01)
likseq <- c()
for (iseq in 1:length(alphaseq) ){
  likseq[iseq] <- dmglik(alphaseq[iseq])
}

plot(alphaseq,likseq,type = 'l')
# If thresh = .9, it means that all the observed data are nondamaged.
# The plot only demontrates the changes of  loglik in the group 2. 
# thresh = 0.8, max loglik = -1038.864, it suggests alpha = 0.965,
# thresh = 0.9, max loglik = -925.6903, it suggests alpha = 0.9
# thresh = 0.85, max loglik = -876.8572, it suggests alpha = .85
# choose the thresh that maximizes the loglik. 
logpost <- readRDS('logpost.rds')

plot(thresh_seq[c(10:30)],logpost[c(10:30)],type = 'l')



dmglik(.1)$likvalue
dmglik(.2)$likvalue

# dmglik (.1)
# dmglik (.05)
dmglik (1)
dmglik (0.98)
optimize(dmglik,c(0,1.2))$minimum

# record the likelihood of every piece, 
# a table of likelihood for each piece 
# estimate y from the observed y^* and evaluate whether y are reasonable.
# make a plot y vs y^*, given alpha and threshold




# l = 4500, suppose threshold = .8, we can choose an alpha = 0.9
# the damaged y*, y^* = y*I(thresh*y > l) + alpha *y*I(thresh*y < l)
# y^* =  y*I(thresh*y > l) +  (1/alpha* y + l/thresh - alpha*l/thresh)*I(thresh*y < l)



yseq <- seq(from = 1, to = 8000, by =1 )
#ystarseq <- (.7*yseq)*ifelse(.8*yseq<4500,1,0) + yseq*ifelse(.8*yseq>=4500,1,0)
ystarseq <- (.7*yseq + 4500/.8 - .7*4500/.8)*ifelse(.8*yseq<4500,1,0) + yseq*ifelse(.8*yseq>=4500,1,0)
#ystarseq <- ((1/0.1)*yseq + 4500/.8 - 1/0.1*4500/.8)*ifelse(.8*yseq<4500,1,0) + yseq*ifelse(.8*yseq>=4500,1,0)

plot(yseq,ystarseq, type = 'l', ylim =c(0,8000))
# the dash line is no damaged plot 
abline(yseq, yseq, lty = 2)
# so we see the discontinuity 
# we cannot observed any data ranging from alpha*thresh and thresh, that's why of the discontinuity
# do we have this discontinuity in the bivariate dataset?
# all the pieces are either broken or non-damaged 
# maintain the linear model and describe how far y is higher than the threshold
# the plot of y and y* should be more smooth, such that there is no contradiction 



yseq <- seq(from = 1, to = 100, by =.1 )
#ystarseq <- (.7*yseq)*ifelse(.8*yseq<32,1,0) + yseq*ifelse(.8*yseq>=32,1,0)
#ystarseq <- (.7*yseq +32/.8 - .7*32/.8)*ifelse(.8*yseq<32,1,0) + yseq*ifelse(.8*yseq>=32,1,0)
ystarseq <- ((1/0.5)*yseq + 32/.8 - 1/0.5*32/.8)*ifelse(.8*yseq<32,1,0) + yseq*ifelse(.8*yseq>=32,1,0)

plot(yseq,ystarseq, type = 'l', ylim =c(0,100))
# the dash line is no damaged plot 
abline(yseq, yseq, lty = 2)

