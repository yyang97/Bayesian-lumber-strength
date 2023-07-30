##------- source python function--------
library(reticulate)
use_python("/Users/yyang/Library/r-miniconda/envs/r-reticulate/bin/python")
#use_python("/Users/yunfeng.yang/Library/r-miniconda/envs/r-reticulate/bin/python")
#py_config()
source_python("gkdist.py")


#'@descrption the transformation function to transform the r.v.from the std normal to the r.v. from the g-and-k distribution
#'@param z the std normal r.v. ~ N(0,1)
#'@param theta the theta of g-and-k distribution, theta = c(A, B, g, k)
#'@return the transformed g-k r.v. 
gktrans_r <- function(z,theta){
  A <- theta[1]
  B <- theta[2]
  g <- theta[3]
  k <- theta[4]
  expterm <- exp(-g*z)
  A + B*(1 + .8*(1 - expterm)/(1 + expterm))*(1 + z^2)^k*z
  
}

#'@param N is the sample size 
#'@param theta is the theta of g-and-k distribution, theta = c(A, B, g, k)
#'the defaulted theta (true value) is c(3,1,2,.5)
#'@return the generated g-k r.v. with sample size of N and parameter of theta = c(A, B, g, k)
gk_generate <- function(N, theta){
  # generate the std norm r.v.
  z <- rnorm(N)
  # transform the std norm r.v. to the g-k dist r.v.
  gktrans_r(z = z,theta)
}

  
#'@description The G function (moment condition) in the empirical likelihood
#'@param yobs_sumstat the vector of size `4`, the summary statistics observations from the g-k dist: (mean, 25% quantile, 50% quantile,75% quantile)
#'@param z the matrix of size `N x m`, the basic r.v. from the std norm N(0,1)
#'@param z_quant the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#'@param theta the theta of g-k distribution, theta = c(A, B, g, k)
#'@return the G matrix of size `m x 4`
#'G[,1] the mean summary stat, is mean(y_obs) - colmean(simulated g-k r.v.)
#'G[,2] the 25% quantile summary stat, is ...
#'G[,3] the 50% quantile summary stat, is ...
#'G[,4] the 75% quantile summary stat, is ...

Gfun_r <- function(yobs_sumstat, z, z_quant, theta){
  # transform the basic std norm r.v. to the gk dist r.v.
  x <- gktrans_r(z, theta)
  
  G <- as.matrix(cbind(colMeans(x) - yobs_sumstat[1],
                       gktrans_r(z_quant[,1],theta) - yobs_sumstat[2],
                       gktrans_r(z_quant[,2],theta) - yobs_sumstat[3],
                       gktrans_r(z_quant[,3],theta) - yobs_sumstat[4]))
}
  

#'@description The logel of the g-and-k distribution
#'@param yobs_sumstat the vector of size `4`, the summary statistics observations from the g-k dist: (mean, 25% quantile, 50% quantile,75% quantile)
#'@param z the matrix of size `N x m`, base r.v. used to simulate g-k r.v., of 
#'@param z_quant the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#'@param theta the theta of g-k distribution, theta = c(A, B, g, k)
#'@param supp_adj The logel param: If TRUE, turn on the support correction. If FALSE, turn off it. 
#'@return the logel of the gk example
gklogel <- function(yobs_sumstat, z, z_quant, theta, supp_adj = TRUE){
  G <- Gfun_r(yobs_sumstat, z, z_quant, theta)
  flexEL::logEL(G = G, supp_adj = supp_adj)
}





#'@description The gradient of logel of the g-and-k distribution, dlogel/dtheta
#'@param yobs_sumstat the vector of size `4`, the summary statistics observations from the g-k dist: (mean, 25% quantile, 50% quantile,75% quantile)
#'@param z the matrix of size `N x m`, base r.v. used to simulate g-k r.v., of 
#'@param z_quant the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#'@param theta the theta of g-k distribution, theta = c(A, B, g, k)
#'@param supp_adj The logel param: If TRUE, turn on the support correction. If FALSE, turn off it. 
#'@return The gradient of logel of the gkdist, dlogel/dtheta
gklogel_dldt <- function(yobs_sumstat, z, z_quant, theta, supp_adj = TRUE){
  # calculate the G matrix
  G_r <- Gfun_r(yobs_sumstat, z, z_quant, theta)
  # get the gradient dlogel/dG
  dldG <- flexEL::logEL(G = G_r, supp_adj = supp_adj,grad = TRUE)$grad
  # get the gradient dG/dtheta
  dGdt <- Gfun_grad(yobs_sumstat, z, z_quant, theta)
  dldt_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:nrow(dldG)) {
    # dGdb 40 x 4 x 4
    # dldG 40 x 4
   # dldt_mat[ii,] <- dGdt[ii,,] %*% dldG[ii,]
    dldt_mat[ii,] <- dldG[ii,] %*% dGdt[ii,,]
  }
  dldt <- colSums(dldt_mat)
  dldt
}



#'@description The negative logel of the g-and-k distribution
#'@param yobs_sumstat the vector of size `4`, the summary statistics observations from the g-k dist: (mean, 25% quantile, 50% quantile,75% quantile)
#'@param z the matrix of size `N x m`, base r.v. used to simulate g-k r.v., of 
#'@param z_quant the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#'@param theta the theta of g-k distribution, theta = c(A, B, g, k)
#'@param supp_adj The logel param: If TRUE, turn on the support correction. If FALSE, turn off it. 
#'@return the negative logel of the gk example, sent to be minimized
neggklogel <- function(theta, supp_adj = TRUE){
  G <- Gfun_r(yobs_sumstat, z, z_quant, theta)
  -flexEL::logEL(G = G, supp_adj = supp_adj )
}





#'@description The negative gradient of logel of the g-and-k distribution, dlogel/dtheta
#'@param yobs_sumstat the vector of size `4`, the summary statistics observations from the g-k dist: (mean, 25% quantile, 50% quantile,75% quantile)
#'@param z the matrix of size `N x m`, base r.v. used to simulate g-k r.v., of 
#'@param z_quant the matrix of size ` m x 3`,  25%, 50%, 75% quantiles of z
#'@param theta the theta of g-k distribution, theta = c(A, B, g, k)
#'@param supp_adj The logel param: If TRUE, turn on the support correction. If FALSE, turn off it. 
#'@return The gradient of negative logel of the gkdist, dlogel/dtheta
neggklogel_dldt <- function(theta, supp_adj = TRUE){
  # calculate the G matrix
  G_r <- Gfun_r(yobs_sumstat, z, z_quant, theta)
  # get the gradient dlogel/dG
  dldG <- flexEL::logEL(G = G_r, supp_adj = supp_adj,grad = TRUE)$grad
  # get the gradient dG/dtheta
  dGdt <- Gfun_grad(yobs_sumstat, z, z_quant, theta)
  dldt_mat <- matrix(NA, nrow = nrow(dldG), ncol = ncol(dldG))
  for (ii in 1:nrow(dldG)) {
    # dGdb 40 x 4 x 4
    # dldG 40 x 4
    # dldt_mat[ii,] <- dGdt[ii,,] %*% dldG[ii,]
    dldt_mat[ii,] <- dldG[ii,] %*% dGdt[ii,,]
  }
  dldt <- colSums(dldt_mat)
  -dldt
}



