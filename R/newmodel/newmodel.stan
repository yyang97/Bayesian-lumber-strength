//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {
  
  
 real int_R(real ystar, real[] mu, real[] sigma, real rho,real l){
   real a_l;
   real int_value;
     a_l = (l - mu[1]-rho*(sigma[1]/sigma[2])*(ystar-mu[2]))/(sigma[1]*sqrt(1-rho^2));
  int_value = exp(normal_lpdf(ystar | mu[2],sigma[2]))*(1 - normal_cdf(a_l,0,1));
  return(int_value);
 }
 
 
  real int_T(real xstar, real[] mu, real[] sigma, real rho,real l){
   real a_l;
   real int_value;
     a_l = (l - mu[2]-rho*(sigma[2]/sigma[1])*(xstar-mu[1]))/(sigma[2]*sqrt(1-rho^2));
  int_value = exp(normal_lpdf(xstar | mu[1],sigma[1]))*(1 - normal_cdf(a_l,0,1));
  return(int_value);
 }
 
   real DmgR_lpdf(real[] x, real l,real[] mu, real[] sigma,real rho, real alpha,real thresh) {
    real loglik;
    real a_l;
    if(x[3] == 1){
      //loglik = (-(x[1] - mu[1])^2/(2*sigma[1]^2))-log(sigma[1]);
      loglik = normal_lpdf(x[1]|mu[1],sigma[1]);
    }
    else{
      loglik = log(1/alpha * (int_R(x[2]/alpha,mu,sigma, rho,l) - int_R(x[2]/alpha,mu,sigma, rho,l/thresh)) + int_R(x[2],mu,sigma, rho,l/thresh));
    }
    return loglik;
  }


   real DmgT_lpdf(real[] x, real l,real[] mu, real[] sigma,real rho, real alpha,real thresh) {
    real loglik;
    real a_l;
    if(x[3] == 1){
      //loglik = (-(x[1] - mu[1])^2/(2*sigma[1]^2))-log(sigma[1]);
      loglik = normal_lpdf(x[1]|mu[2],sigma[2]);
    }
    else{
      loglik = log(1/alpha * (int_T(x[2]/alpha,mu,sigma, rho,l) - int_T(x[2]/alpha,mu,sigma, rho,l/thresh)) + int_T(x[2],mu,sigma, rho,l/thresh));
    }
    return loglik;
  }
  
  
}


// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N_R20;
  int<lower=0> N_R40;
  int<lower=0> N_R60;
  int<lower=0> N_T20;
  int<lower=0> N_T40;
  int<lower=0> N_T60;
  int<lower = 0> N_y;
  int<lower = 0> N_x;
  real X_R20[N_R20,3];
  real X_R40[N_R40,3];
  real X_R60[N_R60,3];
  real X_T20[N_T20,3];
  real X_T40[N_T40,3];
  real X_T60[N_T60,3];
  real t_y[N_y];
  real t_x[N_x];
  real l_R20;
  real l_R40;
  real l_R60;
  real l_T20;
  real l_T40;
  real l_T60;
  real thresh;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu[2];
  real<lower=0> sigma[2];
  real<lower=-1,upper = 1> rho;
  real<lower=0, upper = 5> alpha[6];
}

//generated quantities{
 //vector[N_R20] log_lik;
  //for (n in 1:N_R20){
  //  log_lik[n] = DmgR_lpdf(X_R20[n,]|l_R20,mu,sigma,rho,alpha_R20);
  //}
//}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
    for (i in 1:N_R20){
    target += DmgR_lpdf(X_R20[i,]|l_R20,mu,sigma,rho,alpha[1],thresh);
  }
    for (i in 1:N_R40){
    target += DmgR_lpdf(X_R40[i,]|l_R40,mu,sigma,rho,alpha[2],thresh);
  }
    for (i in 1:N_R60){
    target += DmgR_lpdf(X_R60[i,]|l_R60,mu,sigma,rho,alpha[3],thresh);
  }
    for (i in 1:N_T20){
    target += DmgT_lpdf(X_T20[i,]|l_T20,mu,sigma,rho,alpha[4],thresh);
  }
    for (i in 1:N_T40){
    target += DmgT_lpdf(X_T40[i,]|l_T40,mu,sigma,rho,alpha[5],thresh);
  }
    for (i in 1:N_T60){
    target += DmgT_lpdf(X_T60[i,]|l_T60,mu,sigma,rho,alpha[6],thresh);
  }
    for (i in 1:N_y){
    target += normal_lpdf(t_y[i]|mu[2],sigma[2]);
    }
    for (i in 1:N_x){
    target += normal_lpdf(t_x[i]|mu[1],sigma[1]);
    }    
}

