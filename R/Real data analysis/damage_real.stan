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
  // the _lpdf syntax has a special meaning (see below)
  real DmgR_lpdf(real[] x, real l,real[] mu, real[] sigma,real rho, real alpha) {
    real loglik;
    real a_l;
    real x_real;
    real derivative;
    if(x[3] == 1){
      //loglik = (-(x[1] - mu[1])^2/(2*sigma[1]^2))-log(sigma[1]);
      loglik = normal_lpdf(x[1]|mu[1],sigma[1]);
    }
    else{
      x_real = (x[2] + sqrt(x[2]^2 + 4*alpha))*0.5;
      derivative = 0.5+x[2]/(2*sqrt(x[2]^2 + 4* alpha));
      a_l = (l - mu[1]-rho*(sigma[1]/sigma[2])*(x_real-mu[2]))/(sigma[1]*sqrt(1-rho^2));
      loglik = normal_lpdf(x_real|mu[2],sigma[2]) + normal_lccdf(a_l | 0, 1) + log(derivative);
    }
    return loglik;
  }



  // the _lpdf syntax has a special meaning (see below)
  real DmgT_lpdf(real[] x, real l,real[] mu, real[] sigma,real rho, real alpha) {
    real loglik;
    real a_l;
    real x_real;
    real derivative;
    if(x[3] == 1){
      loglik = normal_lpdf(x[1]|mu[2],sigma[2]);
    }
    else{
      x_real = (x[2] + sqrt(x[2]^2 + 4*alpha))*0.5;
      derivative = 0.5+x[2]/(2*sqrt(x[2]^2 + 4* alpha));
      a_l = (l - mu[2]-rho*(sigma[2]/sigma[1])*(x_real-mu[1]))/(sigma[2]*sqrt(1-rho^2));
      loglik = normal_lpdf(x_real|mu[1],sigma[1]) + normal_lccdf(a_l | 0, 1) + log(derivative);
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
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu[2];
  real<lower=0> sigma[2];
  real<lower=-1,upper = 1> rho;
  real alpha_R20;
  real alpha_R40;
  real alpha_R60;
  real alpha_T20;
  real alpha_T40;
  real alpha_T60;
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
    target += DmgR_lpdf(X_R20[i,]|l_R20,mu,sigma,rho,alpha_R20);
  }
    for (i in 1:N_R40){
    target += DmgR_lpdf(X_R40[i,]|l_R40,mu,sigma,rho,alpha_R40);
  }
    for (i in 1:N_R60){
    target += DmgR_lpdf(X_R60[i,]|l_R60,mu,sigma,rho,alpha_R60);
  }
    for (i in 1:N_T20){
    target += DmgT_lpdf(X_T20[i,]|l_T20,mu,sigma,rho,alpha_T20);
  }
    for (i in 1:N_T40){
    target += DmgT_lpdf(X_T40[i,]|l_T40,mu,sigma,rho,alpha_T40);
  }
    for (i in 1:N_T60){
    target += DmgT_lpdf(X_T60[i,]|l_T60,mu,sigma,rho,alpha_T60);
  }
    for (i in 1:N_y){
    target += normal_lpdf(t_y[i]|mu[2],sigma[2]);
    }
    for (i in 1:N_x){
    target += normal_lpdf(t_x[i]|mu[1],sigma[1]);
    }    
}

// 
// generated quantities {
//   vector[N_R20+N_R40+N_R60+N_T20+N_T40+N_T60+N_x+N_y] log_lik;
//   for (n in 1:N_R20) {
//     log_lik[n] = 
//     DmgR_lpdf(X_R20[n,]|l_R20,mu,sigma,rho,alpha_R20);
//   }
//   for (n in 1:N_R40) {
//     log_lik[n+N_R20] = 
//     DmgR_lpdf(X_R40[n,]|l_R40,mu,sigma,rho,alpha_R40);
//   }
//   for (n in 1:N_R60) {
//     log_lik[n+N_R20+N_R40] =
//     DmgR_lpdf(X_R60[n,]|l_R60,mu,sigma,rho,alpha_R60);
//   }
//   for (n in 1:N_T20) {
//     log_lik[n+N_R20+N_R40+N_R60] =
//     DmgT_lpdf(X_T20[n,]|l_T20,mu,sigma,rho,alpha_T20);
//   }
//   for (n in 1:N_T40) {
//     log_lik[n+N_R20+N_R40+N_R60+N_T20] =
//     DmgT_lpdf(X_T40[n,]|l_T40,mu,sigma,rho,alpha_T40);
//   }
//   for (n in 1:N_T60) {
//     log_lik[n+N_R20+N_R40+N_R60+N_T20+N_T40] =
//     DmgT_lpdf(X_T60[n,]|l_T60,mu,sigma,rho,alpha_T60);
//   }
//   for (n in 1:N_x) {
//     log_lik[n+N_R20+N_R40+N_R60+N_T20+N_T40+N_T60] =
//     normal_lpdf(t_x[n]|mu[1],sigma[1]);
//   }
//   for (n in 1:N_y) {
//     log_lik[n+N_R20+N_R40+N_R60+N_T20+N_T40+N_T60+N_x] =
//     normal_lpdf(t_y[n]|mu[2],sigma[2]);
//   }
// }


