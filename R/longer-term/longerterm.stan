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

// The input data is a vector 'y' of length 'N'.
functions{
    real Dmg_lpdf(real x, real l,real mu, real sigma, real alpha,real thresh) {
    real loglik;
    if (x/alpha < l/thresh){
      loglik = -log(alpha) + normal_lpdf(x/alpha|mu,sigma);
    }
    if (x/alpha > l/thresh && x > l/thresh){
      loglik =  normal_lpdf(x|mu,sigma);
    }
    if (x/alpha > l/thresh && x < l/thresh){
      loglik =  -10000;
    }
    return loglik;
  }
}


data {
  int<lower=0> N_R20;
  int<lower=0> N_R100;
  int<lower=0> N_R20R100;
  real X_R20[N_R20];
  real X_R100[N_R100];
  real X_R20R100[N_R20R100];
  real l;
  real thresh;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> alpha;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
    for (i in 1:N_R20R100){
    target += Dmg_lpdf(X_R20R100[i]|l,mu,sigma,alpha,thresh);
  }
    for (i in 1:N_R20){
    target += normal_lpdf(X_R20[i]|mu,sigma);
    }
    for (i in 1:N_R100){
    target += normal_lpdf(X_R100[i]|mu,sigma);
    }
}

