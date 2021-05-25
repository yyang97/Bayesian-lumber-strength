data {
int n; // number of observation
real y[n]; // the observed crimerate
real homeincome[n,2]; // distance to haz site
real<lower = 0> wmat[n, n]; // weight matrix w_ij
}


parameters {
real u[n];
real beta[2];
real<lower=0> sigma2_u;
real<lower=0> sigma2_v;
}


model {
  sigma2_u ~ inv_gamma(0.001, 0.001);
  sigma2_v ~ inv_gamma(0.001, 0.001);
  beta ~ uniform(-10,10);
  target += -0.5*n*log(sigma2_u);
  for (i in 1:n) {
    for (j in 1:n) {
      target += -(u[i] - u[j])^2 * wmat[i,j] / (2 * sigma2_u);
    }
  }
  for (i in 1:n) {
    y[i] ~ normal(beta[1]*homeincome[i,1] + beta[2]*homeincome[i,2]+
    u[i], sqrt(sigma2_v));
  }
}
