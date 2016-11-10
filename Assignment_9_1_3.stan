
data {
  int<lower=0> N; 
  int<lower=0> y[N];
}
parameters {
  real<lower=0> theta;
}
model {
	theta ~ gamma(1, 1);   // Prior
  y ~ poisson(theta);   // Likelihood
}

