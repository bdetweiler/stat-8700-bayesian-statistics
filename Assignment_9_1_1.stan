
data {
  int<lower=0> N; 
  real<lower=0> y[N];
}
parameters {
  real<lower=0> theta;
}
model {
	theta ~ normal(1000, 200); // Prior
  y ~ normal(theta, 40);    // Likelihood
}


