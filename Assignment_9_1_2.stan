
data {
  int<lower=0> N; 
  real<lower=0> y[N];
}
parameters {
  real<lower=0, upper=100> theta;
}
model {
	theta ~ normal(20, 5);   // Prior
  y ~ normal(theta, 4);    // Likelihood
}


