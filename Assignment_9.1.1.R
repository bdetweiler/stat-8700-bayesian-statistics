#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_9.1.1.R
# 2016-11-04
#########################################
library(rstan)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_9_1_1"
# 1. Suppose we have a population described by a Normal Distribution with 
#    known variance $\sigma^2 = 1600$ and unknown mean $\mu$. 
#    4 observations are collected from the population and the corresponding values were: 
#    940, 1040, 910, and 990.
data <- c(940, 1040, 910, 990)
N <- length(data)
y.bar <- mean(data)

# Write the JAGS model file (only need to do this part once)
modelString = 
"
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

"

stan.file <- paste0(fileName, ".stan")
writeLines(modelString, con=stan.file)

model.data <- list(y = data, N = N)
fit <- stan(file = stan.file, 
            data = model.data, 
            #init = inits, 
            iter = 120000, 
            warmup = 2000, 
            chains = 4, 
            cores = 4, 
            open_progress = F)

print(fit)

fit.mcmc <- As.mcmc.list(fit)
diagMCMC(fit.mcmc)

fit.mcmc.M <- as.matrix(fit)

hist(fit.mcmc.M[,"theta"], breaks=80, freq=F)
x <- seq(0, 1200, length = 1000)
y <- dnorm(x, mean = y.bar, sd = sqrt(1600 / N))
lines(x, y, type="l", lwd=1, col="red")

CI <- quantile(fit.mcmc.M[, "theta"], c(0.025, 0.975))
CI[[1]]
CI[[2]]

abline(v=CI[[1]], col="red")
abline(v=CI[[2]], col="red")
abline(v=y.bar, col="red")