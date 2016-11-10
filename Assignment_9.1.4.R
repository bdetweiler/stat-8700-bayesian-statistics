#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_9.1.4.R
# 2016-11-04
#########################################
library(rstan)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_9_1_4"

years <- c(1976:1985)

fatal.accidents <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)
death.rate <- c(0.19, 0.12, 0.15, 0.16, 0.14, 0.06, 0.13, 0.13, 0.03, 0.15)

passenger.deaths.data <- c(734, 516, 754, 877, 814, 362, 764, 809, 223, 1066)
N <- length(passenger.deaths.data)

modelString = 
  "
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
"

stan.file <- paste0(fileName, ".stan")
writeLines(modelString, con=stan.file)

model.passenger.deaths.data <- list(y = passenger.deaths.data, N = N)

fit <- stan(file = stan.file, 
            data = model.passenger.deaths.data, 
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
x <- seq(200, 1200, length = 1000)

y <- dgamma(x, shape = sum(data) + 1, rate = N + 1)
lines(x, y, type="l", lwd=1, col="red")

CI <- quantile(fit.mcmc.M[, "theta"], c(0.025, 0.975))
CI[[1]]
CI[[2]]

abline(v=CI[[1]], col="red")
abline(v=CI[[2]], col="red")
abline(v=y.bar, col="red")