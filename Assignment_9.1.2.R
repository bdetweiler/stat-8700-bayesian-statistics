#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_9.1.2.R
# 2016-11-04
#########################################

library(rstan)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_9_1_2"

# 2. Suppose we consider a Normal population with a variance of 16, and we collect 15 
#    observations from this population with values: 
#         26.8, 26.3, 28.3, 28.5, 26.3, 31.9, 28.5, 27.2, 20.9, 27.5, 28.0, 18.6, 22.3, 25.0, 31.5.
#     (a) If we choose a Normal(20, 25) prior, 
#         Use R to find the posterior distribution for the population mean.
var <- 16
data <- c(26.8, 26.3, 28.3, 28.5, 26.3, 31.9, 28.5, 27.2, 20.9, 27.5, 28.0, 18.6, 22.3, 25.0, 31.5)
y.bar <- mean(data)
N <- length(data)
sigma_2 <- var

# Write the JAGS model file (only need to do this part once)
modelString = 
"
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

"


stan.file <- paste0(fileName, ".stan")
writeLines(modelString, con=stan.file)

#inits <- list(list(theta=500), list(theta=800), list(theta=1000), list(theta=1200))

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
x <- seq(20, 40, length = 1000)
y <- dnorm(x, mean = y.bar, sd = sqrt(var / N))
lines(x, y, type="l", lwd=1, col="red")

CI <- quantile(fit.mcmc.M[, "theta"], c(0.025, 0.975))
CI[[1]]
CI[[2]]

abline(v=CI[[1]], col="red")
abline(v=CI[[2]], col="red")
abline(v=y.bar, col="red")