#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_9.3.0.R
# 2016-11-04
#########################################
#library(rstan)
library(rjags)
library(runjags)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_9_3_0"

# Read in the rat data - You may need to specify the full path to the file.
ratdata <- read.table("rats.txt", header=TRUE)

# Log Posterior (u, v space)
log.post2 <- function(u, v) {
  
  
  deaths <- ratdata$y
  rats <- ratdata$N
  alpha <- exp(u + v) / (1 + exp(u))
  beta <- exp(v) / (1 + exp(u))
  ldens <- 0
  
  # Loop over each of the 71 experiments
  for(i in 1:length(rats)) {
    
    # There is no gamma function in R - have to use Log Gamma, so the density is logged
    # This is why it's addative rather than multiplicative
    # deaths[i] is the same as y_i
    # rats[i] is the same as n_i
    ldens <- (ldens
             + (lgamma(alpha + beta) + lgamma(alpha + deaths[i]) + lgamma(beta + rats[i] - deaths[i]))
             - (lgamma(alpha) + lgamma(beta) + lgamma(alpha + beta + rats[i])))
  }

  # Return the final posterior density, which is still in logged form
  ldens - 5 / 2 * log(alpha + beta) + log(alpha) + log(beta)
}

# Just defines the size of each contour: 0.05, 0.15, 0.25, ..., 0.95
contours <- seq(0.05, 0.95, 0.1)

# Do the same steps as above, but refine the grid space
# Also, I changed the length to 200, because 2001 was slowing my computer to a crawl
u2 <- seq(-2.3, -1.3, length = 200)
v2 <- seq(1, 5, length = 200)
logdens2 <- outer(u2, v2, log.post2)
# dens2 is a 200x200 matrix of probabilities for the u2, v2 values of alpha and beta
# For instance, u2[1] = -2.3, and v2[1] = 1. dens2[1, 1] = 1.978023e-14
# This is equivalent to saying p(alpha = -2.3, beta = 1) = 1.978023e-14
dens2 <- exp(logdens2 - max(logdens2))
contour(u2, v2, dens2, levels = contours, drawlabels = FALSE, xlim=c(-2.2, -1), ylim=c(1, 5))



modelString ="
model {

  for (j in 1:count) {
    y[j] ~ dbin(theta[j], N[j])
    theta[j] ~ dbeta(alpha, beta)
  }

  y72 ~ dbin(theta72, 30)
  theta72 ~ dbeta(alpha, beta)

  lnx <- log(alpha / beta)
  lny <- log(alpha + beta)

  alpha <- u / pow(v, 2)
  beta <- (1 - u) / pow(v, 2)

  u ~ dunif(0.09, 0.22)
  v ~ dunif(0.08, 0.61)
}
"

writeLines(modelString, con="rats.jags")

#ratdata$theta.bar <- ratdata$y / ratdata$N

ratsModel = jags.model(file="rats.jags", 
                             data=list(y=ratdata$y,
                                       N=ratdata$N, 
                                       count=length(ratdata$N)),
                                       #theta=ratdata$theta.bar),
                             n.chains=4)

update(ratsModel, n.iter=2000)


ratsSamples <- coda.samples(ratsModel, n.iter=50000, variable.names=c("alpha", "beta", "theta", "y", "lnx", "lny", "y72", "theta72"), thin=5)

ratsSamples.M <- as.matrix(ratsSamples)
head(ratsSamples.M)

diagMCMC(codaObject = ratsSamples, parName="alpha")
diagMCMC(codaObject = ratsSamples, parName="beta")
diagMCMC(codaObject = ratsSamples, parName="y72")
diagMCMC(codaObject = ratsSamples, parName="theta72")


hist(ratsSamples.M[,"lnx"], breaks=50, freq=F)
hist(ratsSamples.M[,"beta"], breaks=50, freq=F)
points(ratsSamples.M[,"lnx"], ratsSamples.M[,"lny"], col="red", pch=".")

# 9.3 b.) Simulating a 72nd experiment:
# 9.3 c.) Experiment will contain 30 rats

summary(ratsSamples.M[,"theta72"])
summary(ratsSamples.M[,"y72"])
CI <- quantile(ratsSamples.M[,"y72"], probs = c(0.025, 0.975))
CI
hist(ratsSamples.M[,"y72"], breaks=30, freq=F)
abline(v=CI, col="red")
