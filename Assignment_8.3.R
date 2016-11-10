#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_8.3.R
# 2016-11-04
#########################################

library(rjags)
library(runjags)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_8.3"
years <- c(1976:1985)

fatal.accidents <- c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)
passenger.deaths <- c(734, 516, 754, 877, 814, 362, 764, 809, 223, 1066)
death.rate <- c(0.19, 0.12, 0.15, 0.16, 0.14, 0.06, 0.13, 0.13, 0.03, 0.15)

airline.deaths <- as.data.frame(cbind(years, fatal.accidents, passenger.deaths, death.rate))

data <- fatal.accidents
N <- length(data)

# Write the JAGS model file (only need to do this part once)
modelString = 
  "
model{

for (i in 1:N){
  y[i] ~ dpois(theta) # Likelihood
}

theta ~ dgamma(1, 1) # Prior
}
"

writeLines(modelString, con="pois1.jags")

# There are 3 steps to getting JAGS to generate samples from the posterior distribution. 

# Step 1 is to get all the model information into JAGS, so that it can figure out which algorithm to use to generate the MC. 

poisModel = jags.model(file = "pois1.jags", 
                       data = list(y=data, N=N), 
                       inits = list(list(theta=0), list(theta=10), list(theta=30), list(theta=40)), 
                       n.chains = 4)

#Step 2 is to run the chains for the burn-in/warm-up period

update(poisModel, n.iter=2000)


#Step 3 is to runs and records the MCMC samples that we will use. 

poisSamples <- coda.samples(poisModel, n.iter=10000, variable.names="theta")


# Plot some diagnostics and save the plots to disk

diagMCMC(codaObject = poisSamples, parName="theta")
saveGraph(file = paste0(fileName, "-ThetaPost"), type="pdf" )

# To draw a histogram, we need to combine all the chain 1 matrix
poisSamples_Matrix <- as.matrix(poisSamples)


hist(poisSamples_Matrix[,"theta"], breaks=100, freq=F)
# If you want to see the mean vs. the Sample mean
# abline(v=y.bar, col="red")
mean(poisSamples_Matrix[,"theta"])
sqrt(sd(poisSamples_Matrix[,"theta"]))

CI <- quantile(poisSamples_Matrix[,"theta"], c(0.05, 0.95))
CI[[1]]
CI[[2]]
abline(v=CI[[1]], col="red")
abline(v=CI[[2]], col="red")
# Now we have a credible interval of (24.57245, 27.89124)