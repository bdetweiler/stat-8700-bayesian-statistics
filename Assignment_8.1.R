#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_8.1.R
# 2016-11-04
#########################################

### Some JAGS Examples - Make sure you have step up JAGS based on the instructions on Blackboard
library(rjags)
library(runjags)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_8.1"
# 1. Suppose we have a population described by a Normal Distribution with 
#    known variance $\sigma^2 = 1600$ and unknown mean $\mu$. 
#    4 observations are collected from the population and the corresponding values were: 
#    940, 1040, 910, and 990.
# Bernoulli Data 
data <- c(940, 1040, 910, 990)
N <- length(data)
y.bar <- mean(data)

outcome <- rnorm(N, y.bar, 400)
outcome

# Write the JAGS model file (only need to do this part once)
modelString = 
"
  model{
  
    sigma1_2 = pow(400, 2)
    tau1 = 1 / sigma1_2
  	for (i in 1:N){
  		y[i] ~ dnorm(theta, tau1) # Likelihood
  	}
  
    sigma2_2 = pow(200, 2)
    tau2 = 1 / sigma2_2
  	theta ~ dnorm(1000, tau2) # Prior
  }
"

writeLines(modelString, con="norm.jags")

# There are 3 steps to getting JAGS to generate samples from the posterior distribution. 

# Step 1 is to get all the model information into JAGS, so that it can figure out which algorithm to use to generate the MC. 

normModel <- jags.model(file = "norm.jags", 
                        data = list(y=data, N=N), 
                        inits = list(list(theta=0), list(theta=500), list(theta=1000), list(theta=1600)), 
                        n.chains = 4)

#Step 2 is to run the chains for the burn-in/warm-up period

update(normModel, n.iter=2000)


#Step 3 is to runs and records the MCMC samples that we will use. 

normSamples <- coda.samples(normModel, n.iter=10000, variable.names="theta")


# Plot some diagnostics and save the plots to disk

diagMCMC(codaObject = normSamples, parName="theta")
saveGraph(file = paste0("", "-ThetaPost"), type="pdf" )

# To draw a histogram, we need to combine all the chains into 1 matrix
normSamples_Matrix <- as.matrix(normSamples)


hist(normSamples_Matrix[,"theta"], breaks=100, freq=F)
# If you want to see the mean vs. the Sample mean
# abline(v=y.bar, col="red")

x <- seq(0, 1200, length = 1000)
y <- dnorm(x, mean = y.bar, sd = sqrt(sigma_2 / N))
lines(x, y, type="l", lwd=1, col="red")

summary(normSamples_Matrix)

CI <- quantile(normSamples_Matrix[, "theta"], c(0.05, 0.95))
CI[[1]]
CI[[2]]
abline(v=CI[[1]], col="red")
abline(v=CI[[2]], col="red")
# Now we have a credible interval of (937, 1003)
