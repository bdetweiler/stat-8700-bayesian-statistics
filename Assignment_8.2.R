#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_8.2.R
# 2016-11-04
#########################################

### Some JAGS Examples - Make sure you have step up JAGS based on the instructions on Blackboard
library(rjags)
library(runjags)
library(Bolstad)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_8.2"
# 1. Suppose we have a population described by a Normal Distribution with 
#    known variance $\sigma^2 = 1600$ and unknown mean $\mu$. 
#    4 observations are collected from the population and the corresponding values were: 
#    940, 1040, 910, and 990.
# Bernoulli Data 
var <- 16
data <- c(26.8, 26.3, 28.3, 28.5, 26.3, 31.9, 28.5, 27.2, 20.9, 27.5, 28.0, 18.6, 22.3, 25.0, 31.5)
y.bar <- mean(data)
N <- length(data)
sigma_2 <- var

# Write the JAGS model file (only need to do this part once)
modelString = 
"
  model{
  
    sigma1_2 = pow(4, 2)
    tau1 = 1 / sigma1_2
  	for (i in 1:N){
  		y[i] ~ dnorm(theta, tau1) # Likelihood
  	}
  
    sigma2_2 = pow(5, 2)
    tau2 = 1 / sigma2_2
  	theta ~ dnorm(20, tau2) # Prior
  }
"

writeLines(modelString, con="norm2.jags")

# There are 3 steps to getting JAGS to generate samples from the posterior distribution. 

# Step 1 is to get all the model information into JAGS, so that it can figure out which algorithm to use to generate the MC. 

normModel = jags.model(file = "norm2.jags", 
                       data = list(y=data, N=N), 
                       inits = list(list(theta=0), list(theta=10), list(theta=30), list(theta=40)), 
                       n.chains = 4)

#Step 2 is to run the chains for the burn-in/warm-up period

update(normModel, n.iter=2000)


#Step 3 is to runs and records the MCMC samples that we will use. 

normSamples <- coda.samples(normModel, n.iter=10000, variable.names="theta")


# Plot some diagnostics and save the plots to disk

diagMCMC(codaObject = normSamples, parName="theta")
saveGraph(file = paste0(fileName, "-ThetaPost"), type="pdf" )

# To draw a histogram, we need to combine all the chain 1 matrix
normSamples_Matrix <- as.matrix(normSamples)


hist(normSamples_Matrix[,"theta"], breaks=100, freq=F)
# If you want to see the mean vs. the Sample mean
# abline(v=y.bar, col="red")
mean(normSamples_Matrix[,"theta"])
sqrt(sd(normSamples_Matrix[,"theta"]))

x <- seq(20, 30, length = 1000)
y <- dnorm(x, mean = 26.2404092, sd = 1.0114435)
lines(x, y, type="l", lwd=1, col="red")

CI <- quantile(normSamples_Matrix[,"theta"], c(0.05, 0.95))
CI[[1]]
CI[[2]]
abline(v=CI[[1]], col="red")
abline(v=CI[[2]], col="red")
# Now we have a credible interval of (24.57245, 27.89124)