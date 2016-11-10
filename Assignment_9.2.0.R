#########################################
# Brian Detweiler
# Bayesian Statistics
# Assignment_9.2.0.R
# 2016-11-04
#########################################
#library(rstan)
library(rjags)
library(runjags)

source("DBDA2Eprograms/DBDA2E-utilities.R")

fileName <- "Assignment_9_2_0"

x <- c(-0.86, -0.3, -0.05, 0.73)
n <- c(5, 5, 5, 5)
y <- c(0, 1, 3, 5)
y.mod <- c(0.5, 1, 3, 4.5)

data <- cbind(y, n - y)

bivariate_normal <- function(u, v) {
  mu_U <- 0
  mu_V <- 10
  
  sigma_2_U <- 4
  sigma_2_V <- 100
  
  sigma_U <- sqrt(sigma_2_U)
  sigma_V <- sqrt(sigma_2_V)
  
  rho <- 0.5
  
  first <- 1 / (2 * pi * sigma_U * sigma_V * sqrt(1 - rho^2))
  
  second <- exp(-(1 / (2 * (1 - rho^2))) 
                * (  ((u - mu_U)^2 / sigma_2_U) 
                    + ((v - mu_V)^2 / sigma_2_V) 
                    - (2 * rho * (u - mu_U) * (v - mu_V)) / (sqrt(sigma_2_U) * sqrt(sigma_2_V))))
  rval <- first * second
  return(rval)
}

bioassay <- data.frame(cbind(x, n, y, y.mod))

inv_logit <- function(x){
  exp(x) / (1 + exp(x))
}

posterior <- function(a, b) {
  temp <- 1
	x <- bioassay$x
	y <- bioassay$y
	n <- bioassay$n
	for (i in 1: length(x)) {
	  temp <- (temp * (inv_logit(a + (b * x[i]))^y[i]) 
                  * ((1 - inv_logit(a + (b * x[i])))^(n[i] - y[i])))
  }
  
  bivariate_normal(a, b) * temp
}

posterior_contour <- function(alpha_min, 
                              alpha_max, 
                              grid_size_alpha, 
                              beta_min, 
                              beta_max, 
                              grid_size_beta, 
                              drawlabels = TRUE) {
  
  # Generate a list of alpha values
  alpha <- seq(alpha_min, alpha_max, length = grid_size_alpha)
  
  # Generate a list of beta values
  beta <- seq(beta_min, beta_max, length = grid_size_beta)
  
  # Evaluate the posterior density and all possible combinations of alpha and beta values.  
  post.dens <- outer(alpha, beta, 'posterior')
  
  # Rescale the posterior so values are now relative to height of posterior mode. 
  scaled.dens <- post.dens / max(post.dens)
  
  # Draw contour plot. 
  contour(alpha,
          beta,
          scaled.dens, 
          levels = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95),
          xlab="alpha",
          ylab="beta",
          drawlabels = drawlabels)
}

posterior_contour(-2, 4, 2001, -5, 25, 2001)





















N <- length(x)
C <- 10000    # this just has to be large enough to ensure all phi[i]'s < 1
              # Note that JAGS will NOT WARN you if C is too small, but will sample from a truncated distribution

ones <- rep(1, N)
model.data <- list(y=y, x=x, n=n, N=N)

modelString = 
"
model {

  for (i in 1:N) {
    line[i] <- alpha + (beta * x[i])
    p[i] <- ilogit(line[i])
    y[i] ~ dbin(p[i], n[i])
  }

  # 9.2 d.) - Could make the dosage and number of animals part of the data, but we'll keep it simple for now
  px <- alpha + (beta * -0.2)
  y_pred ~ dbin(ilogit(px), 20)

  alpha ~ dunif(-5, 10)    # Prior
  beta ~ dunif(-10, 40)    # Prior

  # Question 9.2 b.)
  # This will contain the percentage of alphas and betas > 0
  beta0 <- step(beta)
  alpha0 <- step(alpha)

  ld50 <- -alpha / beta
}
"

jags.file <- paste0(fileName, ".jags")
writeLines(modelString, con=jags.file)

# There are 3 steps to getting JAGS to generate samples from the posterior distribution. 

# Step 1 is to get all the model information into JAGS, so that it can figure out which algorithm to use to generate the MC. 

bioassayModel = jags.model(file = jags.file,
                           data = model.data,
                           n.chains = 4)

#Step 2 is to run the chains for the burn-in/warm-up period

update(bioassayModel, n.iter=2000)

#Step 3 is to runs and records the MCMC samples that we will use. 

bioassaySamples <- coda.samples(bioassayModel, n.iter=10000, variable.names=c("alpha", "beta", "beta0", "alpha0", "ld50", "y_pred"))

summary(bioassaySamples)

# Plot some diagnostics and save the plots to disk
diagMCMC(codaObject = bioassaySamples, parName=c("alpha"))
diagMCMC(codaObject = bioassaySamples, parName=c("beta"))

# To draw a histogram, we need to combine all the chains into 1 matrix
alpha.M <- as.matrix(bioassaySamples[,"alpha"])
beta.M <- as.matrix(bioassaySamples[,"beta"])
alpha0.M <- as.matrix(bioassaySamples[,"alpha0"])
beta0.M <- as.matrix(bioassaySamples[,"beta0"])
ld50.M <- as.matrix(bioassaySamples[,"ld50"])
y_pred.M <- as.matrix(bioassaySamples[,"y_pred"])

sims <- data.frame(alpha.M, beta.M)
head(sims, 100)
points(sims, xlim = c(-2, 4), ylim = c(-5, 25), pch='.', col="red")


# Here's the LD50 plot 
hist(ld50.M, breaks=100, freq=F, xlim=c(-0.6, 0.4))

# Here's the predictive results for y given 20 animals and an x_i of -0.2
hist(y_pred.M, breaks=20, freq=F)
ci <- quantile(y_pred.M, c(0.025, 0.975))
ci
abline(v=ci, col="red")

