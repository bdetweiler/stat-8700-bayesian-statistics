# STAT 8700 Final Question 3
# Brian Detweiler
# Thursday, December 15th

library(rjags)
set.seed(0)

components <- c(51, 3, 17, 13, 5, 4, 17, 1, 5, 3, 8, 22, 1, 1, 13, 8, 15, 3, 1, 13)

##################################################################################
#                   Question 3 (a)
##################################################################################

fileName <- "Final_3.a.1.jags"

modelString ="
model{
  for (i in 1:count) {
    y[i] ~ dgamma(sh, ra)
  }

  y.pred ~ dgamma(sh, ra)

  gmean <- sh / ra

  # parameterized by mode (m) and standard deviation (sd):
  sh <- 1 + m * ra
  ra <- ( m + sqrt( m^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )


  m ~ dunif(0,100)
  sd ~ dunif(0,100)
}
"

writeLines(modelString, con=fileName)

components.model.1 = jags.model(file=fileName, 
                              data=list(y=components,
                                        count=length(components)),
                              n.chains=4)
        
update(components.model.1, n.iter=1200000)

components.samples.1 <- coda.samples(model = components.model.1,
                                        variable.names = c("y.pred", "y", "gmean", "sh", "ra", "m", "sd"), 
                                        n.iter = 200000, 
                                        thin = 50)

summary(components.samples.1)
diagMCMC(components.samples.1)

components.dic.1 <- dic.samples(model = components.model.1, n.iter = 200000, thin = 50)
components.dic.1

components.samples.M.1 <- as.matrix(components.samples.1)
hist(components.samples.M.1[,"y.pred"], breaks=100, freq=FALSE)




fileName <- "Final_3.a.1.jags"

modelString ="
model{
  for (i in 1:count) {
    y[i] ~ dgamma(sh, ra)
  }

  y.pred ~ dgamma(sh, ra)

  gmean <- sh / ra

  # parameterized by mode (m) and standard deviation (sd):
  sh <- 1 + m * ra
  ra <- ( m + sqrt( m^2 + 4 * sd^2 ) ) / ( 2 * sd^2 )


  m ~ dunif(0,100)
  sd ~ dunif(0,100)
}
"

writeLines(modelString, con=fileName)

components.model.1 = jags.model(file=fileName, 
                              data=list(y=components,
                                        count=length(components)),
                              n.chains=4)
        
update(components.model.1, n.iter=1200000)

components.samples.1 <- coda.samples(model = components.model.1,
                                        variable.names = c("y.pred", "y", "gmean", "sh", "ra", "m", "sd"), 
                                        n.iter = 200000, 
                                        thin = 50)

summary(components.samples.1)
diagMCMC(components.samples.1)

components.dic.1 <- dic.samples(model = components.model.1, n.iter = 200000, thin = 50)
components.dic.1

components.samples.M.1 <- as.matrix(components.samples.1)
hist(components.samples.M.1[,"y.pred"], breaks=100, freq=FALSE)





fileName <- "Final_3.a.2.jags"
modelString ="
model{
  for (i in 1:count) {
    y[i] ~ dweib(a, b)
  }

  y.pred ~ dweib(a, b)

  a ~ dgamma(0.01, 0.01)
  b ~ dgamma(0.01, 0.01)
}
"

writeLines(modelString, con=fileName)

components.model.2 = jags.model(file=fileName, 
                              data=list(y=components,
                                        count=length(components)),
                              n.chains=4)
        
update(components.model.2, n.iter=1200000)

components.samples.2 <- coda.samples(model = components.model.2,
                                        variable.names = c("y.pred", "y", "a", "b"), 
                                        n.iter = 200000, 
                                        thin = 50)

summary(components.samples.2)
diagMCMC(components.samples.2)

components.dic.2 <- dic.samples(model = components.model.2, n.iter = 200000, thin = 50)
components.dic.2

components.samples.M.2 <- as.matrix(components.samples.2)
hist(components.samples.M.2[,"y.pred"], breaks=800, freq=FALSE, xlim=c(0, 150))





fileName <- "Final_3.a.3.jags"
modelString ="
model{
  for (i in 1:count) {
    y[i] ~ dlnorm(mu, tau)
  }

  y.pred ~ dlnorm(mu, tau)

  mu ~ dnorm(mu_0, tau_0)
  tau ~ dunif(0, 1)
 
  mu_0 <- 10.2
  tau_0 <- .01

  sigma2 <- 1 / tau
}
"

writeLines(modelString, con=fileName)

components.model.3 = jags.model(file=fileName, 
                              data=list(y=components,
                                        count=length(components)),
                              n.chains=4)
       
update(components.model.3, n.iter=1200000)

components.samples.3 <- coda.samples(model = components.model.3,
                                        variable.names = c("y.pred", "mu", "tau", "sigma2"), 
                                        n.iter = 200000, 
                                        thin = 50)

summary(components.samples.3)
diagMCMC(components.samples.3)

components.dic.3 <- dic.samples(model = components.model.3, n.iter = 200000, thin = 50)
components.dic.3

components.samples.M.3 <- as.matrix(components.samples.3)
hist(components.samples.M.3[,"y.pred"], breaks=800, freq=FALSE, xlim=c(0, 150))




diffdic(components.dic, components.dic.1)
diffdic(components.dic.1, components.dic.2)
diffdic(components.dic.1, components.dic.3)

##################################################################################
#                   Question 3 (b)
##################################################################################

hist(exp(components.samples.M.3[,"mu"]), breaks=150, freq=FALSE)
CI <- quantile(exp(components.samples.M.3[,"mu"]), probs = c(0.025, 0.975))
abline(v=CI, col="red")


##################################################################################
#                   Question 3 (c)
##################################################################################

hist(components.samples.M.3[,"y.pred"], breaks=900, freq=FALSE, xlim = c(0, 150))
CI <- quantile(components.samples.M.3[,"y.pred"], probs = c(0.025, 0.975))
abline(v=CI, col="red")