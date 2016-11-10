library(rjags)
library(runjags)

#Also, download the DBDA2Eprograms.zip file from Blackboard, uncompress it and move the produced folder into your working directory
source("DBDA2Eprograms/DBDA2E-utilities.R")

# Bernoulli Data 
flips=100
outcome=rbinom(flips,1,0.75)

Y=sum(outcome)

binomModel = jags.model(file="binom.jags", data=list(Y=Y, N=flips), inits=list(list(theta=0.3), list(theta=0.5), list(theta=0.7), list(theta=0.9)), n.chains=4)
update(binomModel, n.iter=2000)
binomSamples = coda.samples(binomModel, n.iter=10000, variable.names="theta")
diagMCMC(codaObject = binomSamples, parName="theta")

binomSamples_Matrix=as.matrix(binomSamples)
hist(binomSamples_Matrix[,"theta"], breaks=50, freq=F)

#Same Thing, Using Zeros Trick

modelString ="
model{
	C=10000
	zero ~ dpois(lambda)
	lambda = -l+C
	l = logfact(N) - logfact(Y) - logfact(N-Y) +Y*log(theta) + (N-Y)*log(1-theta) # Likelihood
	theta~dbeta(1,1) # Prior
	}
"

writeLines(modelString, con="binom_zeros.jags")

binomZerosModel = jags.model(file="binom_zeros.jags", data=list(Y=Y, N=flips, zero=0), inits=list(list(theta=0.3), list(theta=0.5), list(theta=0.7), list(theta=0.9)), n.chains=4)

update(binomZerosModel, n.iter=2000)


binomZerosSamples = coda.samples(binomZerosModel, n.iter=10000, variable.names="theta")

diagMCMC(codaObject = binomZerosSamples, parName="theta")


binomZerosSamples_Matrix=as.matrix(binomZerosSamples)
hist(binomZerosSamples_Matrix[,"theta"], breaks=50, freq=F)

x=seq(0,1,length=10001)
lines(x, dbeta(x, 1+Y, 1+flips-Y), col="red")
