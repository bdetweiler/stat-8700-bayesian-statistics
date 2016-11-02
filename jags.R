### Some JAGS Examples - Make sure you have step up JAGS based on the instructions on Blackboard
library(rjags)
library(runjags)

#Also, download the DBDA2Eprograms.zip file from Blackboard, uncompress it and move the produced folder into your working directory
source("DBDA2Eprograms/DBDA2E-utilities.R")

# Bernoulli Data 
flips=100
outcome=rbinom(flips,1,0.75)


# Write the JAGS model file (only need to do this part once)

modelString ="
model{
	for (i in 1:N){
		y[i] ~ dbern(theta) # Likelihood
	}
	theta~dbeta(1,1) # Prior
	}
"

writeLines(modelString, con="bern.jags")

# There are 3 steps to getting JAGS to generate samples from the posterior distribution. 

# Step 1 is to get all the model information into JAGS, so that it can figure out which algorithm to use to generate the MC. 

bernModel = jags.model(file="bern.jags", data=list(y=outcome, N=flips), inits=list(list(theta=0.3), list(theta=0.5), list(theta=0.7), list(theta=0.9)), n.chains=4)

#Step 2 is to run the chains for the burn-in/warm-up period

update(bernModel, n.iter=2000)


#Step 3 is to runs and records the MCMC samples that we will use. 

bernSamples = coda.samples(bernModel, n.iter=10000, variable.names="theta")


# Plot some diagnostics and save the plots to disk

diagMCMC(codaObject = bernSamples, parName="theta")
saveGraph( file=paste0(fileNameRoot,"-ThetaPost") , type="pdf" )

# To draw a histogram, we need to combine all the chains into 1 matrix
bernSamples_Matrix=as.matrix(bernSamples)
hist(bernSamples_Matrix[,"theta"], breaks=50, freq=F)

#We know that this posterior should be a Beta Distribution with parameters equal to the number of successes+1 and the number of failures + 1
x=seq(0,1,length=10001)
lines(x, dbeta(x, 1+sum(y), 1+N-sum(y)), col="red")

summary(bernSamples)

# We could've also modeled this same problem directly using the Binomial Distribution:
Y=sum(y)

# Write the JAGS model file (only need to do this part once)

modelString ="
model{
	Y ~ dbinom(theta, N) # Likelihood
	theta~dbeta(1,1) # Prior
	}
"

writeLines(modelString, con="binom.jags")

binomModel = jags.model(file="binom.jags", data=list(Y=Y, N=N), inits=list(list(theta=0.3), list(theta=0.5), list(theta=0.7), list(theta=0.9)), n.chains=4)

update(binomModel, n.iter=2000)

binomSamples = coda.samples(binomModel, n.iter=10000, variable.names="theta")


diagMCMC(codaObject = binomSamples, parName="theta")
saveGraph( file=paste0(fileNameRoot,"-ThetaPost") , type="pdf" )

binomSamples_Matrix=as.matrix(binomSamples)
hist(binomSamples_Matrix[,"theta"], breaks=50, freq=F)

#We know that this posterior should be a Beta Distribution with parameters equal to the number of successes+1 and the number of failures + 1
x=seq(0,1,length=10001)
lines(x, dbeta(x, 1+Y, 1+N-Y), col="red")
summary(binomSamples)



# Rather than using the jags.model command, you can use the run.jags command instead.  One advantage of this is that you can compute the chains in parallel if your CPU has multiple cores. You should, at most, have one less chain than the number of cores you have available.

# How many cores do you have?
library(parallel)
detectCores()

# Notice there are subtle differences between jags.model and run.jags, file vs model for example.  Also notice that there is no need for a seperate update or coda.samples command. 
binomModel2 = run.jags(model="binom.jags", data=list(Y=Y, N=N), inits=list(list(theta=0.3), list(theta=0.5), list(theta=0.7), list(theta=0.9)), n.chains=4, method="parallel", monitor=c("theta"), burnin=2000, sample=10000)

# To convert to a Coda format. 
binomSamples2=as.mcmc.list(binomModel2)
diagMCMC(binomSamples2)


#Compare two coins
N1=100
Y1=rbinom(1, N1, 0.1)

N2=80
Y2=rbinom(1,N2,0.9)



modelString ="
model{
	for (i in 1:2){
	Y[i] ~ dbinom(theta[i], N[i])# Likelihood
	theta[i]~dbeta(1,1) }# Prior

	diff=theta[1]-theta[2]
	test=step(diff)
	}
"
writeLines(modelString, con="binom2.jags")

binom2Model = jags.model(file="binom2.jags", data=list(Y=c(Y1,Y2), N=c(N1,N2)), n.chains=4)

update(binom2Model, n.iter=2000)

binom2Samples = coda.samples(binom2Model, n.iter=10000, variable.names=c("theta", "diff", "test"))

binom2Samples_Matrix=as.matrix(binom2Samples)

# This will loop through all 4 parameters of the above model (theta[1], theta[2], diff, test) and draw the diagnostic plots for each.
for (i in 1:dim(binom2Samples_Matrix)[2]){ 
	diagMCMC(binom2Samples, parName=c(colnames(binom2Samples_Matrix)[i]))
	
	}
