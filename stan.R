source("DBDA2Eprograms/DBDA2E-utilities.R")
library(rstan)


# Generate some random Bernoulli Data from coin with P(Heads)=0.75
flips=100
outcome=rbinom(flips,1,0.75)


# Write the Stan model file (only need to do this part once)

modelString ="
data {  int<lower=0> N; // N >= 0  int<lower=0,upper=1> y[N]; // y[n] in { 0, 1 }}parameters {  real<lower=0,upper=1> theta; // theta in [0, 1]}model {
	  theta ~ beta(10,10); // prior  y ~ bernoulli(theta); // likelihood}
"

writeLines(modelString, con="bern.stan")

bern_data = list(N=flips, y=outcome)
bern_inits=list(list(theta=0.3), list(theta=0.5), list(theta=0.7), list(theta=0.9))


bern_fit=stan(file='bern.stan', data=bern_data, init=bern_inits, iter=12000000, warmup=2000, chains=4, cores=4, open_progress=T)

print(bern_fit)

bern_mcmc=As.mcmc.list(bern_fit)

#diagMCMC(bern_mcmc)

bernSamples_Matrix=as.matrix(bern_fit)
hist(bernSamples_Matrix[,"theta"], breaks=80, freq=F)
x=seq(0,1,length=10001)
lines(x, dbeta(x, 10+sum(outcome), 10+flips-sum(outcome)), col="red")
