### Some JAGS Examples - Make sure you have step up JAGS based on the instructions on Blackboard
library(rjags)
library(runjags)
library(rstan)

#Also, download the DBDA2Eprograms.zip file from Blackboard, uncompress it and move the produced folder into your working directory
source("DBDA2Eprograms/DBDA2E-utilities.R")

#Generate some data
y=rnorm(100,4, 7)
N=length(y)


modelString ="
model{
for ( i in 1:N){
y[i] ~ dnorm(mu, tau) # Likelihood
}
mu ~ dnorm(mu_0, tau_0)
tau ~ dgamma(nu_0/2, nu_0*sigma20 / 2)

tau_0 <- kappa_0 * tau 

kappa_0 <- 2
mu_0 <- 5
nu_0 <- 3
sigma20 <- 6

sigma2 <- 1/tau

ynew ~ dnorm(mu, tau)
}
"

writeLines(modelString, con="class1.jags")

class1=jags.model(file="class1.jags", data=list(y=y, N=N), n.chains=4)

update(class1, n.iter=2000)

class1_samples = coda.samples(class1, n.iter=10000, variable.names=c("mu", "sigma2","ynew"))

diagMCMC(codaObject = class1_samples, parName="mu")
diagMCMC(codaObject = class1_samples, parName="sigma2")

modelString ="
data {
  int<lower=0> N; // N >= 0
  real y[N];
}
parameters {
real mu;
real<lower=0> sigma2;
}
transformed parameters {
real<lower=0> sigma;
real<lower=0> sd_mu;
sigma = sqrt(sigma2);
sd_mu = sqrt(sigma2 / 2);

}

model {
y ~ normal(mu, sigma); // likelihood
mu ~ normal(5, sd_mu); // prior for mu
sigma2 ~ scaled_inv_chi_square(3,6); // prior for sigma^2
}
"

writeLines(modelString, con="class1.stan")

stan_data1 = list(N=N, y=y)

stan_fit1 = stan(file='class1.stan', data=stan_data1, iter=12000, warmup=2000, chains=4, cores=4)

print(stan_fit1)

accidents=c(24, 25, 31, 31, 22, 21, 26, 20, 16, 22)
deaths = c(734, 516, 754, 877, 814, 362, 764, 809, 223, 1066)
death_rate=c(0.19, 0.12, 0.15, 0.16, 0.14, 0.06, 0.13, 0.13, 0.03, 0.15)

miles=(deaths/death_rate)*(10^8)

years=length(accidents)
miles_1986 = 8*10^11

modelString ="
model{
for (i in 1:N){
y[i] ~ dpois(lambda[i])
lambda[i] <- x[i]*theta
}
theta ~ dgamma(1,1)
y1986 ~ dpois(lambda1986)
lambda1986 <- miles_1986*theta
}
"
writeLines(modelString, con="class2.jags")

class2=jags.model(file="class2.jags", data=list(N=years, y=accidents, x=miles, miles_1986=miles_1986), n.chains=4)

update(class2, n.iter=2000)

class2_samples = coda.samples(class2, n.iter=10000, variable.names=c("lambda","theta","y1986"))

diagMCMC(codaObject = class2_samples, parName="theta")
diagMCMC(codaObject = class2_samples, parName="y1986")

