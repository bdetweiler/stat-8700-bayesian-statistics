library(rjags)
library(runjags)
library(rstan)
library(car)
library(mvtnorm)
setwd("~/Documents/Classes/Year 2016-2017/Fall 2016/STAT 8700")
source("DBDA2Eprograms/DBDA2E-utilities.R")

data(Davis)
d=Davis[Davis$sex=="F",] # Filter out the males
d=d[complete.cases(d[,c(2,4)]),]
d[4,2] = 57
d[4,3] = 166

plot(d$repwt,d$weight)

m=lm(d$weight~d$repwt)
summary(m)

m2=lm(d$weight~-1+d$repwt)
summary(m2)


X=cbind(rep(1, length(d$repwt)), d$repwt)

y=d$weight
n=length(d$weight)
k=dim(X)[2] -1k


beta_hat = solve(t(X) %*% X) %*% t(X) %*% y
V_beta = solve(t(X) %*% X)

s2 = as.numeric(t(y - (X %*% beta_hat)) %*% (y - (X %*% beta_hat)) / (n - k - 1))
sqrt(s2)


#Posteror Simulations
nsim=10000


chi=rchisq(nsim, n - k - 1)

sigma2 = (n-k-1)*s2/chi

beta=matrix(NA, nrow=nsim, ncol=k+1)
for(i in 1:nsim){
beta[i,] = rmvnorm(1,mean=beta_hat, sigma=V_beta * sigma2[i])}

#Using JAGS instead
modelString ="
model{
	for (i in 1:N){
		y[i] ~ dnorm(mu[i], tau) # Likelihood
		mu[i] <- beta0 + beta1*repwt[i]
	}
	beta0~dnorm(0,1.0E-4) # Prior
	beta1~dnorm(0,1.0E-4)
	tau ~ dgamma(0.01, 0.01)
	}
"

writeLines(modelString, con="c14_1.jags")

c14_1.Model = jags.model(file="c14_1.jags", data=list(y=y, N=n, repwt=d$repwt), n.chains=4)

update(c14_1.Model, n.iter=30000)

c14_1.samples=coda.samples(c14_1.Model, n.iter=4000000, variable.names=c("beta0", "beta1"), thin=400)

diagMCMC(codaObject = c14_1.samples, parName="beta0")
diagMCMC(codaObject = c14_1.samples, parName="beta1")

modelString ="
model{
	for (i in 1:N){
		y[i] ~ dnorm(mu[i], tau) # Likelihood
		mu[i] <- inprod(X[i,],beta[])
	}
	beta[1]~dnorm(0,1.0E-4) # Prior
	beta[2]~dnorm(0,1.0E-4)
	tau ~ dgamma(0.01, 0.01)
	}
"

writeLines(modelString, con="c14_2.jags")
c14_2.Model = jags.model(file="c14_2.jags", data=list(y=y, N=n, X=X), n.chains=4)

update(c14_2.Model, n.iter=30000)

c14_2.samples=coda.samples(c14_2.Model, n.iter=4000000, variable.names=c("beta"), thin=400)

diagMCMC(codaObject = c14_2.samples, parName="beta[1]")

c14_3.Model = run.jags(model="c14_2.jags", data=list(y=y, N=n, X=X), n.chains=4, method="parallel", monitor=c("beta"), burnin=30000, sample=10000, thin=400)

modelString ="
model{
	for (i in 1:N){
		y[i] ~ dnorm(mu[i], tau) # Likelihood
		mu[i] <- inprod(X[i,],beta[])
	}
	for (i in 1:21){
	beta[i] ~ dnorm(0, 1.0E-4)}
	tau ~ dgamma(0.01, 0.01)
	
ynew ~ dnorm(mu_new, tau)
mu_new <- inprod(xnew, beta[])}
"

writeLines(modelString, con="c14_3.jags")
c14_4.Model = run.jags(model="c14_3.jags", data=list(y=y, N=n, X=X, xnew=xnew), n.chains=4, method="parallel", monitor=c("beta","ynew"), burnin=30000, sample=10000, thin=1)



#Election
d1986=read.table("1986.asc.txt", header=F)
d1988=read.table("1988.asc.txt", header=F)
dem_prop1986 = d1986$V4 / (d1986$V4 + d1986$V5)
rep_prop1986 = d1986$V5 / (d1986$V4 + d1986$V5)
win_prop1986 =apply(cbind(dem_prop1986, rep_prop1986), 1, "max")
dem_prop1988 = d1988$V4 / (d1988$V4 + d1988$V5)
rep_prop1988 = d1988$V5 / (d1988$V4 + d1988$V5)
win_party1988 =1*(dem_prop1988>0.5) - 1*(dem_prop1988<0.5)
inc_party = 1*(dem_prop1986==win_prop1986) - 1*(rep_prop1986==win_prop1986)
inc_run=abs(d1988[,3])
inc_prop = ((1-inc_party)*rep_prop1988 + (1+inc_party)*dem_prop1988)/2

d=data.frame(win_prop1986, inc_party, inc_run, inc_prop)
d=d[d$inc_party!=0,]
d=d[d$inc_run!=9,]
d=d[d$win_prop1986<1,]
d=d[d$inc_prop<1,]


N=length(d$inc_prop)

modelString ="
model{
	for (i in 1:N){
		y[i] ~ dnorm(mu[i], tau) # Likelihood
		mu[i] <- beta0 + beta1*inc_run[i] + beta2*win_prop1986[i]+beta3*inc_party[i]
	}
	beta0~dnorm(0,1.0E-4) # Prior
	beta1~dnorm(0,1.0E-4)
	beta2~dnorm(0,1.0E-4)
	beta3~dnorm(0,1.0E-4)
	tau ~ dgamma(0.01, 0.01)
	}
"

writeLines(modelString, con="c14_4.jags")

c14_4.Model = jags.model(file="c14_4.jags", data=list(y=d$inc_prop, N=N, inc_run=d$inc_run, win_prop1986=d$win_prop1986, inc_party=d$inc_party), n.chains=4)
update(c14_4.Model, n.iter=20000)

c14_4=coda.samples(c14_4.Model, n.iter=200000, variable.names=c("beta0", "beta1","beta2","beta3"), thin=100)

diagMCMC(codaObject = c14_4, parName="beta1")


#One-Way ANOVA 
scores=c(84,58,100,51,28,89,97,50,76,83,45,42,83,64,47,83,81,83,34,61,77,69,94,80,55,79)
candidates=c(rep(1,6), rep(2,7), rep(3,7), rep(4,6))
N=length(scores)

modelString ="
model{
	for (i in 1:N){
		y[i] ~ dnorm(mu[i], tau) # Likelihood
		mu[i] <- m + alpha[candidate[i]]}
	m~dnorm(0,1.0E-4) # Prior
	for (j in 1:3){
	alpha[j]~dnorm(0,1.0E-4)	
	}
	alpha[4] <- -sum(alpha[1:3])
	tau ~ dgamma(0.01, 0.01)
	}
"

writeLines(modelString, con="anova.jags")

anova.Model = jags.model(file="anova.jags", data=list(y=scores, N=N, candidate=candidates), n.chains=4)
update(anova.Model, n.iter=20000)

anova.samples=coda.samples(anova.Model, n.iter=20000, variable.names=c("m","alpha"), thin=1)

diagMCMC(codaObject = c14_4, parName="beta1")


#Two-Way ANOVA

y=array(NA, c(2,3,2))
y[1,1,1]=9
y[2,1,1]=14
y[1,2,1]=25
y[2,2,1]=26
y[1,3,1]=23
y[2,3,1]=24
y[1,1,2]=19
y[2,1,2]=29
y[1,2,2]=22
y[2,2,2]=25
y[1,3,2]=12
y[2,3,2]=13

modelString ="
model{
	for (j in 1:3){
	for (k in 1:2){
	for (i in 1:2){
		y[i,j,k] ~ dnorm(mu[j,k], tau)} # Likelihood
		mu[j,k] <- m + alpha[j]+beta[k]+g[j,k]}}
	m~dnorm(0,1.0E-4) # Prior
	alpha[1]~dnorm(0,1.0E-4)
	alpha[2]~dnorm(0,1.0E-4)	
	alpha[3] <- -alpha[1]-alpha[2]
	beta[1]~dnorm(0,1.0E-4)
	beta[2]<- -beta[1]
	g[1,1]~dnorm(0,1.0E-4)
	g[2,1]~dnorm(0,1.0E-4)
	g[1,2] <- -g[1,1]
	g[2,2] <- -g[2,1]
	g[3,1] <- -g[1,1]-g[2,1]
	g[3,2] <- -g[3,1]
	tau ~ dgamma(0.01, 0.01)
	}
"
writeLines(modelString, con="anova2.jags")

anova2.Model = jags.model(file="anova2.jags", data=list(y=y), n.chains=4)
update(anova2.Model, n.iter=20000)

anova2.samples=coda.samples(anova2.Model, n.iter=20000, variable.names=c("m","alpha","beta","g"), thin=1)

anova2_dic.samples=dic.samples(anova2.Model, n.iter=20000, variable.names=c("m","alpha","beta","g"), thin=1)

modelString ="
model{
	for (j in 1:3){
	for (k in 1:2){
	for (i in 1:2){
		y[i,j,k] ~ dnorm(m, tau)}}} # Likelihood

	m~dnorm(0,1.0E-4) # Prior

	tau ~ dgamma(0.01, 0.01)
	}
"
writeLines(modelString, con="anova3.jags")

anova3.Model = jags.model(file="anova3.jags", data=list(y=y), n.chains=4)
update(anova3.Model, n.iter=20000)

anova3_dic.samples=dic.samples(anova3.Model, n.iter=20000, variable.names=c("m","alpha","beta","g"), thin=1)



#Another Regression Example
d = read.table("pres.txt",header=T)
d2=d[-c(1:50),]
d2=d2[complete.cases(d2),]

modelString ="
model{
	for (i in 1:511){
		
		y[i] ~ dnorm(mu[i], tau) # Likelihood
		mu[i] <- inprod(X[i,], beta[])
	}
	for (j in 1:20){
		beta[j] ~ dnorm(0, 1.0E-4)
		
	}
	tau ~ dgamma(0.01, 0.01)

	}
"

writeLines(modelString, con="elect.jags")
elect.Model=jags.model(file="elect.jags", data=list(y=d2$Dvote, X=d2[,5:24]), n.chains=4)

update(elect.Model, n.iter=20000)

elect.samples=coda.samples(elect.Model, n.iter=20000, variable.names=c("beta","yrep"), thin=1)
diagMCMC(codaObject = elect.samples, parName="beta[1]")