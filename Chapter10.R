#Rejection Sampling examples

triangle=function(x){
	x=as.vector(x)
	if (length(x)==1){
	if (x>-1 & x<1){dens=1-sign(x)*x} else {dens=0}
	dens}
	else{
		for (i in 1:length(x)){x[i]=triangle(x[i])}
		return(x)
		
	}	
}

x=seq(-2,2,length=10001)
plot(x, triangle(x),"l")

# 3 choices for g

#Choice 1: Uniform(-1,1)
#Choice 2: Normal(0,1)
#Choice 3: Normal(0,0.5)

#Find M in each case

#Choice 1:
lines(x, dunif(x,-1,1),col="red") 

plot(x, triangle(x),"l")
lines(x, 2*dunif(x,-1,1),col="red")
#Thus, for choice 1, M=2

#Choice 2:
plot(x, triangle(x),"l")
lines(x, dnorm(x,0,1),col="red") 

plot(x, triangle(x),"l")
lines(x, sqrt(2*pi)*dnorm(x,0,1),col="red") 
#Thus, for choice 2, M=sqrt(2*pi)

#Choice 3
plot(x, triangle(x),"l")
lines(x, dnorm(x,0,0.5),col="red") 

plot(x, triangle(x),"l")
lines(x, 0.5*sqrt(2*pi)*dnorm(x,0,0.5),col="red") 

nsim=10000
n=0
count=0
while (n<nsim){
	sim=runif(1,-1,1)
	count=count+1
	prob=runif(1,0,1)
	if (prob<=(triangle(sim)/(2*dunif(sim,-1,1)))){
		if (n==0){sims=sim} else {sims=c(sims,sim)}
		n=n+1
		}
	}

count
hist(sims, breaks=50, freq=F)
lines(x, triangle(x), col="red")

nsim=10000
n=0
count=0
while (n<nsim){
	sim=rnorm(1,0,1)
	count=count+1
	prob=runif(1,0,1)
	if (prob<=(triangle(sim)/(sqrt(2*pi)*dnorm(sim,0,1)))){
		if (n==0){sims=sim} else {sims=c(sims,sim)}
		n=n+1
		}
	}

count
hist(sims, breaks=50, freq=F)
lines(x, triangle(x), col="red")


nsim=10000
n=0
count=0
while (n<nsim){
	sim=rnorm(1,0,0.5)
	count=count+1
	prob=runif(1,0,1)
	if (prob<=(triangle(sim)/(0.5*sqrt(2*pi)*dnorm(sim,0,1)))){
		if (n==0){sims=sim} else {sims=c(sims,sim)}
		n=n+1
		}
	}

count
hist(sims, breaks=50, freq=F)
lines(x, triangle(x), col="red")

# Importance Sampling
theta=rnorm(10000,0,0.5)
w = triangle(theta)/dnorm(theta,0,0.5)
sum(theta*w)/sum(w)

normalized_w = w/sum(w)

s_eff=1/sum(normalized_w^2)


#Importance Resampling (SIR)
theta=rnorm(100000,0,0.5)
w = triangle(theta)/dnorm(theta,0,0.5)
theta_SIR_wo=sample(theta, 10000,replace=F, prob=w)
theta_SIR_w=sample(theta, 10000,replace=T, prob=w)

hist(theta_SIR_wo, breaks=30, freq=F)
lines(x, triangle(x), col="red")
hist(theta_SIR_w, breaks=30, freq=F)
lines(x, triangle(x), col="red")
