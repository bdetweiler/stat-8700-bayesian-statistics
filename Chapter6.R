#Simon Newcomb's measurements of the speed of light, from Stigler(1977).  The data are recorded as deviations from $24,\!800$ nanoseconds.  Table 3.1 of Bayesian Data Analysis.

light.data<-c(28,26,33,24,34,-44,27,16,40,-2,29,22,24,21,25,30,23,29,31,19,24,20,36,32,36,28,25,21,28,29,37,25,28,26,30,32,36,26,30,22,36,23,27,27,28,27,31,27,26,33,26,32,32,24,39,28,24,25,32,25,29,27,28,29,16,23)

#From section 3.2 with a non-informative prior we have that sigma^2|y  is Inv-Chi(n-1,s^2) and mu|sigma^2 is Normal(ybar, sigma^2/n)

light.mean<-mean(light.data)
light.var<-(sd(light.data))^2
n<-length(light.data)

#Simulate 20 values for sigma^2
temp<-rchisq(20,65)
sim.sigma2<-((n-1)*light.var)/temp

#Simulating corresponding values of mu
sim.mu<-rnorm(20,light.mean,sqrt(sim.sigma2/n))

#Simulate replicated data for each mu,sigma^2 combination
replicates<-matrix(rnorm(20*n,sim.mu,sqrt(sim.sigma2)),20,n)

#Draw histograms of replicated data and compare to real data
par(mfrow=c(3,7))
hist(light.data,br=30,col="red")
for (i in (1:20))
{hist(replicates[i,],br=30,xlim=c(-40,60))}

#Draw histogram of minimum data value for each replicate
min.y<-apply(replicates,1,min)
par(mfrow=c(1,1))
hist(min.y,xlim=c(-45,20),br=c(-20,-15,-10,-5,0,5,10,15,20))
lines(rep(min(light.data),2),c(0,10))

#Repeat for a new 200 simulations
temp<-rchisq(200,65)
sim.sigma2<-(65*light.var)/temp
sim.mu<-rnorm(200,light.mean,sqrt(sim.sigma2/n))
replicates<-matrix(rnorm(200*n,sim.mu,sqrt(sim.sigma2)),200,n)

#Draw histogram of minimum data value for each replicate
min.y<-apply(replicates,1,min)
hist(min.y,xlim=c(-45,20),br=c(-20,-15,-10,-5,0,5,10,15,20))
lines(rep(min(light.data),2),c(0,100),col="red")

#Calculate p-value
1-mean(min.y<min(light.data))
1-mean(min.y>min(light.data))

#Draw histogram of sample variance
var.y<-apply(replicates,1,sd)^2
hist(var.y,br=30,xlim=c(50,250))
lines(rep(sd(light.data)^2,2),c(0,100),col="red")

#Calculate p-value
1-mean(var.y<sd(light.data)^2)

#Now using the test value defined on the bottom of page 146
yrep6=apply(replicates,1,sort)[6,]
yrep61=apply(replicates,1,sort)[61,]
y6=sort(light.data)[6]
y61=sort(light.data)[61]
T=abs(y61-sim.mu)-abs(y6-sim.mu)
Trep=abs(yrep61-sim.mu)-abs(yrep6-sim.mu)
plot(T, Trep, xlim=c(-15,15), ylim=c(-15,15))
abline(a=0,b=1)
sum(Trep>T)/200

## Binomial Trials (bottom of page 147)
sim.theta=rbeta(10000,8,14)
yrep=rbinom(10000*20,1,sim.theta)
yrep=matrix(yrep, 10000,20)

switches=function(v){
	temp=0
	for (i in 2:length(v)){
		if (v[i]!=v[i-1]){temp=temp+1}
		
	}
	temp
}

Trep=apply(yrep, 1, switches)
hist(Trep, breaks=seq(-0.5,17.5, by=1))
abline(v=3, col="red")

#p-value
1-sum(Trep<3)/10000
# reversed p-value
1-sum(Trep>3)/10000




#Schools example from Chapter 5, make sure you run the Chapter 5 code first to obtain the 200 simulated values for each of the 8 thetas. These are stored in sim.theta

#Use the posterior predictive distribution to simulate new y values for each school (repeated 200 times)
yrep=rnorm(8*200, sim.theta, schooldata.sigmaj)
yrep=matrix(yrep, 200,8)

T_max=apply(yrep,1,max)

T_min=apply(yrep,1,min)

T_mean=apply(yrep,1,mean)

T_sd=apply(yrep,1,sd)

par(mfrow=c(2,2))
hist(T_max)
abline(v=max(schooldata.y), col="red")
1-sum(T_max<max(schooldata.y))/200
hist(T_min)
abline(v=min(schooldata.y), col="red")
1-sum(T_min<min(schooldata.y))/200
hist(T_mean)
abline(v=mean(schooldata.y), col="red")
1-sum(T_mean<mean(schooldata.y))/200
hist(T_sd)
abline(v=sd(schooldata.y), col="red")
1-sum(T_sd<sd(schooldata.y))/200

