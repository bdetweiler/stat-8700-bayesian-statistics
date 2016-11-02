library(rjags)
library(runjags)

sigma<-c(15,10,16,11, 9,11,10,18)
schoolobs<-c(28,8, -3, 7,-1, 1,18,12)
sat.jags=jags.model(file="8schools.jags", data=list('sigma'=sigma, 'schoolobs'=schoolobs, 'N'=length(schoolobs)), n.adapt = 1000)
samps.coda <- coda.samples(sat.jags, c('mu','tau', 'schoolmean'), n.iter=10000, thin=10)
summary(samps.coda)
                           