## An examploe of a Gibbs Sampler
require(mvtnorm, quietly=TRUE)
g<-function(a,b,m,s){dmvnorm(x=cbind(a,b),mean=m,sigma=s)}
x=seq(-5,5,length=1001)
y=x
z=outer(x,y,"g",m=c(0,0),s=matrix(c(1,0.8,0.8,1),2,2))
z=z/max(z)
contours<-seq(0.05,0.95,.1)
contour(x,y,z,levels=contours,drawlabels=FALSE,col="red")

iterations<-500
theta<-matrix(NA,2*iterations+1,2)
theta[1,]<-c(-4,-4)
y1<-0
y2<-0
rho<-0.8
points(theta)

for (i in 1:iterations){
mean1<-y1 + rho*(theta[(2*i)-1,2] - y2)
var1<- 1 - (rho^2)
theta[2*i,1]<-rnorm(1,mean1,sqrt(var1))
theta[2*i,2]<-theta[(2*i)-1,2]
lines(theta)
mean2<-y2 + rho*(theta[2*i,1] - y1)
var2<- 1 - (rho^2)
theta[(2*i)+1,2]<-rnorm(1,mean2,sqrt(var2))
theta[(2*i)+1,1]<-theta[2*i,1]
lines(theta)
}
conv<-matrix(NA, iterations/2,2)
for (j in 1:(iterations/2)){
conv[j,]<-theta[iterations+1+(2*j),]}
contour(x,y,z,levels=contours,drawlabels=FALSE)
points(conv,col="blue")

sigma<-matrix(c(1,rho,rho,1),2,2)
mu<-c(0,0)
test<-rmvnorm(250,mu,sigma)
points(test,col="red")

#Metropolis Example
require(mvtnorm, quietly=TRUE)
g<-function(a,b,m,s){dmvnorm(x=cbind(a,b),mean=m,sigma=s)}
x=seq(-3,3,length=1001)
y=x
z=outer(x,y,"g",m=c(0,0),s=matrix(c(1,0,0,1),2,2))
z=z/max(z)
contours<-seq(0.05,0.95,.1)
contour(x,y,z,levels=contours,drawlabels=FALSE)

iterations<-2000
theta1<-matrix(NA,iterations+1,2)
theta2<-matrix(NA,iterations+1,2)
theta3<-matrix(NA,iterations+1,2)
theta4<-matrix(NA,iterations+1,2)

theta1[1,]<-c(-2.5,-2.5)
theta2[1,]<-c(-2.5,2.5)
theta3[1,]<-c(2.5,-2.5)
theta4[1,]<-c(2.5,2.5)

mu<-c(0,0)
Sigma<-matrix(c(1,0,0,1),2,2)
points(theta1,col="red")
points(theta2,col="blue")
points(theta3,col="green")
points(theta4,col="gray")

for (i in 1:iterations){
thetastar1<-NULL
thetastar2<-NULL
thetastar3<-NULL
thetastar4<-NULL
require(mvtnorm)
thetastar1<-rmvnorm(n=1,theta1[i,],0.04*Sigma)
thetastar2<-rmvnorm(n=1,theta2[i,],0.04*Sigma)
thetastar3<-rmvnorm(n=1,theta3[i,],0.04*Sigma)
thetastar4<-rmvnorm(n=1,theta4[i,],0.04*Sigma)

tau1<-((dmvnorm(thetastar1,mu,Sigma))/(dmvnorm(theta1[i,],mu,Sigma)))
tau2<-((dmvnorm(thetastar2,mu,Sigma))/(dmvnorm(theta2[i,],mu,Sigma)))
tau3<-((dmvnorm(thetastar3,mu,Sigma))/(dmvnorm(theta3[i,],mu,Sigma)))
tau4<-((dmvnorm(thetastar4,mu,Sigma))/(dmvnorm(theta4[i,],mu,Sigma)))
temp1<-runif(1,0,1)
temp2<-runif(1,0,1)
temp3<-runif(1,0,1)
temp4<-runif(1,0,1)
if (temp1 < tau1){theta1[i+1,]<-thetastar1}
if (temp1 > tau1){theta1[i+1,]<-theta1[i,]}
if (temp2 < tau2){theta2[i+1,]<-thetastar2}
if (temp2 > tau2){theta2[i+1,]<-theta2[i,]}
if (temp3 < tau3){theta3[i+1,]<-thetastar3}
if (temp3 > tau3){theta3[i+1,]<-theta3[i,]}
if (temp4 < tau4){theta4[i+1,]<-thetastar4}
if (temp4 > tau4){theta4[i+1,]<-theta4[i,]}

lines(theta1,col="red")
lines(theta2,col="blue")
lines(theta3,col="green")
lines(theta4,col="gray")
}

conv1<-matrix(NA, iterations/4,2)
conv2<-matrix(NA, iterations/4,2)
conv3<-matrix(NA, iterations/4,2)
conv4<-matrix(NA, iterations/4,2)

for (j in 1:(iterations/4)){
conv1[j,]<-theta1[(3*iterations/4)+1+j,]
conv2[j,]<-theta2[(3*iterations/4)+1+j,]
conv3[j,]<-theta3[(3*iterations/4)+1+j,]
conv4[j,]<-theta4[(3*iterations/4)+1+j,]}
plot(conv1, xlim=c(-3,3), ylim=c(-3,3),col="red")
points(conv2,col="blue")
points(conv3,col="green")
points(conv4,col="gray")

conv1<-matrix(NA, iterations/2,2)
conv2<-matrix(NA, iterations/2,2)
conv3<-matrix(NA, iterations/2,2)
conv4<-matrix(NA, iterations/2,2)

for (j in 1:(iterations/2)){
conv1[j,]<-theta1[(iterations/2)+1+j,]
conv2[j,]<-theta2[(iterations/2)+1+j,]
conv3[j,]<-theta3[(iterations/2)+1+j,]
conv4[j,]<-theta4[(iterations/2)+1+j,]}

contour(x,y,z,levels=contours,drawlabels=FALSE, xlim=c(-3,3), ylim=c(-3,3))
points(conv1)
points(conv2)
points(conv3)
points(conv4)




points(rmvnorm(2*iterations,mu,Sigma), col="red")
