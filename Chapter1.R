football=read.table("football.txt", header=T)

spread=football$spread[1:672]+runif(672, -0.1, 0.1)
outcome=football$favorite[1:672] - football$underdog[1:672] + runif(672, -0.2, 0.2)

plot(spread, outcome, cex=0.5, pch=15, xlab="point spread", ylab="outcome") # Figure 1.1


plot(spread, outcome-spread, cex=0.5, pch=15, xlab="point spread", ylab="outcome-points spread") # Figure 1.2(a)

hist(football$favorite[1:672] - football$underdog[1:672]-football$spread[1:672], breaks=50, freq=F, xlab="outcome - point spread")
x=seq(-100,100, length=1000001)
lines(x, dnorm(x,0,14))  #Figure 1.2(b)