# Read in the rat data - You may need to specify the full path to the file.
ratdata <- read.table("rats.txt", header=TRUE)

# Log Posterior (u, v space)
log.post2 <- function(u, v) {
  
  
  deaths <- ratdata$y
  rats <- ratdata$N
  alpha <- exp(u + v) / (1 + exp(u))
  beta <- exp(v) / (1 + exp(u))
  ldens <- 0
  
  # Loop over each of the 71 experiments
  for(i in 1:length(rats)) {
    
    # There is no gamma function in R - have to use Log Gamma, so the density is logged
    # This is why it's addative rather than multiplicative
    # deaths[i] is the same as y_i
    # rats[i] is the same as n_i
    ldens <- (ldens
             + (lgamma(alpha + beta) + lgamma(alpha + deaths[i]) + lgamma(beta + rats[i] - deaths[i]))
             - (lgamma(alpha) + lgamma(beta) + lgamma(alpha + beta + rats[i])))
  }

  # Return the final posterior density, which is still in logged form
  ldens - 5 / 2 * log(alpha + beta) + log(alpha) + log(beta)
}

# Just defines the size of each contour: 0.05, 0.15, 0.25, ..., 0.95
contours <- seq(0.05, 0.95, 0.1)

# Defines the grid range of alpha and beta (in u,v-space) over which to simulate
u2 <- seq(-2.5, -1, length = 200)
v2 <- seq(1.5, 3, length = 200)

# Run the posterior over each combination of our values for alpha and beta (again, in u,v-space) 
logdens2 <- outer(u2, v2, log.post2)

# Since the density is logged, we can "unlog" it by "e"ing it (I'm using very sophisticated math terms here)
# Note: I don't quite understand subtracting e^{max(logdens2)}.
dens2 <- exp(logdens2 - max(logdens2))

# Draw the contours in our defined space u2, v2
contour(u2, v2, dens2, levels = contours, drawlabels = FALSE)


# Do the same steps as above, but refine the grid space
# Also, I changed the length to 200, because 2001 was slowing my computer to a crawl
u2 <- seq(-2.3, -1.3, length = 200)
v2 <- seq(1, 5, length = 200)
logdens2 <- outer(u2, v2, log.post2)
# dens2 is a 200x200 matrix of probabilities for the u2, v2 values of alpha and beta
# For instance, u2[1] = -2.3, and v2[1] = 1. dens2[1, 1] = 1.978023e-14
# This is equivalent to saying p(alpha = -2.3, beta = 1) = 1.978023e-14
dens2 <- exp(logdens2 - max(logdens2))
contour(u2, v2, dens2, levels = contours, drawlabels = FALSE)

# 10,000 simulations
nsim <- 10000
# apply the "sum" function to the first dimension of the matrix (rows) - so basically a row sum, which is the 
# density of the u's
# Note: I'm not sure why it's not alphas at this point since we're no longer in log space. Little confused here.
dens.u <- apply(dens2, 1, sum)

# Take 10,000 samples, with replacement, of a number from 1 to 200 with a probability of alpha (u)
uindex <- sample(1:length(u2), nsim, replace = TRUE, prob = dens.u)
# Now use those numbers (1 - 200) to pull out values of alpha from the alpha (u) vector
sim.u <- u2[uindex]

# Create an empty vector 10,000 elements long
sim.v <- rep(NA, nsim)

# Loop 10,000 times
for (i in (1:nsim)) {
  
  # Sample a beta (v) and store it into sim.v at location i. 
	sim.v[i] <- sample(v2, 1, prob = dens2[uindex[i], ])
	
}

# Overlay points on the contour plot from the simulated alphas and betas (u's and v's)
points(sim.u, sim.v, col="red", pch='.')

# NOTE: I don't see this anywhere in my notes - Definitely confused at this point
#       Apparently this finally gives us simulated alphas and betas?
sim.alpha <- exp(sim.u + sim.v) / (1 + exp(sim.u))
sim.beta <- exp(sim.v) / (1 + exp(sim.u))

# Well, at that point, we can simulate thetas
# Create an empty matrix 10,000x71 to hold the simulated thetas
theta.sim <- matrix(NA, nrow = nsim, ncol = 71)

# Loop 71 times (there are 71 experiments, so we want simulated thetas for each experiment)
for (j in 1:71) {
  # Simulate random values from a Beta(alpha + y_j, beta + n_j - y_j) - which is our thetas
	theta.sim[,j] <- rbeta(nsim, sim.alpha + ratdata$y[j], sim.beta + ratdata$N[j] - ratdata$y[j])
}

# Build a dataframe with the following variables:
# obs - observed death to rat proportion (y/N)
# q025 - Apply the quantile function to the 10,000 columns of theta.sim to get the 0.025th quantile
# qm - Apply the quantile function to the 10,000 columns of theta.sim to get the median
# q975 - Apply the quantile function to the 10,000 columns of theta.sim to get the 0.975th quantile
aa <- data.frame(obs = ratdata$y / ratdata$N, 
                 q025 = apply(theta.sim, 2, quantile, 0.025), 
                 qm = apply(theta.sim, 2, quantile, 0.5), 
                 q975 = apply(theta.sim, 2, quantile, 0.975))


# Add jitter to the (y/N) proportions - jitter is just a visualizing technique so that identical
# values don't completely obscure each other.
j_obs = jitter(aa$obs, amount=.005)

# Plot it!
plot(j_obs , aa$qm, xlim = c(0, 0.4), ylim = c(0, 0.4), pch = 15, cex = 0.75)

# Draw a 45 degree diagonal to show how the thetas have been pulled halfway in between the middle and the diagonal
abline(a = 0, b = 1)

# Draw the confidence intervals for each theta (this comes from the dataframe where we computed the quantiles)
for (i in 1:71){
	lines(c(j_obs[i], j_obs[i]), c(aa$q025[i], aa$q975[i]), col ="grey")	
}



#Schools

schooldata.y<-c(28,8,-3,7,-1,1,18,12)
schooldata.sigmaj<-c(15,10,16,11,9,11,10,18)

mu.hat<-function(tau){
sum((1/(schooldata.sigmaj^2 + tau^2))*(schooldata.y)) / sum((1/(schooldata.sigmaj^2 + tau^2)))}

V.mu.inv<-function(tau)
{sum((1/(schooldata.sigmaj^2 + tau^2)))}

post.tau<-function(tau){
V.mu.inv(tau)^(-(1/2))*prod(((schooldata.sigmaj^2 + tau^2)^(-1/2))*exp(-((schooldata.y - mu.hat(tau))^2)/(2*(schooldata.sigmaj^2 + tau^2))))}

tau<-seq(0.01,30,length=1000)
post.tau.dens<-apply(as.array(tau),1,FUN="post.tau")
plot(tau,post.tau.dens,"l")

sim.tau<-sample(tau,200,prob=post.tau.dens,replace=TRUE)
sim.mu.hat<-apply(as.array(sim.tau),1,FUN="mu.hat")
sim.Vinv<-apply(as.array(sim.tau),1,FUN="V.mu.inv")
sim.mu<-rnorm(200,sim.mu.hat,(sim.Vinv)^(-1/2))

v.j<-function(tau){1/((1/schooldata.sigmaj^2)+(1/tau^2))}

theta.hat.j<-function(mu,tau){
((schooldata.y/schooldata.sigmaj^2)+(mu/(tau^2)))/((1/schooldata.sigmaj^2)+(1/tau^2))}

# Figure 5.6
plot(tau, theta.hat.j(mu.hat(tau), tau), cex=0.1, col=c("black","red","blue","purple","brown", "darkgreen", "orange","cyan"))

##### Figure 5.7
plot(tau, sqrt(v.j(tau) + (1/((tau^2/schooldata.sigmaj^2 + 1)^2))*(1/(sapply(tau, V.mu.inv)))), cex=0.1, ylim=c(0,20), col=c("black","red","blue","purple","brown", "darkgreen", "orange","cyan"))


sim.vj<-apply(as.array(sim.tau),1,FUN="v.j")
sim.thetahat.j<-matrix(NA,ncol=8,nrow=length(sim.mu))

for(j in (1:8)){
for (i in (1:length(sim.mu))){
sim.thetahat.j[i,j]<-theta.hat.j(sim.mu[i],sim.tau[i])[j]}}

sim.theta<-matrix(NA,ncol=8,nrow=length(sim.mu))
for(j in (1:8)){
for (i in (1:length(sim.mu))){
sim.theta[i,j]<-rnorm(1,sim.thetahat.j[i,j],(sim.vj[i])^(1/2))}}