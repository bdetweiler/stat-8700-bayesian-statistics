#Normal Data with unknown mean and variance
#Non-Informatiive Prior

n=20
y_bar=50
s2=16

sim=1000000

ch=rchisq(sim, df=n-1)
sig2 = (n-1)*s2/ch

mu=rnorm(sim, mean=y_bar, sd = sqrt(s2/n))

hist(mu, breaks=50, freq=F)


z=seq(30,70,length=1000001)
x=(z-y_bar)/sqrt(s2/n)
lines(z, dt(x, df=n-1)/sqrt(s2/n))

y_new=rnorm(sim, mean=mu, sd=sqrt(sig2))
hist(y_new, breaks=50, freq=F)

x2=(z-y_bar)/sqrt(s2*(1+1/n))
lines(z, dt(x2, df=n-1)/sqrt(s2*(1+1/n)))


#Multinomial Data
#Dirichlet Prior
#Election Example

library(gtools)
sim=rdirichlet(100000, c(728, 584, 138))
diff=sim[,1]-sim[,2]
hist(diff, breaks=30)
sum(sim[,1]-sim[,2]>0)/100000


#Bioassay

# Store Data
bioassay <- data.frame(cbind(x = c(-0.86, -0.30, -0.05, 0.73),
                             n = c(5, 5, 5, 5),
                             y = c(0, 1, 3, 5),
                             y.mod = c(0.5, 1, 3, 4.5)))


logit <- function(x) {
  log(x / (1 - x))
}
 
inv_logit <- function(x) {
  (exp(x)) / (1 + exp(x))
}

posterior <- function(a, b) {
	temp <- 1
	x <- bioassay$x
	y <- bioassay$y
	n <- bioassay$n
	for (i in 1: length(x)) {
	  temp <- temp * (inv_logit(a + (b * x[i]))^y[i]) * ((1 - inv_logit(a + (b * x[i])))^(n[i] - y[i]))
  }
	temp
}


posterior_contour <- function(alpha_min, 
                              alpha_max, 
                              grid_size_alpha, 
                              beta_min, 
                              beta_max, 
                              grid_size_beta, 
                              drawlabels = TRUE) {
  alpha <- seq(alpha_min, alpha_max, length = grid_size_alpha) # Generate a list of alpha values
  beta <- seq(beta_min, beta_max, length = grid_size_beta) # Generate a list of beta values
  post.dens <- outer(alpha, beta, "posterior") # Evaluate the posterior density and all possible combinations of alpha and beta values.  
  scaled.dens <- post.dens / max(post.dens) # Rescale the posterior so values are now relative to height of posterior mode. 
  contour(alpha,
          beta,
          scaled.dens,
          levels = c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95),
          xlab = "alpha", 
          ylab = "beta",
          drawlabels = drawlabels) # Draw contour plot. 
}

posterior_contour(-100, 100, 2001, -100, 100, 2001)
posterior_contour(-25, 25, 2001, -25, 50, 2001)
posterior_contour(-5, 10, 2001, -10, 40, 2001)
posterior_contour(-5, 10, 2001, -10, 40, 2001, FALSE)


##Sampling from the joint posterior

posterior_sample=function(alpha_min, alpha_max, grid_size_alpha, beta_min, beta_max, grid_size_beta, nsim){
	alpha<-seq(alpha_min,alpha_max,length=grid_size_alpha) # Generate a list of alpha values
	beta<-seq(beta_min,beta_max,length=grid_size_beta) # Generate a list of beta values
	post.dens<-outer(alpha,beta,"posterior") # Evaluate the posterior density and all possible combinations of alpha and beta values.  
	scaled.dens<-post.dens / max(post.dens) # Rescale the posterior so values are now relative to height of posterior mode. 
	normalized.dens = scaled.dens / sum(scaled.dens) # Makes sure our discrete probability density sums to 1
	marg.alpha<-rowSums(normalized.dens)
	#Sample 'nsim' values of alpha and beta
	sampled.alpha<-NULL #Define a vector to hold our alpha values
	sampled.beta<-NULL #Define a vector to hold out beta values
	alpha_width=(alpha_max-alpha_min)/(grid_size_alpha-1)
	beta_width=(beta_max-beta_min)/(grid_size_beta-1)
	for(i in 1:nsim){
		temp.alpha<-sample(alpha,1,marg.alpha,replace=TRUE) #Generate a value of alpha from the marginal of alpha
		j<-((temp.alpha - alpha_min)/alpha_width) + 1 #Translate the obtained alpha to the appropriate row of our grid
		conditional.beta<-normalized.dens[j,] / sum(normalized.dens[j,]) # Calculate the condition distribution of beta, conditional on the obtained value of alpha
		temp.beta<-sample(beta,1,conditional.beta,replace=TRUE) # Obtain a beta from the conditional distribution of beta
		temp.alpha<-temp.alpha + runif(1,-alpha_width/2,alpha_width/2) # Add some noise to try to undiscretize the posterior
		temp.beta<-temp.beta + runif(1,-beta_width/2,beta_width/2)
		sampled.alpha<-c(sampled.alpha, temp.alpha) # Add the generated alpha to the list of previously generated alpha
		sampled.beta<-c(sampled.beta, temp.beta)}
	data.frame(sampled.alpha, sampled.beta)
}
sims=posterior_sample(-5, 10, 2001, -10, 40, 2001, 10000)
plot(sims, xlim=c(-5,10), ylim=c(-10,40), pch='.')


# Posterior Distribution of LD50 # LD50 is the dose at which the death rate is 50%, only makes sense for beta>0.
LD50 <- (-sims$sampled.alpha[sims$sampled.beta > 0] / sims$sampled.beta[sims$sampled.beta>0])
hist(LD50,br=35)

# Posterior Predictive distribution for a dose level x = 0.25, n = 5
# y~Bin(5,theta)
# logit(theta) = alpha + beta*x
theta <- inv_logit(sims$sampled.alpha + (sims$sampled.beta * 0.25))
y<-rbinom(1000,5,theta)
barplot(table(y))


