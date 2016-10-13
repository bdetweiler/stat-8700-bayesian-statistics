#DIC calculations for the Schools Example

# Run the Schools code from chapter 5 first

# First case, all schools independent (tau = infinity)

lpd_ind=log(prod(dnorm(schooldata.y, schooldata.y, schooldata.sigmaj)))
dic_ind = -2*lpd_ind+2*8

#Second case, all schools the same (tau = 0)

#Need to find posterior mean of theta first
a=1/(schooldata.sigmaj^2)/(sum(1/(schooldata.sigmaj^2)))
m=weighted.mean(schooldata.y,a)

lpd_same = log(prod(dnorm(schooldata.y, m, schooldata.sigmaj)))
dic_same = -2*lpd_same + 2*1

#Third Case, Hierarcical Model
#First, have to estimate the posterior means of the thetas from simulations
m2 = apply(sim.theta,2,mean)

lpd_hier = log(prod(dnorm(schooldata.y, m2, schooldata.sigmaj)))

#to find the effective number of parameters, we need to calculate the log predictive density at each on of the simulations and then average out.
lpd_sim=function(thetas){
	log(prod(dnorm(schooldata.y, thetas, schooldata.sigmaj)))
}

p_dic=2*(lpd_hier-(sum(apply(sim.theta,1,lpd_sim))/200))

dic_hier = -2*lpd_hier + 2*p_dic


#WAIC calcuations for school example - hierarchical model
#We need to calculate the log likelihood for each combination of school and simulated theta
log_like=matrix(NA,200,8)
for (i in 1:200){
	for (j in 1:8){
	log_like[i,j]=log(dnorm(schooldata.y[j], sim.theta[i,j], schooldata.sigmaj[j]))	
}}

lppd=sum(log(apply(exp(log_like), 2, mean)))
var_post = apply(log_like,2,var)
pwaic2=sum(var_post)
waic=-2*lppd+2*pwaic2

#LOO-CV calculations for the school example - hierarchical model



