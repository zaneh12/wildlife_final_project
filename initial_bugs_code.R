library(statsecol)
library(jagsUI)
library(MCMCvis)
library(ggplot2)
data("voles")

#First we will create a model for Bugs to integrate

#Setting the priors

psi1 ~ U(0,1)
mu_p ~ U(0,1)
sig_p ~ U(0,1)
alpha_0 ~ U(0,1)
alpha_1 ~ U(0,1)

#State Process

#First we will loop through the rows which are the sites

for(i in 1:n.sites){
  z[i,1] ~ bernoulli(psi1)
  
  #Next we will loop through the columns, which are the primary visits
  #We need to use 2:time because we set the initial value column 1 already
  
  for(t in 2:time){
    
    #we are going to loop through again the sites intermittently to get the
    #gamma and epsilon parameters in order to implement them to the next part of
    #The state process/likelihood model
    
    #Now looping through sites
    for(j in 1:n.sites){
      #This is the collonizaiton with connectivity parameter
      gamma_params[i,j,t-1] = 1-rho_0*voles$Connectivity[j]*z[j,t-1]
    }
    #This is back in the t in 2:time for loop and where we use the gamma param
    #Ans we create the two gamma and epsilon vectors
    gamma[i,t-1] = prod(gamma_params[i,1:n.sites,t-1])
    
    logit(epsilon[i,t-1]) = alpha_0 + alpha_1*voles$Length[i]
    
    z[i,t] ~ dbern(psi[i,t])
    psi[i,t] = (1-z[i,t-1])*gamma[t-1] + z[i,t-1]*(1-epsilon[t-1])
  }
}

#Now to define the observation process 

# Likelihood - Observation process
for (i in 1:n.sites) {
  for(t in 1:time){
    q[i,t] ~ N(mu_p, sig_p)
    logit(p[i,t]) = q[i,t]
    y[i,t] ~ binomial(voles[,8:11][i,t], p[i,t]*z[i,t])
  }
}




















