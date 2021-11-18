library(statsecol)
library(jagsUI)
library(MCMCvis)
library(ggplot2)
data("voles")

sink("Project_3.txt")
cat("
model{
  #First we will create a model for Bugs to integrate
  
  #Setting the priors
  psi1 ~ dunif(0,1)
  mu_p ~ dunif(0,1)
  sig_p ~ dunif(0,1)
  alpha_0 ~ dunif(0,1)
  alpha_1 ~ dunif(0,1)
  beta_0 ~ dunif(0,1)
  beta_1 ~ dunif(0,1)
  
  #State Process
  
  #First we will loop through the rows which are the sites
  
  for(site in 1:sites){
    z[site,1] ~ dbern(psi1)
    logit(gamma[site]) = beta_0 + beta_1 * connectivity[site]
    logit(epsilon[site]) = alpha_0 + alpha_1*length[site]
    
    
    #Next we will loop through the columns, which are the primary visits
    #We need to use 2:time because we set the initial value column 1 already
    
    #Need to define connectivity, area etc
    
    for(time in 2:times){
    
      
      
      z[site,time] ~ dbern(psi[site,time])
      psi[site,time] = (1-z[site,time-1])*gamma[time-1] + z[site,time-1]*(1-epsilon[time-1])
    }
  }
  
  #Now to define the observation process 
  
  # Likelihood - Observation process
  for (site in 1:sites) {
    for(time in 1:times){
      q[site,time] ~ dnorm(mu_p, sig_p)
      logit(p[site,time]) = q[site,time]
      y[site,time] ~ dbinom(p[site,time]*z[site,time], voles_revis[site,time])
    }
  }
}
",fill = TRUE)
sink()

#giving the data its inputs so it is able to use them
voles_data = list(connectivity = voles$Connectivity, 
                  length = voles$Length,
                  sites = nrow(voles),
                  times = 4,
                  voles_revis = voles[,8:11])
# Initial values
voles_inits = function(){
  list(
  psi1 = runif(1, 0, 1) 
  )
}


voles_params <- c("gamma", "epsilon", "p")

# MCMC settings
ni <- 5000
nt <- 4
nb <- 1000
nc <- 3

voles_out <- jags(data = voles_data,
                  inits = voles_inits,
                  parameters.to.save = voles_params,
                  model.file = "Project_3.txt",
                  n.chains = nc,
                  n.iter = ni,
                  n.burnin = nb,
                  n.thin = nt)
MCMCtrace(voles_out,                 #the fitted model
          params = voles_params[1],   #core model parameters
          iter = ni,                 #plot all iterations
          pdf = FALSE,               #DON'T write to a PDF
          type = "trace") 

print(MCMCsummary(voles_out,
            params = voles_params)) #out parameters of interest



