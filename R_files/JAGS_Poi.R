
### R code for Hierarchical Poisson GLMM model implementation using JAGS

require(R2jags)

## The .RData file below contains the simulated data used for the paper example.
## You can also create new ones using the file 'Data_sim.R' in the folder.


#load("DataPoi.RData")     #uncomment for use


## Model formula for JAGS. Input for 'jags.parallel' function below. The hierarchical
## model structure contains one intercept coefficient 'alpha', a 10-dimensional
## random effect encoded in 're' and one hyperparameter for the random effect
## precision 'tau'.


my.model <- function() {
  for(i in 1:N) {
    y[i] ~ dpois(mu[i])
    log(mu[i]) <- alpha + re[grps[i]]
  }
  for(j in 1:num.grps) {
    re[j] ~ dnorm(0,tau)
  }
  alpha ~ dnorm(0,0.001)
  tau ~ dgamma(0.1,0.1)
  
  sd <- 1/sqrt(tau)
}


## JAGS Implementation with parallel run using 20 Markov Chains. 


mod.jags <- jags.parallel(data = c("y", "grps", "N", "num.grps"),
                          n.chains = 20,
                          n.iter = 6000000,
                          n.burnin = 1000000,
                          n.thin = 100,
                          parameters = c("mu", "alpha", "re", "tau"),
                          inits = list(list(alpha = 1,
                                            re = rep(0, num.grps),
                                            tau = 1.5)),
                          model.file = my.model)$BUGSoutput


## We save the resulting samples in a list named 'output'. 


output <- list()

alpha <- mod.jags$sims.list$alpha
re <- mod.jags$sims.list$re
mu <- mod.jags$sims.list$mu
tau <- mod.jags$sims.list$tau
output$alpha <- alpha
output$re <- re
output$mu <- mu
output$tau <- tau

save(output, file = "output_jags.RData")

