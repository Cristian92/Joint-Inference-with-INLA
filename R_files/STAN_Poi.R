
### R code for Hierarchical Poisson GLMM model implementation using STAN


require(rstan)

## The .RData file below contains the simulated data used for the paper example.
## You can also create new ones using the file 'Data_sim.R' in the folder.


#load("DataPoi.RData")   #uncomment for use

## Need to rename 'num.grps' variable for STAN implementation.

numgrps = num.grps

## Model formula for STAN. Input for 'stan' function below. The hierarchical
## model structure contains one intercept coefficient 'alpha', a 10-dimensional
## random effect encoded in 're' and one hyperparameter for the random effect
## precision 'tau'.

formula.Stan <- "
  data {
    int < lower = 1 > N;
    int < lower = 1 > numgrps;
    int < lower = 1 > grps[N];
    int < lower = 0 > y[N];
  }
  
  parameters {
    real alpha;
    vector[numgrps] re;
    real < lower = 0 > tau_ran;
  }
  
  transformed parameters {
    vector[N] mu;
    real < lower = 0 > sd_ran;
    sd_ran = 1/sqrt(tau_ran);
    for(i in 1:N) {
      mu[i] = exp(alpha + re[grps[i]]);
    }
  }

  model {
    alpha ~ normal(0, 1000);
    tau_ran ~ gamma(0.1, 0.1);
    re ~ normal(0, sd_ran);
    for(i in 1:N) {
      y[i] ~ poisson(mu[i]);
    }
  }
  "

## STAN Implementation with parallel run using 20 Markov Chains. 



initial.pars <- function() {
  list(alpha = 1, re = rep(0, numgrps), tau_ran = 1.5)
} 

dat.new <- list(N = N, y = y, grps = grps, numgrps = numgrps)

mod.stan <- stan(model_code = formula.Stan, 
                 data = dat.new,
                 pars = c("mu", "alpha", "re", "tau_ran"),
                 iter = 60000, 
                 warmup = 10000,
                 chains = 20,
                 cores = 20,
                 init = initial.pars)


## We save the resulting samples in the object 'output.stan'. 
## NOTE: No need to use thinning here as we use a 'smarter strategy' given by
## HMC. STAN outcomes exactly match with JAGS ones under similar conditions.

output.stan <- extract(mod.stan)

save(output.stan, file = "output_stan.RData")

