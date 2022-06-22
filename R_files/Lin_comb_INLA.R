
### R code for producing sets of posterior linear combinations using INLA.

## Here we construct linear combinations based on 'mod.sla' obtained from
## 'INLA_Poi.R' code. The same can be done for 'mod.la'. To do so we need
## to add one more option in the 'inla' function.

library(latex2exp)
library(INLA)

#load("DataPoi.RData")      #uncomment for use

# INLA implementation using the Simplified Laplace Strategy (SLA) to reproduce
# Figure 10 in the paper.

lg.prior <- list(prec = list(prior = "loggamma",
                             param = c(0.1, 0.1)))

f.inla <- y ~ 1 + f(grps, hyper = lg.prior)

mod.sla = inla(formula = f.inla, family = "poisson", 
               data = dat, 
               control.compute = list(config = TRUE, 
                                      return.marginals.predictor = TRUE), 
               control.predictor = list(link = 1),
               selection = list('Predictor' = 9:13),   # selection option
               control.fixed = list(mean.intercept = 0, 
                                    prec.intercept = 0.001))


# Wrapper function to generate joint posterior samples from the selected
# linear predictors of the model (the function exploits 'inla.posterior.sample' function). 

mskew = inla.rjmarginal(1e05, jmarginal = mod.sla$selection)

# Create matrix A of indexes for the linear combinations to approximate and compute.

A = rbind(c(rep(1, 2), rep(0, 3)), c(rep(1, 3), rep(0, 2)), 
          c(rep(1, 4), rep(0, 1)), c(rep(1, 5)))


# Compute deterministic posterior approximations for the defined linear combinations in A.
# We exploit the new functions 'inla.tjmarginal' and 'inla.1djmarginal'.

mdskew = inla.tjmarginal(jmarginal = mod.sla$selection, A = A)

mdskew.res = inla.1djmarginal(jmarginal = mdskew)


# Compute the joint posterior approximations of the defined linear combinations
# using sampling.

mskew.t1 = inla.rjmarginal.eval(function() {Predictor[1]+Predictor[2]}, mskew)
mskew.t2 = inla.rjmarginal.eval(function() {
                                Predictor[1]+Predictor[2]+Predictor[3]
                                }, mskew)
mskew.t3 = inla.rjmarginal.eval(function() {
                                Predictor[1]+Predictor[2]+Predictor[3]+Predictor[4]
                                }, mskew)
mskew.t4 = inla.rjmarginal.eval(function() {
                                Predictor[1]+Predictor[2]+Predictor[3]+Predictor[4]+Predictor[5]
                                }, mskew)


# Figure 10 in the paper. Comparison between the two different approximations.

par(mfrow = c(2,2))


dt1 = density(mskew.t1)
plot(dt1, main = "", col = "blue", lty = 2, lwd = 2, xlab = TeX("$\\eta_9+\\eta_{10}$"))
pt1 = seq(1, length(dt1$x), len = 20)
points(dt1$x[pt1], dt1$y[pt1], pch = 2, col='blue')
dl1 = mdskew.res$`Lin:1`
lines(dl1, col='red', lty = 3, lwd = 2)
pl1 = seq(1, length(dl1[,"x"]), len = 20)
points(dl1[,"x"][pl1], dl1[,"y"][pl1], pch = 3, col='red')
legend("topleft",c("Sampling", "Deter."), col = c('blue', 'red'),
       lty =c(2,3), cex = 0.7)

dt2 = density(mskew.t2)
plot(dt2, main = "", col = "blue", lty = 2, lwd = 2, xlab = TeX("$\\eta_9+\\eta_{10}+\\eta_{11}$"))
pt2 = seq(1, length(dt2$x), len = 20)
points(dt2$x[pt2], dt2$y[pt2], pch = 2, col='blue')
dl2 = mdskew.res$`Lin:2`
lines(dl2, col='red', lty = 3, lwd = 2)
pl2 = seq(1, length(dl2[,"x"]), len = 20)
points(dl2[,"x"][pl2], dl2[,"y"][pl2], pch = 3, col='red')
legend("topleft",c("Sampling", "Deter."), col = c('blue', 'red'),
       lty =c(2,3), cex = 0.7)

dt3 = density(mskew.t3)
plot(dt3, main = "", col = "blue", lty = 2, lwd = 2, xlab = TeX("$\\eta_9+\\eta_{10}+\\eta_{11}+\\eta_{12}$"))
pt3 = seq(1, length(dt3$x), len = 20)
points(dt3$x[pt3], dt3$y[pt3], pch = 2, col='blue')
dl3 = mdskew.res$`Lin:3`
lines(dl3, col='red', lty = 3, lwd = 2)
pl3 = seq(1, length(dl3[,"x"]), len = 20)
points(dl3[,"x"][pl3], dl3[,"y"][pl3], pch = 3, col='red')
legend("topleft",c("Sampling", "Deter."), col = c('blue', 'red'),
       lty =c(2,3), cex = 0.7)

dt4 = density(mskew.t4)
plot(dt4, main = "", col = "blue", lty = 2, lwd = 2, xlab = TeX("$\\eta_9+\\eta_{10}+\\eta_{11}+\\eta_{12}+\\eta_{13}$"))
pt4 = seq(1, length(dt4$x), len = 20)
points(dt4$x[pt4], dt4$y[pt4], pch = 2, col='blue')
dl4 = mdskew.res$`Lin:4`
lines(dl4, col='red', lty = 3, lwd = 2)
pl4 = seq(1, length(dl4[,"x"]), len = 20)
points(dl4[,"x"][pl4], dl4[,"y"][pl4], pch = 3, col='red')
legend("topleft",c("Sampling", "Deter."), col = c('blue', 'red'),
       lty =c(2,3), cex = 0.7)



## Speed Performance (between the two posterior approximations)


# Function for sampling from the joint posterior approximations
# of the sets of linear combinations indexed by A

fun_sam = function(ns){
  aa = inla.rjmarginal(ns, jmarginal = mod.sla$selection)
  
  
  
  fun1 = function(...) {Predictor[1]+Predictor[2]}
  fun2 = function(...) {Predictor[1]+Predictor[2]+Predictor[3]}
  fun3 = function(...) {Predictor[1]+Predictor[2]+Predictor[3]+Predictor[4]}
  fun4 = function(...) {Predictor[1]+Predictor[2]+Predictor[3]+Predictor[4]+Predictor[5]}
  
  
  sa1 = inla.rjmarginal.eval(fun1 , aa)
  sa2 = inla.rjmarginal.eval(fun2 , aa)
  sa3 = inla.rjmarginal.eval(fun3 , aa)
  sa4 = inla.rjmarginal.eval(fun4 , aa)
}



library(microbenchmark)

microbenchmark(inla.tjmarginal(jmarginal = mod.sla$selection, A = A),
               fun_sam(1e03), fun_sam(1e04), times = 1e02, unit = 'ms')

