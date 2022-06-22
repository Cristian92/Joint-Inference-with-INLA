

## Poisson Data simulation for the paper example. One can simulate different
## datasets using the code below or simply use the available simulated dataset
## used for the paper, stored in 'DataPoi.RData'.

N <- 50
num.grps <- 10
grps <- sample(1:num.grps, N, rep = T)
re <- rnorm(num.grps, sd = 1.5)
eta <- -1+re[grps]
y <- rpois(N, lambda = exp(eta))
dat <- data.frame(y, grps)