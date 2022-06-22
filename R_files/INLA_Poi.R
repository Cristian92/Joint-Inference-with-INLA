
### R code for Hierarchical Poisson GLMM model implementation using INLA

## Check 'https://www.r-inla.org/download-install' for download and installation
## of INLA package as it is not on CRAN.

library(latex2exp)
library(INLA)

## Below we provide the R codes for reproducing the plots in the paper.

## Figure 1 in paper.

library(sn)

sn_val <- function(s, mu, sigma)
{
  #delta <- alpha/sqrt(1+alpha^2)
  s[which(is.na(s))] <- 0
  delta <- sign(s)*sqrt((pi/2)*(abs(s)^(2/3)/(((4-pi)/2)^(2/3)+abs(s)^(2/3))))
  alpha <- delta/sqrt(1-delta^2)
  xi <- mu-delta*sqrt((2*sigma^2)/(pi-2*delta^2))
  omega <- sqrt((pi*sigma^2)/(pi-2*delta^2))
  return(list(alpha = alpha, xi = xi, omega = omega))
}

x = seq(-4, 4, len =1000)
xx = seq(-4, 4, len = 20)
skew_val = c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8)
col_seq = c('green', 'yellow', 'blue', 'red', 'purple', 'orange')
for (i in 1:length(skew_val)) {
  sn_list = sn_val(skew_val[i], 0, 1)
  x2 = 0
  x3 = 0
  alpha_x <- sn_list$alpha
  xi_x <- sn_list$xi
  omega_x <- sn_list$omega
  x2 <- qsn(pnorm(x), xi = xi_x,
                      omega = omega_x, 
                      alpha = alpha_x)
  x3 <- qsn(pnorm(xx), xi = xi_x,
                       omega = omega_x, 
                       alpha = alpha_x)
  if (i==1) plot(x, x2, type = "l", col = col_seq[i],
                 xlab = TeX("Standardized values $z_i$"),
                 ylab = '', xlim = c(-4,4), ylim = c(-4,4), lty = i+1, lwd = 2)
  lines(x, x2, col = col_seq[i], xlab = TeX("Standardized values $z_i$"), 
        ylab = '', xlim = c(-4,4), ylim = c(-4,4), lwd = 2, lty = i+1)
  points(xx, x3, pch = i, col = col_seq[i])
}

abline(a = 0, b = 1, lwd = 2)
legend("topleft", legend = c("-0.8", "-0.5", "-0.2", "0.2", "0.5", "0.8"),
       col = col_seq, lty = 2:7, cex = 0.9, title = TeX("Skewness $\\gamma$"), 
       text.font = 4)
mtext(TeX("$\\tilde{F}_i^{-1}(\\Phi(z_i))$"), side = 2, line = 2)
text(-1, 0,TeX("$p_l$"), cex = 2)
text(1, -0.3,TeX("$p_r$"), cex = 2)


## Joint Posterior samples using 'inla.posterior.sample' function in INLA.


# Load Poisson simulated dataset and samples from JAGS implementation.

#load("DataPoi.RData")          #uncomment for use
#load("output_jags.RData")      #uncomment for use

store_mu <- output$mu

store_alpha <- output$alpha

store_re <- output$re

store_tau <- output$tau

lmu <- log(store_mu)

# INLA implementation using the Simplified Laplace Strategy (SLA).

lg.prior <- list(prec = list(prior = "loggamma",
                             param = c(0.1, 0.1)))

f.inla <- y ~ 1 + f(grps, hyper = lg.prior)

mod.sla = inla(formula = f.inla, family = "poisson", 
               data = dat, 
               control.compute = list(config = TRUE, 
                                      return.marginals.predictor = TRUE), 
               control.predictor = list(link = 1),
               #selection = list('Predictor' = 9:13),
               control.fixed = list(mean.intercept = 0, 
                                    prec.intercept = 0.001))

M = 1e05    # number of samples

numpred = 50  # number of linear predictors in the example

store_inl <- inla.posterior.sample(M, mod.sla, skew.corr = FALSE) # mean corr.

store_inlSK <- inla.posterior.sample(M, mod.sla, skew.corr = TRUE)  # skew corr.

sm = sapply(store_inl, function(x) x$latent)

ss = sapply(store_inlSK, function(x) x$latent)


# R code to produce all linear predictor posterior marginals. Only few results
# are reported in the paper to avoid being redundant. The code below only 
# produces right tail plot for posterior marginal linear predictor 9. 
# One may adjust the code to produce any posterior linear predictor marginal
# based on the previous fitted model.


linpred.index <- c(9)    # linear predictor 9

for (i in linpred.index){
  mp <- mod.sla$marginals.linear.predictor[[i]]
  pst = seq(1, length(mp[, 'x']), len = 20)
  jd <- density(lmu[,i], na.rm = T)
  pjd = seq(1, length(jd$x), len = 20)
  id <- density(sm[i, ], na.rm = T)
  pid = seq(1, length(id$x), len = 20)
  idSK <- density(ss[i, ],  na.rm = T)
  pidSK = seq(1, length(idSK$x), len = 20)
  plot(jd, type = "l",  main = "", lty = 2, lwd = 3, 
       ylim = c(0, max(idSK$y)+0.1), xlim = c(0,2))
  points(jd$x[pjd], jd$y[pjd], pch = 2, col = 'black')
  lines(id, col = 'blue', lty = 3, lwd = 3)
  points(id$x[pid], id$y[pid], pch = 3, col = 'blue')
  lines(idSK, col = 'red', lty = 4, lwd = 3)
  points(idSK$x[pidSK], idSK$y[pidSK], pch = 3, col ='red')
  lines(mp, type = 'l', col = 'green', lty = 5, lwd = 3)
  points(mp[,'x'][pst], mp[,'y'][pst], pch = 4, col = 'green')
  legend("topright", c("JAGS", "INLA (Mean Corr.)", "INLA (Skew Corr.)", 
                       "INLA (Marginals)"), 
         col = c('black', 'blue', 'red', 'green'), lty = c(2,3,4,5))
}


# INLA implementation using the Laplace Strategy (LA).


mod.la = inla(formula = f.inla, family = "poisson", 
           data = dat, 
           control.compute = list(config = TRUE, 
                                  return.marginals.predictor = TRUE), 
           control.predictor = list(link = 1),
           control.inla = list(strategy = "laplace"),
           control.fixed = list(mean.intercept = 0, 
                                prec.intercept = 0.001))


#NOTE: To get the comparison with Laplace approximated posterior results you
# only need to replace 'mod.sla' with 'mod.la' above. The code to run is then
# exactly the same. Use again 'inla.posterior.sample' function to generate
# samples from the joint posterior density from the model 'mod.la'.






