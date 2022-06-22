# R codes for paper "Joint Posterior Inference for Latent Gaussian Models with R-INLA"

Here we provide the R codes for reproducing the Hierarchical Poisson example of the paper with its plots using INLA and JAGS (the Binomial simulation
relies on similar assumption but using a logit link function). An alternative implementation using STAN is also added on a separate file for MCMC comparison purposes. STAN outcomes exactly match with JAGS ones. 

# Abstract

Efficient Bayesian inference remains a computational challenge in hierarchical models. Simulation-based approaches such as Markov Chain Monte Carlo methods are still popular but have a large computational cost. When dealing with the large class of Latent Gaussian Models, the INLA methodology embedded in the R-INLA software provides accurate Bayesian inference by computing deterministic mixture representation to approximate the joint posterior, from which marginals are computed. The INLA approach has from the beginning been targeting to approximate univariate  posteriors. In this paper we lay out the development foundation of the tools for also providing joint approximations for subsets of the latent field. These approximations inherit Gaussian copula structure and additionally provide corrections for skewness. The same idea is carried forward also to sampling from the mixture representation, which we now can adjust for skewness. 

# Keywords

Bayesian statistics; Computational statistics; Latent Gaussian Models; Markov Chain Monte Carlo

# Authors

- Cristian Chiuchiolo
- Janet van Niekerk
- Haavard Rue

# R_files Folder

All the source R codes to reproduce the Poisson example simulations under the Generalized Linear Mixed Model (GLMM) structure of the paper using INLA, JAGS and STAN softwares are available within this folder. 
