## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----install, eval=FALSE-------------------------------------------------
#  if (!require(devtools)) {
#    install.packages("devtools")
#    library(devtools)
#  }
#  install_github("kevinmcgregor/micore", dependencies=TRUE)

## ----sim_and_run---------------------------------------------------------
# Loading in sample American Gut dataset included in package
library(micore)
data(amgut)

counts <- amgut$counts
X <- model.matrix(~BMI, amgut$clin.dat)

# Number of burn-in samples and number of MCMC samples to save
n.burn <- 500
n.samp <- 500

# Running micore
mc.fit <- micore(counts, X, n.burn = n.burn, n.samp = n.samp, 
                 n.chain=4, n.cores=4, verbose=TRUE)

## ----extract_chain-------------------------------------------------------
# Chain 1
tmp <- mc.fit[[1]]
attributes(tmp)
# Chain 2
tmp <- mc.fit[[2]]
attributes(tmp)
# etc...

## ------------------------------------------------------------------------
# Extracting the B parameter from chain 3
B.3 <- mc.fit[[3]]$B
dim(B.3)
# Get 101th sample of B in chain 3
B.3[101,,]

## ----merge---------------------------------------------------------------
# Merging all 4 chains into single array
B.merge <- mergeChains(mc.fit, par="B")
# Mean of B over all chains
apply(B.merge, 2:3, mean)

## ----abund---------------------------------------------------------------
# Want to estimate OTU abundances for individual with x=1.3.
# Create the covariate profile with intercept
x.try <- matrix(c(1, 1.3), nrow=1)
# Get predicted abundances on additive log-ratio scale
p1 <- predict(mc.fit, newdata=x.try)
p1$fit
# Credible intervals
p1$quant

# Get predicted abundances on proportions scale
p2 <- predict(mc.fit, newdata=x.try, type="prop")
p2$fit
# Credible intervals
p2$quant

## ----predCov-------------------------------------------------------------
# Get estimated covariance matrix using covariate profile defined earlier...
cov1 <- getPredCov(mc.fit, newdata=x.try)
# Extract the covariance matrix for the first individual in x.try
cov1$fit[1,,]
# Credible intervals
cov1$quant$`2.5%`[1,,]
cov1$quant$`97.5%`[1,,]

# Partial correlation instead
pc1 <- getPredCov(mc.fit, newdata=x.try, type="pcor")
# Extract the partial correlation matrix for the first individual in x.try
pc1$fit[1,,]
# Credible intervals
pc1$quant$`2.5%`[1,,]
pc1$quant$`97.5%`[1,,]

## ----traceplot-----------------------------------------------------------
trplot(mc.fit, par="B", ind=c(1,1))

## ----traceplot_eta-------------------------------------------------------
trplot(mc.fit, par="eta", ind=c(1,1))

## ----change_a, eval=FALSE------------------------------------------------
#  # Running micore with a different value for "a" parameter in adaptive Metropolis.
#  # Setting a = 0.3.
#  mc.fit <- micore(counts, X, n.burn = n.burn, n.samp = n.samp,
#                   n.chain=4, n.cores=4, verbose=TRUE,
#                   adapt.control = list(a=0.3))

