# Testing for micore method
n <- 50
p <- 5
q <- 2

x <- rnorm(n)
X <- cbind(1, x)

counts <- matrix(0, n, p+1)
for (i in 1:n) {
  counts[i,] <- rmultinom(1, size=100, prob=rep(1,p+1))
}

n.burn <- 200
n.samp <- 200

library(micore)
mc.fit <- micore(counts, X, n.burn = n.burn, n.samp = n.samp,
                 n.chain=4, verbose=TRUE)

m.pars <- mergeChainsAll(mc.fit)






