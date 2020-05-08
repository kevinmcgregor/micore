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

n.burn <- 100
n.samp <- 100

library(micore)
mc.fit <- micore(counts, X, n.burn = n.burn, n.samp = n.samp,
                 n.chain=4, n.cores=4, verbose=TRUE)

m.pars <- mergeChainsAll(mc.fit)

new.dat <- cbind(c(1,1,1),c(0,1,0))
pred <- predict(mc.fit, newdata = new.dat)
pred.p <- predict(mc.fit, newdata=cbind(c(1,1,1),c(0,1,0)), "prop")
pred.all <- predict(mc.fit, type = "prop", post.stat="median")

p.cov <- getPredCov(mc.fit)
p.cov <- getPredCov(mc.fit, newdata = new.dat)
p.cor <- getPredCov(mc.fit, newdata = new.dat, type = "cor")
p.prec <- getPredCov(mc.fit, newdata = new.dat, type = "prec")
p.pcor <- getPredCov(mc.fit, newdata = new.dat, post.stat = "median", type = "pcor")
p.pcor <- getPredCov(mc.fit, newdata = new.dat, quant=c(0.2, 0.8))

trplot(mc.fit, par="gamma", 3)

print(mc.fit)

