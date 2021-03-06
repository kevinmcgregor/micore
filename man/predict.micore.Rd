% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{predict.micore}
\alias{predict.micore}
\title{Get predicted OTU abundances}
\usage{
\method{predict}{micore}(
  object,
  newdata = NULL,
  type = c("alr", "proportion"),
  post.stat = c("mean", "median"),
  quant = c(0.025, 0.975),
  ...
)
}
\arguments{
\item{object}{An object of class \code{micore}}

\item{newdata}{An optional numeric matrix containing covariates for new observations to get predictions for.}

\item{type}{Character specifying what scale to get predicted OTU abundances for. Additive log-ratios or proportions?}

\item{post.stat}{Character specifying whether the predictions be based on the posterior mean or median.}

\item{quant}{Numeric vector specifying the quantiles of the posterior to return for the OTU abundances.}

\item{...}{Further arguments passed to or from other methods}
}
\value{
A list containing:
\itemize{
   \item \code{fit}: A matrix containing the posterior mean (or median) of the OTU abundances on the appropriate scale
   \item \code{quant}: A list containing the requested quantiles from the posterior mean for the OTU abundances.
}
}
\description{
Get predicted OTU abundances
}
\examples{
n <- 50
p <- 5
X <- cbind(1, rnorm(n))
counts <- matrix(0, n, p+1)
for (i in 1:n) {
  counts[i,] <- rmultinom(1, size=100, prob=rep(1,p+1))
}

library(micore)
mc.fit <- micore(counts, X, n.samp=100, n.burn=100, n.chain=1)

new.dat <- cbind(c(1,1,1),c(0,1,0))
pred <- predict(mc.fit, newdata = new.dat)

pred.p <- predict(mc.fit, newdata=new.dat, "prop")


}
