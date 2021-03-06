

#' Run Microbiome Covariance Regression
#'
#' @param counts Matrix containing OTU counts.  Rows are samples, columns are OTUs.  The last column is considered the reference OTU for the
#'  additive log-ratio transformation and is not included in the covariance matrix.
#' @param X Design matrix contating covariates of interest. Rows are samples, columns are covariates.  Can include intercept column.
#' @param C0 Prior mean for "C" parameter.  If null, then C0 is set automatically.
#' @param Psi0 Prior inverse-Wishart matrix hyperparameter for "Psi" parameter.  If null, then Psi0 is set automatically.
#' @param Gamma0 Prior inverse-Wishart matrix hyperparameter for "Gamma" parameter.  If null, then Gamma0 is set automatically.
#' @param nuPsi Prior inverse-Wishart degrees of freedom for "Psi" parameter.  If null, then nuPsi is set automatically.
#' @param nuGamma Prior inverse-Wishart degrees of freedom for "Gamma" parameter.  If null, then nuGamma is set automatically.
#' @param target.accept.rate Target acceptance rate for adaptive Metropolis sampler.
#' @param n.samp Number of MCMC samples
#' @param n.burn Number of burn-in samples
#' @param n.chain Number of independent MCMC chains to run in parallel
#' @param n.cores Number of cores to use to parallelize the MCMC chains
#' @param save.eta.cov Logical. Save the proposal distribution variance from the Metropolis sampler?  Caution, this will be large.
#' @param verbose Logical. Output progress?
#' @param adapt.control List contatining tuning parameters for adaptive MCMC part.  See details section below.
#'
#' @return An object of class \code{micore} containing one list element for each MCMC chain.  Each
#' list element contains the following attributes:
#' \itemize{
#'   \item \code{eta}: An array containing the MCMC samples for the eta parameter.  First dimension indexes the MCMC samples.
#'   \item \code{Psi}: An array containing the MCMC samples for the Psi parameter.  First dimension indexes the MCMC samples.
#'   \item \code{A}: An array containing the MCMC samples for the A parameter.  First dimension indexes the MCMC samples.
#'   \item \code{B}: An array containing the MCMC samples for the B parameter.  First dimension indexes the MCMC samples.
#'   \item \code{gamma}: An array containing the MCMC samples for the gamma_i parameters.  First dimension indexes the MCMC samples.
#'   \item \code{Gamma}: An array containing the MCMC samples for the Gamma (matrix) parameter.  First dimension indexes the MCMC samples.
#'   \item \code{eta.accepted}: A matrix indicating which eta samples were accepted (1) and which were rejected (0).  Rows are MCMC iterations, columns are subjects.
#'   \item \code{sigma.zero}: A matrix containing values of the Metropolis variance scaling parameter for all subjects over all MCMC samples.
#'   \item \code{acc.probs}: A matrix containing the Metropolis acceptance probabilities for all subjects over all MCMC samples.
#'   \item \code{counts}: The OTU count matrix.
#'   \item \code{X}: The model matrix.
#' }
#'
#' @details This function allows running the Microbiome Covariance Regession
#' method on multiple chains in parallel.  The user can specify values for
#' the hyperparameters \code{C0}, \code{Psi0}, \code{Gamma0}, \code{nuPsi},
#' and \code{nuGamma}.  If they are not specified, then they will be set
#' automatically.
#'
#' A list of arguments passed to the adaptive Metropolis sampler for \code{eta} can
#' be input using the \code{adapt.control} argument.  These can be modified
#' if the \code{eta} parameters do not converge; the default values are
#' not always optimal.  The parameters (and defaults) that can be edited are:
#' \itemize{
#'   \item \code{init}: (default 0.1) The initial stepsize for the adaptive Metropolis parameters.
#'   \item \code{a}: (default 0.5) The adaptation rate.  A higher value means the adaptations will vanish more quickly with the number of MCMC steps.
#'   \item \code{sigma.zero}: (default 1) The initial value of the adaptive Metropolis variance scaling parameter.
#' }
#'
#' @export
#'
#' @importFrom parallel mclapply
#' @importFrom compositions alr
#'
#' @examples
#' n <- 50
#' p <- 5
#' X <- cbind(1, rnorm(n))
#' counts <- matrix(0, n, p+1)
#' for (i in 1:n) {
#'   counts[i,] <- rmultinom(1, size=100, prob=rep(1,p+1))
#' }
#'
#' library(micore)
#' mc.fit <- micore(counts, X, n.samp=100, n.burn=100, n.chain=1)
#'
micore <- function(counts, X, C0=NULL, Psi0=NULL, Gamma0=NULL, nuPsi=NULL, nuGamma=NULL,
                      n.chain=4, n.cores=n.chain, target.accept.rate=0.23, n.samp=4000, n.burn=4000,
                      adapt.control=NULL, save.eta.cov=FALSE, verbose=FALSE) {

  # Checking input formatting
  if (!is.matrix(counts) | any(counts<0) | any(floor(counts)!=counts)) stop("\"counts\" must be a matrix containing non-negative integers")
  if (!is.matrix(X) | !is.numeric(X)) stop("\"X\" must be a numeric matrix")
  if (n.chain<=0 | floor(n.chain)!=n.chain) stop("\"n.chain\" must be a positive integer")
  if (n.cores<=0 | floor(n.cores)!=n.cores) stop("\"n.cores\" must be a positive integer")
  if (target.accept.rate<=0 | target.accept.rate>=1) stop("\"target.accept.rate\" must be between 0 and 1")
  if (n.samp<1 | floor(n.samp)!=n.samp) stop("\"n.samp\" must be a positive integer")
  if (n.burn<=0 | floor(n.burn)!=n.burn) stop("\"n.burn\" must be a non-negative integer")

  # Checking for dimension mismatches
  if (NROW(counts)!=NROW(X)) stop("\"counts\" and \"X\" should have the same number of rows")

  # Set hyperparameters
  n <- NROW(counts)
  p <- NCOL(counts)-1
  q <- NCOL(X)
  lr.counts <- compositions::alr(counts+1)
  hp <- setHyper(lr.counts, X, C0, Psi0, Gamma0, nuPsi, nuGamma, n, p, q)

  count.in <- lapply(1:n.chain, function(x) return(counts))
  mc.fit <- parallel::mclapply(count.in, mc, lr.counts=lr.counts, Xt=X,
               C0=hp$C0, Psi0=hp$Psi0, Gamma0=hp$Gamma0, nuPsi=hp$nuPsi, nuGamma=hp$nuGamma,
               target.accept.rate=target.accept.rate, n.samp=n.samp, n.burn=n.burn,
               adapt.control=adapt.control, save.eta.cov=save.eta.cov, verbose=verbose)

  class(mc.fit) <- "micore"
  return(mc.fit)
}

#' @importFrom Matrix bdiag
mc <- function(counts, lr.counts, Xt, C0=NULL, Psi0=NULL, Gamma0=NULL, nuPsi=NULL, nuGamma=NULL,
               target.accept.rate=0.23, n.samp=4000, n.burn=4000,
               adapt.control=NULL, save.eta.cov=FALSE, verbose=FALSE) {
  X <- Xt # X parameter clashes in mclapply

  n <- NROW(counts)
  p <- NCOL(counts)-1
  q <- NCOL(X)
  tot.samp <- n.samp + n.burn

  #Initialize parameters
  eta <- lr.counts
  Psi <- cov(eta)
  Psi.inv <- qr.solve(Psi)
  gamma <- rnorm(n)
  Xg <- cbind(X, diag(gamma)%*%X)
  C <- matrix(C0, p, 2*q)
  A <- C[,1:q]
  B <- C[,(q+1):(2*q)]
  obs.prec = qr.solve(crossprod(X))
  Gamma <- n*as.matrix(Matrix::bdiag(obs.prec, obs.prec))
  Gammainv <- qr.solve(Gamma)

  # Containers for sampled values
  eta.s <- array(dim=c(n.samp, n, p))
  Psi.s <- array(dim=c(n.samp, p, p))
  gamma.s <- matrix(0, n.samp, n)
  A.s <- B.s <- array(dim=c(n.samp, p, q))
  sigma.zero.s <- matrix(0, tot.samp, n)
  Gamma.s <- array(dim=c(n.samp, 2*q, 2*q))
  if (save.eta.cov) {
    etacov.s <- array(0, dim=c(tot.samp, n, p, p))
  }

  etacov <- array(0, dim=c(n, p, p))
  etamu <- array(0, dim=c(n, p))
  for (s in 1:n) {
    etacov[s,,] <- 0.01*diag(p)
  }

  if (is.null(adapt.control$a)) {
    a <- 0.5
  } else {
    a <- adapt.control$a
  }

  if (is.null(adapt.control$init)) {
    step.init <- step <- 0.1
  }

  if (is.null(adapt.control$sigma.zero)) {
    sigma.zero <- rep(1, n)
  } else {
    sigma.zero <- adapt.control$sigma.zero
  }

  eta.accepted <- matrix(0, tot.samp, n)
  calc.accept.probs <- matrix(0, tot.samp, n)

  # Sample parameters
  for (i in 1:tot.samp) {
    if (verbose) {
      if (i==1) cat("Beginning burn-in:", "\n")
      if (i==n.burn+1) cat("Beginning sampling:", "\n")
      if (i%%100==0) cat("Sample number: ", i, "\n")
    }

    gamma <- sampgamma_i(eta, X, A, B, Psi.inv)
    Xg <- cbind(X, diag(gamma)%*%X)
    C <- sampC(eta, Xg, Psi, C0, Gammainv)
    A <- C[,1:q]
    B <- C[,(q+1):(2*q)]

    Psi <- sampPsi(eta, Xg, C, C0, Gammainv, Psi0, nuPsi)
    Psi.inv <- qr.solve(Psi)

    Gamma <- sampGamma(C, C0, Psi.inv, Gamma0, nuGamma)
    Gamma.inv <- qr.solve(Gamma)

    etaobject <- sampeta(eta, counts, X, Xg, Psi, C, sigma.zero, etacov)
    eta <- etaobject$samp
    eta.accepted[i,] <- etaobject$accept
    calc.accept.probs[i,] <- etaobject$acc.prob

    if (i>1) {
      yym <- eta - etamu
      sigma.zero <- sigma.zero*exp(step*(etaobject$acc.prob - target.accept.rate))
      etamu <- etamu + step*yym
      for (s in 1:n) {
        etacov[s,,] <- etacov[s,,] + step*(tcrossprod(yym[s,]) - etacov[s,,])
      }

      step <- step.init/(i)^a
    }
    if (save.eta.cov) {
      etacov.s[i,,,] <- etacov
    }
    sigma.zero.s[i,] <- sigma.zero

    if (i > n.burn) {
      idx <- i - n.burn
      eta.s[idx,,] <- eta
      Psi.s[idx,,] <- Psi
      A.s[idx,,] <- C[,1:q]
      B.s[idx,,] <- C[,(q+1):(2*q)]
      gamma.s[idx,] <- gamma
      Gamma.s[idx,,] <- Gamma
    }
  }

  #samp.acc.probs <- colSums(eta.accepted[(n.burn+1):tot.samp,])/n.samp

  ret <- list(eta=eta.s, Psi=Psi.s, A=A.s, B=B.s, gamma=gamma.s, eta.accepted=eta.accepted,
              sigma.zero=sigma.zero.s, Gamma=Gamma.s, acc.probs=calc.accept.probs,
              counts=counts, X=X)
  if (save.eta.cov) {
    ret$etacov <- etacov.s
  }

  return(ret)
}



