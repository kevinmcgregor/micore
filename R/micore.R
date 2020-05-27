

#' Title
#'
#' @param counts
#' @param X
#' @param C0
#' @param Psi0
#' @param Gamma0
#' @param nuGamma
#' @param target.accept.rate
#' @param n.samp
#' @param n.burn
#' @param nuPsi
#' @param n.chain
#' @param n.cores
#' @param save.eta.cov
#' @param verbose
#' @param adapt.control
#'
#' @return
#' @export
#'
#' @importFrom parallel mclapply
#'
#' @examples
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

# Main function for running micore (not exported)
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

  samp.acc.probs <- colSums(eta.accepted[(n.burn+1):tot.samp,])/n.samp

  ret <- list(eta=eta.s, Psi=Psi.s, A=A.s, B=B.s, gamma=gamma.s, eta.accepted=eta.accepted, samp.acc.probs=samp.acc.probs,
              sigma.zero=sigma.zero.s, Gamma=Gamma.s, acc.probs=calc.accept.probs,
              counts=counts, X=X)
  if (save.eta.cov) {
    ret$etacov <- etacov.s
  }

  return(ret)
}



