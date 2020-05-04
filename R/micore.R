

#' Title
#'
#' @param counts
#' @param X
#' @param C0
#' @param V0
#' @param nu0
#' @param Psi0
#' @param Gamma0
#' @param nuGamma
#' @param target.accept.rate
#' @param n.samp
#' @param n.burn
#' @param adapt.control
#' @param save.Ycov
#'
#' @return
#' @export
#'
#' @examples
micore <- function(counts, X, C0, nu0, Psi0, Gamma0, nuGamma,
                      target.accept.rate=0.23, n.samp=1000, n.burn=1000,
                      adapt.control=NULL, save.Ycov=FALSE) {

  #Temp since X argument name clashes in mclapply
  #X <- Xt

  n <- NROW(counts)
  p <- NCOL(counts)-1
  q <- NCOL(X)
  tot.samp <- n.samp + n.burn

  #Initialize parameters
  #Y <- matrix(0, n, p)
  Y <- compositions::alr(counts+1)
  #Y <- matrix(0, n, p) # DEBUG... CHANGE BACK TO LR TRANSFORM AFTER
  Psi <- cov(Y)
  Psi.inv <- qr.solve(Psi)
  gamma <- rnorm(n)
  Xg <- cbind(X, diag(gamma)%*%X)
  C <- matrix(C0, p, 2*q)
  A <- C[,1:q]
  B <- C[,(q+1):(2*q)]
  V0inv <- qr.solve(V0)
  Gamma <- n*as.matrix(Matrix::bdiag(obs.prec, obs.prec))

  # Containers for sampled values
  Y.s <- array(dim=c(n.samp, n, p))
  Psi.s <- array(dim=c(n.samp, p, p))
  gamma.s <- matrix(0, n.samp, n)
  A.s <- B.s <- array(dim=c(n.samp, p, q))
  sigma.zero.s <- matrix(0, tot.samp, n)
  V0.s <- array(dim=c(n.samp, 2*q, 2*q))
  if (save.Ycov) {
    Ycov.s <- array(0, dim=c(tot.samp, n, p, p))
  }

  #Y.warmup <- array(dim=c(n.burn, n, p))
  Ycov <- array(0, dim=c(n, p, p))
  Ymu <- array(0, dim=c(n, p))
  for (s in 1:n) {
    #diag(Ycov[s,,]) <- diag(Psi)
    Ycov[s,,] <- 0.01*diag(p) #Psi
  }
  step.init <- step <- 0.1

  if (is.null(adapt.control$a)) {
    a <- 1e-04
  } else {
    a <- adapt.control$a
  }

  if (is.null(adapt.control$sigma.zero)) {
    sigma.zero <- rep(1, n)
  } else {
    sigma.zero <- adapt.control$sigma.zero
  }

  Y.accepted <- matrix(0, tot.samp, n)
  calc.accept.probs <- matrix(0, tot.samp, n)

  # Sample parameters
  for (i in 1:tot.samp) {
    if (verbose) {
      if (i==1) cat("Beginning burn-in:", "\n")
      if (i==n.burn+1) cat("Beginning sampling:", "\n")
      if (i%%100==0) cat("Sample number: ", i, "\n")
    }

    gamma <- sampGamma(Y, X, A, B, Psi.inv)
    Xg <- cbind(X, diag(gamma)%*%X)
    C <- sampC(Y, Xg, Psi, C0, V0inv, r)
    A <- C[,1:q]
    B <- C[,(q+1):(2*q)]

    Psi <- sampPsi(Y, Xg, C, C0, V0inv, Psi0, nu0)
    Psi.inv <- qr.solve(Psi)

    Gamma <- sampGamma(C, C0, Psi.inv, Gamma0, nuGamma)
    Gamma.inv <- qr.solve(Gamma)

    Yobject <- sampY(Y, counts, X, Xg, Psi, C, sigma.zero, Ycov)
    Y <- Yobject$samp
    Y.accepted[i,] <- Yobject$accept
    calc.accept.probs[i,] <- Yobject$acc.prob

    if (i>1) {
      yym <- Y - Ymu
      sigma.zero <- sigma.zero*exp(step*(Yobject$acc.prob - target.accept.rate))
      Ymu <- Ymu + step*yym
      for (s in 1:n) {
        Ycov[s,,] <- Ycov[s,,] + step*(tcrossprod(yym[s,]) - Ycov[s,,])
      }

      step <- step.init/(i)^a
    }
    if (save.Ycov) {
      Ycov.s[i,,,] <- Ycov
    }
    sigma.zero.s[i,] <- sigma.zero

    if (i > n.burn) {
      idx <- i - n.burn
      Y.s[idx,,] <- Y
      Psi.s[idx,,] <- Psi
      A.s[idx,,] <- C[,1:q]
      B.s[idx,,] <- C[,(q+1):(2*q)]
      gamma.s[idx,] <- gamma
      r.s[idx] <- r
      V0.s[idx,,] <- V0
    }
  }

  samp.acc.probs <- colSums(Y.accepted[(n.burn+1):tot.samp,])/n.samp

  ret <- list(Y=Y.s, Psi=Psi.s, A=A.s, B=B.s, gamma=gamma.s, Y.accepted=Y.accepted, samp.acc.probs=samp.acc.probs,
              sigma.zero=sigma.zero.s, r=r.s, V0=V0.s, acc.probs=calc.accept.probs)
  if (save.Ycov) {
    ret$Ycov <- Ycov.s
  }

  return(ret)
}

# Sample Psi from inverse-Wishart dist
sampPsi <- function(Y, Xg, C, C0, V0inv, Psi0, nu0) {
  n <- NROW(Xg)
  p <- NCOL(Y)
  q <- NCOL(C)/2

  V0inv <- V0inv
  d <- Y - tcrossprod(Xg, C)
  df <- nu0 + n + 2*q
  Wpar <- qr.solve(Psi0 + crossprod(d) + tcrossprod((C-C0)%*%V0inv, (C-C0)))

  return(qr.solve(rWishart(1, df, Wpar)[,,1]))
}

# Sample the gamma_i values from normal dist
sampGamma <- function(Y, X, A, B, Psi.inv) {
  n <- NROW(Y)

  bx <- tcrossprod(B, X)
  pi.bx <- Psi.inv%*%bx
  den <- diag(crossprod(bx, pi.bx)) + 1
  ya <- t(Y) - tcrossprod(A, X)
  num <- diag(crossprod(ya, pi.bx))

  return(rnorm(n, num/den, 1/sqrt(den)))
}


sampGamma <- function(C, C0, Psi.inv, Gamma0, nuGamma) {
  p <- NROW(C)
  cmc <- C - C0
  mat <- qr.solve(crossprod(cmc, Psi.inv)%*%cmc + Gamma0)

  return(rWishart(1, nuGamma+p, mat)[,,1])
}



