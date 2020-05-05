

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
micore <- function(counts, X, C0=NULL, Psi0=NULL, Gamma0=NULL, nuPsi=NULL, nuGamma=NULL,
                      target.accept.rate=0.23, n.samp=4000, n.burn=4000,
                      adapt.control=NULL, save.eta.cov=FALSE, verbose=FALSE) {
  n <- NROW(counts)
  p <- NCOL(counts)-1
  q <- NCOL(X)
  tot.samp <- n.samp + n.burn

  # Set hyperparameters
  lr.counts <- compositions::alr(counts+1)
  hp <- setHyper(lr.counts, X, C0, Psi0, Gamma0, nuPsi, nuGamma)

  #Initialize parameters
  eta <- lr.counts
  Psi <- cov(eta)
  Psi.inv <- qr.solve(Psi)
  gamma <- rnorm(n)
  Xg <- cbind(X, diag(gamma)%*%X)
  C <- matrix(hp$C0, p, 2*q)
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
    C <- sampC(eta, Xg, Psi, hp$C0, Gammainv)
    A <- C[,1:q]
    B <- C[,(q+1):(2*q)]

    Psi <- sampPsi(eta, Xg, C, hp$C0, Gammainv, hp$Psi0, hp$nuPsi)
    Psi.inv <- qr.solve(Psi)

    Gamma <- sampGamma(C, hp$C0, Psi.inv, hp$Gamma0, hp$nuGamma)
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
              sigma.zero=sigma.zero.s, Gamma=Gamma.s, acc.probs=calc.accept.probs)
  if (save.eta.cov) {
    ret$etacov <- etacov.s
  }

  return(ret)
}

setHyper <- function(lr.counts, X, C0, Psi0, Gamma0, nuPsi, nuGamma) {
  p <- NCOL(lr.counts)
  q <- NCOL(X)
  if (is.null(C0)) {
    C0 <- runEM(lr.counts, X, n.iter=10)$C
  }
  if (is.null(Psi0)) {
    Psi0 <- cov(lr.counts)
  }
  if (is.null(Gamma0)) {
    Gamma0 <- diag(2*q)
  }
  if (is.null(nuPsi)) {
    nuPsi <- p
  }
  if (is.null(nuGamma)) {
    nuGamma <- 2*q
  }

  return(list(C0=C0, Psi0=Psi0, Gamma0=Gamma0, nuPsi=nuPsi, nuGamma=nuGamma))
}

sampPsi <- function(eta, Xg, C, C0, Gammainv, Psi0, nuPsi) {
  n <- NROW(Xg)
  p <- NCOL(eta)
  q <- NCOL(C)/2

  d <- eta - tcrossprod(Xg, C)
  df <- nuPsi + n + 2*q
  Wpar <- qr.solve(Psi0 + crossprod(d) + tcrossprod((C-C0)%*%Gammainv, (C-C0)))

  return(qr.solve(rWishart(1, df, Wpar)[,,1]))
}

sampgamma_i <- function(eta, X, A, B, Psi.inv) {
  n <- NROW(eta)

  bx <- tcrossprod(B, X)
  pi.bx <- Psi.inv%*%bx
  den <- diag(crossprod(bx, pi.bx)) + 1
  ya <- t(eta) - tcrossprod(A, X)
  num <- diag(crossprod(ya, pi.bx))

  return(rnorm(n, num/den, 1/sqrt(den)))
}


sampGamma <- function(C, C0, Psi.inv, Gamma0, nuGamma) {
  p <- NROW(C)
  cmc <- C - C0
  mat <- qr.solve(crossprod(cmc, Psi.inv)%*%cmc + Gamma0)

  return(rWishart(1, nuGamma+p, mat)[,,1])
}

sampC <- function(eta, Xg, Psi, C0, Gammainv) {
  xgt.inv <- qr.solve(crossprod(Xg) + Gammainv)
  M <- (crossprod(eta, Xg) + C0%*%Gammainv)%*%xgt.inv

  return(rmatnorm(M, Psi, xgt.inv))
}

rmatnorm <- function(M, rowCov, colCov) {
  n.row <- NROW(M)
  n.col <- NCOL(M)

  cc <- chol(colCov)
  rc <- chol(rowCov)

  s.norm <- matrix(rnorm(n.row*n.col), n.row, n.col)

  return(M + crossprod(rc, s.norm)%*%cc)
}

sampeta <- function(etac, counts, X, Xg, Psi, C, sigma.zero, etacov) {
  n <- NROW(etac)
  p <- NCOL(etac)
  q <- NCOL(X)

  # Mean to go into multivar normal dist
  mean.norm <- tcrossprod(Xg, C)

  B <- C[,(q+1):(2*q)]

  # Generate proposal values of eta
  etasamp <- matrix(0, n, p)
  accepted <- rep(0, n)
  acc.prob <- rep(0, n)
  for (i in 1:n) {
    etaprop <- mvnfast::rmvn(1, etac[i,], sigma.zero[i]*etacov[i,,])

    lr.obj <- evalPostLR(counts[i,], etaprop, etac[i,], Psi, mean.norm[i,])
    post.log.ratio <- lr.obj
    r <- log(runif(1))

    # Pick sampled value based on posterior ratio
    if(r<post.log.ratio) {
      etasamp[i,] <- etaprop
      accepted[i] <- 1
    } else {
      etasamp[i,] <- etac[i,]
      accepted[i] <- 0
    }

    acc.prob[i] <- min(1, exp(post.log.ratio))
  }

  return(list(samp=etasamp, accept=accepted, acc.prob=acc.prob))
}

evalPostLR <- function(counts_i, etai_num, etai_den, Psi, mean.norm) {
  p <- length(etai_num)

  pr.num.unscaled <- c(exp(etai_num), 1)
  pr.den.unscaled <- c(exp(etai_den), 1)

  post_log_ratio <- dmultinom(counts_i, prob=pr.num.unscaled, log = TRUE) -
    dmultinom(counts_i, prob=pr.den.unscaled, log = TRUE) +
    mvnfast::dmvn(etai_num, mean.norm, Psi, log = TRUE) -
    mvnfast::dmvn(etai_den, mean.norm, Psi, log = TRUE)

  return(post_log_ratio)
}

runEM <- function(eta, X, n.iter=10, quiet=TRUE, trace=FALSE) {
  n <- NROW(eta)
  p <- NCOL(eta)
  q <- NCOL(X)

    Psi.inv <- qr.solve(cov(eta))
  C <- matrix(1, p, 2*q)

  if (trace) C.tr <- array(0, dim=c(n.iter, p, 2*q))

  for (i in 1:n.iter) {
    A <- C[,1:q]
    B <- C[,(q+1):(2*q)]

    bx <- tcrossprod(B, X)
    pi.bx <- Psi.inv%*%bx
    v <- 1/(diag(crossprod(bx, pi.bx)) + 1)
    ya <- t(eta) - tcrossprod(A, X)
    num <- diag(crossprod(ya, pi.bx))
    m <- num*v
    s <- sqrt(v)

    Xp <- rbind(cbind(X, m*X), cbind(matrix(0, n, q), s*X))
    etap <- rbind(eta, matrix(0, n, p))
    C <- crossprod(etap, Xp)%*%qr.solve(crossprod(Xp))
    if (!quiet) {print(C); cat("\n")}
    if (trace) C.tr[i,,] <- C
  }

  if (trace) {
    ret <- list(C=C, ms=c(m,s), trace=C.tr)
  } else {
    ret <- list(C=C, ms=c(m,s))
  }

  return(ret)
}

getPsiEst <- function(eta, X, C0, ms) {
  n <- NROW(eta)
  p <- NCOL(eta)
  q <- NCOL(X)

  Xp <- rbind(cbind(X, ms[1]*X), cbind(matrix(0, n, q), ms[2]*X))
  etap <- rbind(eta, matrix(0, n, p))
  YXc <- etap - tcrossprod(Xp, C0)
  return(crossprod(YXc)/n)
}

