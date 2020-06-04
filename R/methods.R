#' @importFrom stats cov
setHyper <- function(lr.counts, X, C0, Psi0, Gamma0, nuPsi, nuGamma, n, p, q) {
  p <- NCOL(lr.counts)
  q <- NCOL(X)
  if (is.null(C0)) {
    C0 <- runEM(lr.counts, X, n.iter=10)$C
  } else {
    if (!is.matrix(C0) | !is.numeric(C0)) stop("\"C0\"  must be a numeric matrix")
    if (NROW(C0)!=p | NCOL(C0)!=2*q) stop("\"C0\"  must have the same number of rows as \"X\" and two times the number of columns as \"X\" ")
  }
  if (is.null(Psi0)) {
    Psi0 <- cov(lr.counts)
  } else {
    if (!is.matrix(Psi0) | !is.numeric(Psi0)) stop("\"Psi0\"  must be a numeric matrix")
    if (NROW(Psi0)!=p | NCOL(Psi0)!=p) stop("\"Psi0\" must be a square matrix with the same number of columns as \"counts\" minus one")
    if (any(eigen(Psi0)$values<=0)) stop("\"Psi0\" must be positive-definite")
  }
  if (is.null(Gamma0)) {
    Gamma0 <- diag(2*q)
  } else {
    if (!is.matrix(Gamma0) | !is.numeric(Gamma0)) stop("\"Gamma0\"  must be a numeric matrix")
    if (NROW(Gamma0)!=2*q | NCOL(Gamma0)!=2*q) stop("\"Gamma0\" must be a square matrix with two times the number of columns as \"X\" ")
    if (any(eigen(Gamma0)$values<=0)) stop("\"Gamma0\" must be positive-definite")
  }
  if (is.null(nuPsi)) {
    nuPsi <- p
  } else {
    if (!is.numeric(nuPsi)) stop("\"nuPsi\" must be numeric")
    if (nuPsi<=p-1) stop("\"nuPsi\" must be greater than the number of columns of \"counts\" minus two")
  }
  if (is.null(nuGamma)) {
    nuGamma <- 2*q
  } else {
    if (!is.numeric(nuGamma)) stop("\"nuGamma\" must be numeric")
    if (nuGamma<=(2*q-1)) stop("\"nuGamma\" must be greater than two times the number of columns of \"X\" minus one")
  }

  return(list(C0=C0, Psi0=Psi0, Gamma0=Gamma0, nuPsi=nuPsi, nuGamma=nuGamma))
}

#' @importFrom stats rWishart
sampPsi <- function(eta, Xg, C, C0, Gammainv, Psi0, nuPsi) {
  n <- NROW(Xg)
  p <- NCOL(eta)
  q <- NCOL(C)/2

  d <- eta - tcrossprod(Xg, C)
  df <- nuPsi + n + 2*q
  Wpar <- qr.solve(Psi0 + crossprod(d) + tcrossprod((C-C0)%*%Gammainv, (C-C0)))

  return(qr.solve(rWishart(1, df, Wpar)[,,1]))
}

#' @importFrom stats rnorm
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


#' @importFrom mvnfast rmvn
#'
#' @importFrom stats runif
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
    r <- log(stats::runif(1))

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

#' @importFrom mvnfast dmvn
#' @importFrom stats dmultinom
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


#' @importFrom stats cov
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


getPredMean <- function(A, X, type=c("alr", "proportion")) {
  type <- match.arg(type)

  Ax <- tcrossprod(X, A)
  if (type=="alr") {
    r <- Ax
  } else {
    exp.Ax <- cbind(exp(Ax), 1)
    r <- exp.Ax/rowSums(exp.Ax)
  }

  return(r)
}

#' Get predicted OTU abundances
#'
#' @param object An object of class \code{micore}
#' @param newdata An optional numeric matrix containing covariates for new observations to get predictions for.
#' @param type Character specifying what scale to get predicted OTU abundances for. Additive log-ratios or proportions?
#' @param post.stat Character specifying whether the predictions be based on the posterior mean or median.
#' @param quant Numeric vector specifying the quantiles of the posterior to return for the OTU abundances.
#' @param ... Further arguments passed to or from other methods
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{fit}: A matrix containing the posterior mean (or median) of the OTU abundances on the appropriate scale
#'    \item \code{quant}: A list containing the requested quantiles from the posterior mean for the OTU abundances.
#' }
#' @export
#' @importFrom stats median quantile
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
#' new.dat <- cbind(c(1,1,1),c(0,1,0))
#' pred <- predict(mc.fit, newdata = new.dat)
#'
#' pred.p <- predict(mc.fit, newdata=new.dat, "prop")
#'
#'
predict.micore <- function(object, newdata=NULL, type=c("alr", "proportion"),
                           post.stat=c("mean", "median"),
                           quant=c(0.025, 0.975), ...) {
  if (class(object)!="micore") stop("object must be of class \"micore\"")
  if (!is.numeric(quant) | any(quant<=0 | quant>=1)) stop("\"quant\" must be a numeric vector with elements between 0 and 1")
  type <- match.arg(type)
  post.stat <- match.arg(post.stat)

  q.orig <- NCOL(object[[1]]$X)

  if (is.null(newdata)) {
    X <- object[[1]]$X
  } else {
    if (!is.numeric(newdata) | !is.matrix(newdata)) stop("\"newdata\" must be a numeric matrix")
    if (NCOL(newdata)!=q.orig) stop("The number of columns of \"newdata\" does not match with the model matrix in \"object\" ")
    X <- newdata
  }

  n <- NROW(X)

  A.s <- mergeChains(object, "A")

  n.s <- dim(A.s)[1]
  p <- NCOL(A.s)
  if (type=="alr") {
    Ax <- array(0, dim=c(n.s, n, p))
  } else {
    Ax <- array(0, dim=c(n.s, n, p+1))
  }

  for (i in 1:n.s) {
    Ax[i,,] <- getPredMean(A.s[i,,], X, type)
  }

  qtile <- vector("list", length = length(quant))
  for (q in 1:length(quant)) {
    qtile[[q]] <- apply(Ax, 2:3, quantile, probs=quant[q])
  }

  fn <- ifelse(post.stat=="mean", mean, median)
  stat <- apply(Ax, 2:3, fn)

  names(qtile) <- paste0(quant*100, "%")

  return(list(fit=stat, quant=qtile))
}

#' @importFrom stats cov2cor
getC <- function(Psi, B, x, type=c("cov","cor","prec", "pcor")) {
  Bx <- B%*%x
  if (type=="prec" | type=="pcor") {
    Psi.inv <- qr.solve(Psi)
    bxp <- tcrossprod(Bx)%*%Psi.inv
    tr.bxp <- sum(diag(bxp))
    ret <- Psi.inv - 1/(1+tr.bxp)*Psi.inv%*%bxp
    if (type=="pcor") {
      ret <- -cov2cor(ret)
      diag(ret) <- 1
    }
  } else {
    ret <- Psi + tcrossprod(Bx)
    if (type=="cor") {
      ret <- cov2cor(ret)
    }
  }

  return(ret)
}

#' Get the predicted covariance matrix based on the \code{micore} model fit. Can also
#' return correlation, precision, and partial correlation matrices.
#'
#' @param obj An object of class \code{micore}.
#' @param newdata n optional numeric matrix containing covariates for new observations to get
#' predicted covariance matrices for.
#' @param quant Numeric vector specifying the quantiles of the posterior to return for predicted covariances.
#' @param type Type of matrix to return: covariance, correlation, precision, or partial correlation.
#' @param post.stat Character specifying whether the predictions be based on the posterior mean or median.
#'
#' @return A list containing:
#' \itemize{
#'    \item \code{fit}: A matrix containing the posterior mean (or median) of the covariance matrix.
#'    \item \code{quant}: A list containing the requested quantiles from the posterior mean for the covariance matrix.
#' }
#' @export
#' @importFrom stats median quantile
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
#' new.dat <- cbind(c(1,1,1),c(0,1,0))
#'
#' c.mat <- getPredCov(mc.fit, new.dat)
#' pc.mat <- getPredCov(mc.fit, new.dat, type="pcor")
#'
getPredCov <- function(obj, newdata=NULL, quant=c(0.025, 0.975),
                   type=c("cov","cor","prec", "pcor"),
                   post.stat=c("mean", "median")) {
  if (class(obj)!="micore") stop("obj must be of class \"micore\"")
  if (!is.numeric(quant) | any(quant<=0 | quant>=1)) stop("\"quant\" must be a numeric vector with elements between 0 and 1")
  type <- match.arg(type)
  post.stat <- match.arg(post.stat)

  q.orig <- NCOL(obj[[1]]$X)

  if (is.null(newdata)) {
    X <- obj[[1]]$X
  } else {
    if (!is.numeric(newdata) | !is.matrix(newdata)) stop("\"newdata\" must be a numeric matrix")
    if (NCOL(newdata)!=q.orig) stop("The number of columns of \"newdata\" does not match with the model matrix in \"obj\" ")
    X <- newdata
  }

  B.s <- mergeChains(obj, "B")
  Psi.s <- mergeChains(obj, "Psi")

  n <- NROW(X)
  n.s <- dim(Psi.s)[1]
  p <- NCOL(Psi.s)
  Sig <- array(0, dim=c(n.s, n, p, p))

  type <- match.arg(type)

  for (i in 1:n.s) {
    for (j in 1:n) {
      x.cur <- X[j,]
      Sig[i,j,,] <- getC(Psi.s[i,,], B.s[i,,], as.vector(x.cur), type)
    }
  }

  qtile <- vector("list", length = length(quant))
  for (q in 1:length(quant)) {
    qtile[[q]] <- apply(Sig, 2:4, quantile, probs=quant[q])
  }

  fn <- ifelse(post.stat=="mean", mean, median)
  stat <- apply(Sig, 2:4, fn)

  names(qtile) <- paste0(quant*100, "%")

  return(list(fit=stat, quant=qtile))
}

getPost <- function(Psi.s, B.s, x, quant=c(0.025, 0.975), type=c("cov","cor","prec", "pcor")) {
  n.s <- dim(Psi.s)[1]
  p <- NCOL(Psi.s)
  Sig <- array(0, dim=c(n.s, p, p))

  type <- match.arg(type)

  for (i in 1:n.s) {
    Sig[i,,] <- getPredCov(Psi.s[i,,], B.s[i,,], as.vector(x), type)
  }

  qtile <- vector("list", length = length(quant))
  for (q in 1:length(quant)) {
    qtile[[q]] <- apply(Sig, 2:3, quantile, probs=quant[q])
  }

  mean <- apply(Sig, 2:3, mean)
  median <- apply(Sig, 2:3, median)

  names(qtile) <- paste0(quant*100, "%")

  return(list(mean=mean, median=median, interval=qtile, samples=Sig))
}


#' Merge chains of MCMC run from \code{micore} object
#'
#' @param obj An object of class \code{micore}
#' @param par Character object specifying which parameter to merge chains for.
#'
#' @return An array with the merged MCMC samples from the specified parameter.  The
#' first dimension of the array is the MCMC iteration.
#' @export
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
#' mc.fit <- micore(counts, X, n.samp=100, n.burn=100, n.chain=2)
#'
#' Amerge <- mergeChains(mc.fit, par="B")
#'
#' @importFrom abind abind
#'
mergeChains <- function(obj, par=c("A", "B", "Psi", "eta", "gamma", "Gamma")) {
  if (class(obj)!="micore") stop("obj must be of class \"micore\"")
  par <- match.arg(par)
  m <- do.call(abind::abind, args=list(lapply(obj, function(x){x[[par]]}), along=1))
  return(m)
}

#' Plots a traceplot for the different MCMC chains of a \code{micore} object
#' for a specific parameter.
#'
#' @param obj An object of class \code{micore}
#' @param par Character object specifying which parameter to make a traceplot for.
#' @param ind For matrix parameters, a 2-dimensional vector specifying the element of
#' the parameter matrix to make the traceplot for.
#' @param xlab Label for the x-axis.
#' @param ylab Label for the y-axis.
#' @param ... Other arguments passed to \code{plot}.
#'
#' @return Makes a traceplot for the specified parameter.
#' @export
#' @importFrom grDevices rainbow
#' @importFrom graphics lines plot
#'
#' @examples
#'
#' n <- 50
#' p <- 5
#' X <- cbind(1, rnorm(n))
#' counts <- matrix(0, n, p+1)
#' for (i in 1:n) {
#'   counts[i,] <- rmultinom(1, size=100, prob=rep(1,p+1))
#' }
#'
#' library(micore)
#' mc.fit <- micore(counts, X, n.samp=100, n.burn=100, n.chain=2)
#'
#' trplot(mc.fit, par="B", ind=c(1,2))
#'
trplot <- function(obj, par=c("A", "B", "Psi", "eta", "gamma", "Gamma"), ind, xlab=NULL,
                   ylab=NULL, ...) {
  if (class(obj)!="micore") stop("obj must be of class \"micore\"")
  par <- match.arg(par)
  n.chain <- length(obj)

  if (par=="gamma") {
    b.only = lapply(obj, function(x) return(x[[par]][,ind[1]]))
  } else {
    b.only = lapply(obj, function(x) return(x[[par]][,ind[1],ind[2]]))
  }
  bind = do.call(c, b.only)
  rge = range(bind)

  if (is.null(xlab)) xlab <- "Sample number"
  if (is.null(ylab)) {
    if (par=="gamma") {
      ylab <- paste0(par, "[", ind[1], "]")
    } else {
      ylab <- paste0(par, "[", ind[1], ", ", ind[2], "]")
    }
  }
  col <- rainbow(n.chain)
  for (ch in 1:n.chain) {
    if (ch==1) {
      if (par=="gamma") {
        plot(obj[[ch]][[par]][,ind[1]], type="l", col=col[ch], ylim=rge,
             xlab=xlab, ylab=ylab, ...)
      } else {
        plot(obj[[ch]][[par]][,ind[1],ind[2]], type="l", col=col[ch], ylim=rge,
             xlab=xlab, ylab=ylab, ...)
      }
    } else {
      if (par=="gamma") {
        lines(obj[[ch]][[par]][,ind[1]], col=col[ch])
      } else {
        lines(obj[[ch]][[par]][,ind[1],ind[2]], col=col[ch])
      }
    }
  }
}

#' Print method for \code{micore} object.
#'
#' @param x An object of class \code{micore}
#' @param ... Further arguments passed to or from other methods
#'
#' @return Prints information about \code{micore} object.
#' @export
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
#' print(mc.fit)
#'
print.micore <- function(x, ...) {
  if (class(x)!="micore") stop("x must be of class \"micore\"")

  n <- NROW(x[[1]]$counts)
  p <- NCOL(x[[1]]$counts)
  q <- NCOL(x[[1]]$X)

  n.chain <- length(x)
  n.samp <- dim(x[[1]]$A)[1]

  cat("MiCoRe: Microbiome Covariance Regression\n")
  cat(" -", n, "observations,", p-1, "OTUs (plus 1 reference OTU),", q-1, "covariates\n")
  cat(" -", n.chain, "MCMC chains: each chain is a list element of this object\n")
  cat(" -", n.samp, "MCMC samples\n")
}


