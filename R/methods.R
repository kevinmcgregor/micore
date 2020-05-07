
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

predict.micore <- function(obj, newdata=NULL, type=c("alr", "proportion"),
                           post.stat=c("mean", "median"),
                           quant=c(0.025, 0.975)) {
  if (class(obj)!="micore") stop("obj must be of class \"micore\"")
  type <- match.arg(type)
  post.stat <- match.arg(post.stat)

  if (is.null(newdata)) {
    X <- obj[[1]]$X
  } else {
    X <- newdata
  }

  n <- NROW(X)

  A.s <- mergeChains(obj, "A")

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

getPredCov <- function(obj, newdata=NULL, quant=c(0.025, 0.975),
                   type=c("cov","cor","prec", "pcor"),
                   post.stat=c("mean", "median")) {
  if (class(obj)!="micore") stop("obj must be of class \"micore\"")
  type <- match.arg(type)
  post.stat <- match.arg(post.stat)

  if (is.null(newdata)) {
    X <- obj[[1]]$X
  } else {
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


mergeChains <- function(obj, par=c("A", "B", "Psi", "eta", "gamma", "Gamma")) {
  if (class(obj)!="micore") stop("obj must be of class \"micore\"")
  par <- match.arg(par)
  m <- do.call(abind::abind, args=list(lapply(obj, function(x){x[[par]]}), along=1))
  return(m)
}

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

