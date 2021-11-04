#' The Loglikelihood Function for Overdispersed Multivariate Count Outcome in Repeated Measurements
#' @param pars The vector of fixed effect parameters followed by logarithm of \code{s.u} and logarithm of the overdispersion parameter \code{theta}.
#' @param dataset A list of two elements, where the first element is the multivariate outcome and the second element is the matrix of categorical covariates.
#' @param Des A list of two elements, where the first element is the left hand side of the matrix equation of the loglinear model and the second element is the design matrix for the regression covariates.
#' @param model A model for the multivariate count outcome. Either \code{DMM}, \code{CNBM} or \code{UNBM}.
#'
#' @return The estimated regression parameters of regression model, the loglikelihood values, and the Hessian matrix of the estimated parameters.
#' This loglikelihood function is used to estimate parameters of the regression of multivariate count outcome on categorical covariate in the repeated measurement setting
#' in which the correlation between observation at two different time points is accounted ($\sigma^2_u \neq 0$).




logliklong <- function(pars, dataset, Des, model){

  require(ecoreg)

  theta <- exp(pars[length(pars)])
  u <- exp(pars[(length(pars) - 1)])
  b <- c(0,pars[1:(length(pars) - 2)])

  Cs <- dataset[[1]]
  Xs <- dataset[[2]]


  if(is.null(dim(Xs))){Xs <- as.matrix(Xs,ncol = 1)}

  S <- rowSums(Cs)

  Q <- ncol(Cs)
  if (is.null(nrow(Xs)) == 1){
    N <- length(Xs)/2
  } else {  N <- nrow(Xs)/2
  }


  n.cov <- ncol(Xs)
  if (is.null(n.cov) == 1){
    lvl.cov <- length(levels(as.factor(Xs)))
    n.cov <- 1
    Xs.c <- matrix(Xs,ncol=1)
  } else {
    lvl.cov <- apply(Xs,2,function (x) length(levels(as.factor(x))))
    Xs.c <- Xs
  }


  if (is.null(Des) == 1){
    D.m <- Des.Mat(Q,lvl.cov)
    D <- as.matrix(D.m$Des.mat)
    Y.form <- D.m$LHS
  } else {
    Y.form <- as.matrix(Des[[1]])
    D <- as.matrix(Des[[2]])
  }

  Des <- list("Y.form" = Y.form, "D" = D)

  L <- NULL
  G <- matrix(NA,nrow = 2,ncol = Q)

  if(model == "DMM"){
    for (i in 1:N)
    { G[1,] <- FixEf(Xs.c[i,], Q = Q, b = b, lvl.cov = lvl.cov, Des = Des)
      G[2,] <- FixEf(Xs.c[(N+i),], Q = Q, b = b, lvl.cov = lvl.cov, Des = Des)
    FF1 <- function(z1) FF(z1,u,G , theta, Cs.vec = rbind(Cs[i,],Cs[(N+i),]))
    opt <- try(optim(0.1,FF1,method="BFGS",control=list(fnscale=-1,maxit=5000),hessian=TRUE))
    L[i] <- log(integrate.gh(function(z1) exp(FF1(z1)),mu = opt$par,scale=1/sqrt(-opt$hessian),points=5))
    }
  }

  if(model == "UNBM"){
    for (i in 1:N){
      G[1,] <- FixEf(X[i,], Q, lvl.cov, b, Des)
      G[2,] <- FixEf(X[(N+i),], Q, lvl.cov, b, Des)

      p1 <- function(a) sum(dnbinom(Cs[i,],size = (1/theta)*exp(G[1,] + a),prob = rep(theta/(1+theta),Q),log=TRUE))
      p2 <- function(a) sum(dnbinom(Cs[(N+i),],size = (1/theta)*exp(G[2,] + a),prob = rep(theta/(1+theta),Q),log=TRUE))

      v.g <- function(a) p1(a) + p2(a) + dnorm(a,mean=0,sd = u,log=TRUE)
      opt <- try(optim(0,v.g,method="BFGS",control=list(fnscale=-1,maxit=5000),hessian=TRUE))
      L[i] <- log(integrate.gh(function(a) exp(v.g(a)),mu = opt$par,scale=1/sqrt(-opt$hessian),points=5))
    }
  }

  if(model == "CNBM"){
    for (i in 1:N){
      G[1,] <- c(log(S[i]) + FixEf(Xs[i,],Q, lvl.cov,b,Des))
      G[2,] <- c(log(S[N+i]) + FixEf(Xs[(N+i),],Q,lvl.cov,b,Des))

      p1 <- function(a) sum(dnbinom(Cs[i,],size = (1/theta)*exp(c(G[1,] + a)), prob = rep((1/theta)/(1+(1/theta)),Q),log=TRUE))
      p2 <- function(a) sum(dnbinom(Cs[(N+i),],size = (1/theta)*exp(c(G[2,] + a)), prob = rep((1/theta)/(1+(1/theta)),Q),log=TRUE))
      p3 <- function(a) dnbinom(S[i], size = (1/theta)*sum(exp(c(G[1,] + a))),prob = (1/theta)/(1+(1/theta)),log = TRUE)
      p4 <- function(a) dnbinom(S[(N+i)], size = (1/theta)*sum(exp(c(G[2,] + a))),prob = (1/theta)/(1+(1/theta)),log = TRUE)

      v.g2 <- function(a) (p1(a) + p2(a)) + dnorm(a,mean=0,sd=u,log=TRUE)
      v.h2 <- function(a) (p3(a) + p4(a)) + dnorm(a,mean=0,sd=u,log=TRUE)
      optn <- try(optim(0.1,v.g2,method="BFGS",control=list(fnscale=-1,maxit=5000),hessian=TRUE))
      optd <- try(optim(0.1,v.h2,method="BFGS",control=list(fnscale=-1,maxit=5000),hessian=TRUE))

      L[i] <- log(integrate.gh(function(a) exp(v.g2(a)),mu = optn$par,scale=1/sqrt(-optn$hessian),points=5)) -
        log(integrate.gh(function(a) exp(v.h2(a)),mu = optd$par,scale=1/sqrt(-optd$hessian),points=5))
    }
  }

  res <- sum(L)
  print(res)
  return(res)
}
