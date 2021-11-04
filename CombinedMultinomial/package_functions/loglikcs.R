#' The loglikelihood function for two independent overdispersed multivariate count outcome
#' @param pars The vector of fixed effect parameters followed by the overdispersion parameter \code{theta}.
#' @param dataset A list of two elements, where the first element is the multivariate outcome and the second element is the matrix of categorical covariates.
#' @param Des A list of two elements, where the first element is the left hand side of the matrix equation of the loglinear model and the second element is the design matrix for the regression covariates.
#' @param model A model for the multivariate count outcome. Either \code{DMM}, \code{CNBM} or \code{UNBM}.
#'
#' @return The estimated regression parameters of regression model, the loglikelihood values, and the Hessian matrix of the estimated parameters.
#' This loglikelihood function is used to estimate parameters of the regression of multivariate count outcome on categorical covariate in the repeated measurement setting
#' in which the correlation between observations at two different time point is ignored  ($\sigma^2_u = 0$).


loglikcs <- function(pars, dataset, Des, model){

  theta <- exp(pars[length(pars)])
  b <- c(0,pars[1:(length(pars) - 1)])

  Cs <- dataset[[1]]
  Xs <- dataset[[2]]

  if(is.null(dim(Xs))) {Xs <- as.matrix(Xs,ncol=1)}


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

    for (i in 1:N){
      G[1,] <- FixEf(Xs[i,], Q, lvl.cov, b, Des)
      G[2,] <- FixEf(Xs[(N+i),], Q, lvl.cov, b, Des)
      eta <- exp(G)
      Eta <- (1/theta)*eta
      L[i] <- CNB(Cs[i,],Eta[1,]) + CNB(Cs[(N+i),], Eta[2,])
    }
  }

  if(model == "UNBM"){
    for (i in 1:N){
      G[1,] <- FixEf(X[i,], Q, lvl.cov, b, Des)
      G[2,] <- FixEf(X[(N+i),], Q, lvl.cov, b, Des)

      p1 <- sum(dnbinom(Cs[i,],size = (1/theta)*exp(G[1,]),prob = rep(theta/(1 + theta),Q),log=TRUE))
      p2 <- sum(dnbinom(Cs[(N+i),],size = (1/theta)*exp(G[2,]),prob = rep(theta/(1 + theta),Q),log=TRUE))

      v.g <- p1 + p2
      L[i] <- v.g
  }}

  if(model == "CNBM"){
      for (i in 1:N){
      G[1,] <- FixEf(Xs[i,],Q, lvl.cov,b,Des)
      G[2,] <- FixEf(Xs[(N+i),],Q,lvl.cov,b,Des)

      p1 <- sum(dnbinom(Cs[i,],size = (1/theta)*exp(G[1,]), prob = rep((1/theta)/(1+(1/theta)),Q),log=TRUE))
      p2 <- sum(dnbinom(Cs[(N+i),],size = (1/theta)*exp(G[2,]), prob = rep((1/theta)/(1+(1/theta)),Q),log=TRUE))
      p3 <- dnbinom(S[i], size = (1/theta)*sum(exp(G[1,])),prob = (1/theta)/(1+(1/theta)),log = TRUE)
      p4 <- dnbinom(S[(N+i)], size = (1/theta)*sum(exp(G[2,])),prob = (1/theta)/(1+(1/theta)),log = TRUE)

      v.g <- (p1 + p2) - (p3+p4)
      L[i] <- v.g}
  }
  res <- sum(L)
  print(res)
  return(res)
}
