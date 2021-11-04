#' Generating Dataset consisting of Multivariate Count Outcome
#'
#' @param N The number of subject in the generated dataset.
#' @param Q The number of categories in the multivariate count outcome.
#' @param S The total number of multivariate count outcome per subject at each time point.
#' @param theta The overdispersion parameter.
#' @param s.u The variance of random effect to account for the correlation within time point.
#' @param n.cov The number of covariates.
#' @param Des Design matrix for fixed effects.
#' @param method An argument of either \code{DMM} or \code{UNBM}
#'
#' @return A list of two elements where the first element is an \eqn{N\times Q} matrix of multivariate count outcome and the second element
#' is \code{N}$\times$ \code{n.cov} categorical covariates.
#' @examples See vignette


DataGen <- function(N,Q,S,beta,theta,s.u,n.cov, lvl.cov, Des, method, T = 2){
  require(HMP)
  U <- rnorm(N,mean=0,sd=s.u)
  X <- cbind(rbinom(2*N,size=1,prob = 0.5))+1
  XB1 <- t(sapply(X[1:N],FixEf,Q = Q, lvl.cov = lvl.cov, b = c(0,beta), Des = Des))
  XB2 <- t(sapply(X[(N+1):(2*N)],FixEf,Q = Q, lvl.cov = lvl.cov, b = c(0,beta), Des = Des))


  if(method == "DMM"){
    XB.tilde1 <- (1/theta)*exp(cbind(XB1+rep(U,Q)))
    C1 <- t(apply(XB.tilde1,1,Dirichlet.multinomial,Nrs = S))

    XB.tilde2 <- (1/theta)*exp(cbind(XB2+rep(U,Q)))
    C2 <- t(apply(XB.tilde2,1,Dirichlet.multinomial,Nrs = S))

    datalong <- data.frame(cbind(rep(1:N, 2),rbind(cbind(C1,X[1:N]),cbind(C2,X[(N+1):(2*N)]))))
    colnames(datalong) <- c("ID",paste0("C",1:Q), paste0("X",1:n.cov))
  }

  if (method == "UNBM"){
        XB.tilde1 <- (1/theta)*exp(log(S/10) + cbind(XB1+rep(U,Q)))
        XB.tilde2 <- (1/theta)*exp(log(S/10) + cbind(XB2+rep(U,Q)))

    C1 <- matrix(NA,nrow=N,ncol=Q)
    C2 <- matrix(NA,nrow=N,ncol=Q)

    for (i in 1:N)
    {for (j in 1:Q)
    {C1[i,j] <- rnbinom(1,size = XB.tilde1[i,j],prob = rep((1/theta)/(1+(1/theta)),Q))
    C2[i,j] <- rnbinom(1,size = XB.tilde2[i,j],prob = rep((1/theta)/(1+(1/theta)),Q))
    }
    }
    datalong <- data.frame(cbind(rep(1:N, 2),rbind(cbind(C1,X[1:N]),cbind(C2,X[(N+1):(2*N)]))))
    colnames(datalong) <- c("ID",paste0("C",1:Q), paste0("X",1:n.cov))
  }

  return(datalong)
}
