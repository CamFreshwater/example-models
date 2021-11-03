#' A loglikelihood value of a subject given a random effect \code{z}
#'
#' @param z A random effect value.
#' @param u A standard deviation of random effect.
#' @param G Fixed effect regression.
#' @param theta An overdispersion parameter.
#' @param Cs.vec A vector of categorical covariate for a specific subject.
#'



FF <- function(z,u,G,theta,Cs.vec){
  Q <- ncol(G)
  Eta1 <- (1/theta)*c(exp(G[1,] + rep(z,Q)))
  Eta2 <- (1/theta)*c(exp(G[2,] + rep(z,Q)))
  val <- CNB(Cs.vec[1,],Eta1) + CNB(Cs.vec[2,],Eta2) + dnorm(z,mean=0,sd = u,log=TRUE)
  return(val)}
