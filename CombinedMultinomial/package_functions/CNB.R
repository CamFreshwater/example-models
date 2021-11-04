#' The marginal loglikelihood of the multivariate count outcome.
#'
#' @param Ct A vector of categorical count.
#' @param g A function of expected values for each category of the multivariate count outcome.
#'
#' @return a value of marginal loglikelihood  of Dirichlet - multinomial distribution.
#' @examples
#' Ct <- c(179, 12,9)
#' g <- c(0.65, 0.25, 0.10)
#' CNB(Ct,g)

CNB <- function(Ct,g){
  gs <- sum(g)
  ys <- sum(Ct)
  val <- lgamma(ys+1) + lgamma(gs) - lgamma(ys+gs) + sum(lgamma(Ct+g) - lgamma(g) - lgamma(Ct+1))
  return(val)
}
