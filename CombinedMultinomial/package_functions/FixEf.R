#' The Fixed Effect of a Regression on Count Data.
#' @param x.vec The observed vector of covariate for each category.
#' @param Q The number of categories in a multivariate count outcome.
#' @param lvl.cov A list containing the level of each covariate.
#' @param b A vector containing estimates of covariate effects.
#' @param Des A list containing two matrices, the first element is the left hand side of the equation. The second element is the design matrix.
#' Could be a result from \code{DesMat}.
#'
#' @return A \code{Q} dimensional vector consisting of the expected values of \code{Q} categories given the \code{x.vec}.
#'

FixEf <- function(x.vec,Q,lvl.cov,b,Des){

  Y.form <- Des[[1]]
  D <- Des[[2]]

  gt <- c(D %*% b)

  idx.val <- x.vec
  idx1 <- c(1:nrow(Y.form))
  rep <- 1
  while(rep <= length(x.vec)){
    idx1 <- idx1[which(Y.form[idx1,(rep+1)] == idx.val[[rep]])]
    rep <- rep+1}

  return(gt[c(idx1)])
}

