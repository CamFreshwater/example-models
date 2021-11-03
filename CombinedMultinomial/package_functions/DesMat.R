#' Design Matrix of a loglinear regression
#'
#'
#' @param Q A number of categories in a multivariate count outcome.
#' @param lvl.cov A list containing the level of each covariate.
#'
#' @return A list of two elements where the first element is the left hand side of loglinear regression and the second element is the design matrix of loglinear regression.
#' The first element is a matrix of all possible combination of multivariate count and categorical covariate levels. See vignette for more detail.




DesMat <- function(Q,lvl.cov){
#
  library(gtools)
  library(rlist)


  nvar <- length(lvl.cov) + 1
  lvl.cov <- append(Q,lvl.cov)
  names(lvl.cov) <- letters[1:length(lvl.cov)]

  lvl.cov1 <- unlist(lvl.cov)
  var.names <- NULL
  idx.names <- NULL

for (i in 1:nvar){
  Temp <- combinations(nvar,i)
  for (j in 1:nrow(Temp)){
    tp <- Temp[j,]
    if (length(tp) == 1) {
    idxused <- c(2:lvl.cov1[tp])
    var.names <- append(var.names,rep(letters[Temp[j,]],length(idxused)))
    idx.names <- append(idx.names,idxused)
    }else {
      tmp.idx <- sapply(lvl.cov1[tp],function(x) c(2:x))
      idxused <- expand.grid(tmp.idx)
        if(ncol(idxused) == 1){
          var.names <- append(var.names,list(letters[tp]))
          idx.names <- append(idx.names,list(as.numeric(tmp.idx)))
          } else {
          x <- as.matrix(idxused)
          xl <- tapply(x,rep(1:nrow(x),ncol(x)),function(i)i)
          var.names <- append(var.names,rep(list(letters[tp]),nrow(idxused)))
        idx.names <- append(idx.names, xl) }
    }
    }
    }


## Filling in the matrix entries

  lhs.out <- expand.grid(lapply(lvl.cov,function(i) c(1:i)))
  outMat <- lhs.out#order(lhs.out[,1]),]

  MM <- matrix(0,nrow = nrow(outMat), ncol = length(var.names)+1)
  MM[,1] <- 1



for (i in 2:((length(var.names))+1)){
  vname <- var.names[[i - 1]]
  idx <- idx.names[[i - 1]]

      idxcolsel <- match(vname,names(outMat))
      idx1 <- c(1:nrow(outMat))
      rep <- 1
      while(rep <= length(idxcolsel)){
        idx1 <- idx1[which(outMat[idx1,idxcolsel[rep]] == idx[rep])]
        rep <- rep+1}
    MM[idx1,i] <- 1
    #}
}

var.names.col <- lapply(var.names, function(x) paste0(x, collapse = ""))
idx.names.col <- lapply(idx.names, function(x) paste0(x, collapse = ""))

CN <- NULL
for (i in 1:length(var.names)){
  CN[i] <- paste0(var.names.col[[i]],"_",idx.names.col[[i]])
}

colnames(MM) <- c("l_0",CN)
return(list("LHS" = outMat, "Des.mat" = MM))
}
