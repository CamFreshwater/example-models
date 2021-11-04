FFcs <- function(i,b, theta, Cs, Xs, Q, lvl.cov, n.cov, Des){ ## cross-sectional setting
  x.vec <- Xs[i,]
  Cs.i <- Cs[i,]
  G <- FixEf(x.vec,Q,lvl.cov,b,Des)
  eta <- exp(G)
  Eta <- (1/theta)*c(eta)
  val<- CNB(Cs.i,Eta)
  return(val)}
