#' Moran's I Basis Function
#'
#' This code computes the Moran's I basis functions
#' @param X An nxp matrix of covariates. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param r The number of basis functions.
#' @param A An nxn adjacency matrix.
#' @return Psi nxr matrix of basis functions
#' @export
MoransI.Basis<-function(X,r,A){


  n = dim(X)[1]

  PsiOper = (diag(n) - X%*%solve(t(X)%*%X)%*%t(X))%*%A%*%(diag(n) - X%*%solve(t(X)%*%X)%*%t(X))
  output2<-eigen(PsiOper)
  Psi = output2$vectors[,1:r]

  return(Psi)
}
