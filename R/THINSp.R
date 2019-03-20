#' Thin-Plate Spline
#'
#' This code computes the thin-plate spline function f(s) = s^2 log(s) for 1-dimensional spatial locations
#' @param locs n-dimensional vector of locations
#' @param knots r-dimensional vector of knots
#' @param tol thresholds small values of the elements of psi to be zero. Default is no threshold.
#' @examples
#' #example two dimensional separable thin-plate spline
#' points1 = seq(0,1,length.out=1001)
#' points1=points1[2:1001]
#' r = 10
#' knots = seq(0,1,length.out=r)
#' G1 = THINSp(as.matrix(points1,m,1),as.matrix(knots,r,1))
#'
#' G=G1 %x% G1
#'
#' @return psi nxr matrix of basis functions
#' @export
THINSp<-function(locs,knots,tol=0){

  r = dim(knots)[1]
  n = dim(locs)[1]

  psi = matrix(0,n,r)
  for (i in 1:n){
    for (j in 1:r){
      psi[i,j] = (abs(knots[j]-locs[i]))^2*log(abs(knots[j]-locs[i])+0.01)
   if (abs(psi[i,j])<tol){
     psi[i,j] = 0
   }
  }
  }

  return(psi)
}
