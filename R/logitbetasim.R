#' Simulate logitBeta
#'
#' This code simulates logit beta random variables
#' @param alpha A n-dimensional vector of shape parameters
#' @param kappa A second n-dimensional vector of shape parameters. kappa must be greater than alpha.
#' @return W An n-dimensional vector of simulated values
#' @export
logitbetasim<-function(alpha,kappa){

  temp1=rbeta(length(alpha),alpha,kappa-alpha)

  indics = temp1<=1e-12
  indics2 = (1-temp1)<=1e-12

  W = log(temp1/(1-temp1))

  bet2 = kappa-alpha

  #small shape and scale stuff
  if (sum(indics)>0){
    X1 = -(1/alpha[indics==1])*rgamma(sum(indics),matrix(1,sum(indics),1),matrix(1,sum(indics),1))
    Y1 = rgamma(sum(indics),bet2[indics==1],matrix(1,sum(indics),1))
    W[indics==1] = X1 - log(exp(X1)+Y1);
  }

  if (sum(indics2)>0){
    X1 = -(1/bet2[indics2==1])*rgamma(sum(indics2),matrix(1,sum(indics2),1),matrix(1,sum(indics2),1))
    Y1 = rgamma(sum(indics2),alpha[indics2==1],matrix(1,sum(indics2),1))
    W[indics2==1] = -X1 + log(exp(X1)+Y1);
  }

  return(W)
}
