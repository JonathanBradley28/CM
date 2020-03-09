#' LGP
#'
#' This code implements a standard Gibbs sampler for normal data using a mixed effects model
#' @param B The number of iterations of the Gibbs sampler
#' @param X An nxp matrix of covariates.
#' @param G An nxr matrix of basis functions.
#' @param Z An n dimensional vector data.
#' @param betas2 An optional argument to define the initialized value for the fixed effects.
#' @param etas2 An optional argument to define the initialized value for the correlated random effects.
#' @param xi An optional argument to define the initialized value for the uncorrelated random effects.
#' @param taus2 An optional argument to define the initialized value for the variance of the data.
#' @param Hphib2 An optional argument to define the initialized value for the covariance of the random effects.
#' @param sig2xi An optional argument to define the initialized value for the variance of the random effects.
#' @param printevery Option to print the iteration number of the MCMC.
#' @import actuar MfUSampler stats coda MASS progress Matrix LearnBayes
#' @return ans A list of updated parameter values.
#' @examples
#'
#'#load the necessary packages
#'library(CM)
#'
#'set.seed(123)
#'#define a test function
#'#A non-linear test function
#'lambda <- function(t) exp(1.1 + sin(2*pi*t))
#'
#'#define some 1-d locations
#'points = seq(0,1,length.out=1001)
#'points=points[2:1001]
#'m = dim(as.matrix(points))[1]
#'
#'#get the true mean at these locations
#'truemean<-matrix(0,m,1)
#'for (j in 1:length(points)){
#'  truemean[j,1] = lambda(points[j])
#'}
#'
#'#simulate the data
#'data = matrix(0,m,1)
#'for (i in 1:m){
#'  data[i] = rnorm(1,truemean[i],1)
#'}
#'
#'#plot the data
#'plot(data,xlab="time",ylab="Poisson counts",main="Response vs. time")
#'
#'#covariate intercept-only
#'X = matrix(1,m,1)
#'p <- dim(X)[2]
#'
#'##compute the basis function matrix
#'#compute thin-plate splines
#'r = 8
#'knots = seq(0,1,length.out=r)
#'
#'#orthogonalize G
#'G = THINSp(as.matrix(points,m,1),as.matrix(knots,r,1))
#'outG<-qr(G)
#'G<-qr.Q(outG)
#'
#'#orthogonalize X
#'outX<-qr(X)
#'X<-qr.Q(outX)
#'
#'#Run the MCMC algorithm
#'output<-GibbsNormalGAU(2000,X,G,Z=data)
#'
#'#trace plots (without burnin)
#'plot(as.mcmc(output$betas[1000:2000]))
#'plot(as.mcmc(output$etas[1,1000:2000]))
#'plot(as.mcmc(output$etas[8,1000:2000]))
#'plot(as.mcmc(output$xis[10,1000:2000]))
#'#estimates (remove a burnin)
#'f_est = apply(X%*%output$betas[,1000:2000]+G%*%output$etas[,1000:2000],1,mean)
#'f_lower= apply(X%*%output$betas[,1000:2000]+G%*%output$etas[,1000:2000],1,quantile,0.025)
#'f_upper= apply(X%*%output$betas[,1000:2000]+G%*%output$etas[,1000:2000],1,quantile,0.975)
#'
#'#plot estimates and truth
#'plot(1:m,truemean,ylim = c(0,max(f_upper)+1))
#'lines(1:m,f_est,col="red")
#'lines(1:m,f_lower,col="blue")
#'lines(1:m,f_upper,col="blue")
#' @export
GibbsNormalGAU<-function(B,X,G,Z,etas2=NA,betas2=NA,xi=NA,taus2=NA,Hphib2=NA,sig2xi=NA,sig2beta2=NA,printevery = 100){

p = dim(X)[2]
r = dim(G)[2]
n = length(Z)


if (is.na(betas2[1])==1){
storeit= 1
betas2 = mvrnorm(1, matrix(0,p,1), diag(p), tol = 1e-3)
etas2 = mvrnorm(1, matrix(0,r,1), diag(r), tol = 1e-3)
xi= mvrnorm(1, matrix(0,n,1), diag(n), tol = 1e-3)
taus2 = 1/rgamma(1,1)
Hphib2 = (1/rgamma(1,1))*diag(r)
sig2xi=1/rgamma(1,1)
sig2beta2=1
}


betas3 = betas2
etas3= etas2
xi3 = xi
taus3=taus2
Hphib3 = Hphib2[1,1]
sig2xi3= sig2xi
sig2beta3 = sig2beta2

  #Gibbs Sampling
  for (i in 2:B){


    #update betas
    BetAcovar<-solve((1/taus2)*t(X)%*%X + (1/sig2beta2)*diag(p))
    BetAmean<-(1/taus2)*BetAcovar%*%(t(X)%*%(Z-G%*%etas2-xi))
    betas2=mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)

betas3 = cbind(betas3,betas2)

    #update etas
    BetAcovar<-solve((1/taus2)*t(G)%*%G + solve(Hphib2))
    BetAmean<-(1/taus2)*BetAcovar%*%(t(G)%*%(Z-X%*%betas2-xi))
    etas2=mvrnorm(1, BetAmean, BetAcovar, tol = 1e-3)

etas3 = cbind(etas3,etas2)


    #update xi
    BetAcovar<-1/((1/taus2) + 1/sig2xi)
    BetAmean<-(1/taus2)*BetAcovar*((Z-X%*%betas2-G%*%etas2))
    xi=BetAmean+sqrt(BetAcovar)*rnorm(n)


xi3 = cbind(xi3 ,xi)


    #update sig2xi
    alphaTau =1+n/2
    betaTau = 1 + 0.5*(t(xi)%*%(xi))
    sig2xi=rigamma(1,alphaTau,betaTau)


sig2xi3= cbind(sig2xi3 ,sig2xi)


    #update taus
    alphaTau =1+n/2
    betaTau = 1 + 0.5*(t(Z-X%*%betas2-G%*%etas2-xi)%*%(Z-X%*%betas2-G%*%etas2-xi))
    taus2=rigamma(1,alphaTau,betaTau)


taus3= cbind(taus3,taus2)


    #update sigma2s
    alphaTau =1+r/2
    betaTau = 1 + 0.5*(t(etas2)%*%(etas2))
    sigtemp2=rigamma(1,alphaTau,betaTau)
     Hphib2=sigtemp2*diag(r)


Hphib3= cbind(Hphib3,sigtemp2)

    alphaTau =1+p/2
    betaTau = 1 + 0.5*(t(betas2)%*%(betas2))
    sigtemp2=rigamma(1,alphaTau,betaTau)
     sigbeta2=sigtemp2*diag(p)
sig2beta3 = cbind(sig2beta3 ,sigtemp2)

if (i%%printevery==0){
print(c("Iteration Number:", i))
}

}



ans<-list(as.matrix(betas3),as.matrix(etas3),as.matrix(xi3),as.matrix(taus3),as.matrix(Hphib3),as.matrix(sig2xi3),as.matrix(sig2beta3))

names(ans) <- c("betas","etas","xis","sig2eps","sig2eta","sig2xi","sig2beta")

return(ans)
}



