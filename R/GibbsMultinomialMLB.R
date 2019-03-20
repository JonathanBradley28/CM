#' Multinomial LCM
#'
#' This code implements a version of the collapsed Gibbs sampler from Bradley et al. (2019). The model assumes that data follow a multinomial distribution with a logit link to a mixed effects model. The priors and random effects are assumed to follow a multivariate logit-beta distribution.
#' @param report Option to print the iteration number of the MCMC.
#' @param B The number of iterations of the collapsed Gibbs sampler.
#' @param y An n dimensional vector consisting of multinomial count-valued observations.
#' @param nn An n dimensional vector consisting of totals associated with multinomial count-valued observations.
#' @param X An nxp matrix of covariates. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param Psi An nxr matrix of basis functions. Each column represents a basis function. Each row corresponds to one of the n replicates.
#' @param nslice The burnin for the slice sampler
#' @import actuar MfUSampler stats devtools roxygen2 coda MASS
#' @return beta A pxB matrix of MCMC values for beta. The vector beta corresponds to the covariates X.
#' @return eta A rxB matrix of MCMC values for eta. The vector eta corresponds to the basis functions G.
#' @return xi A nxB matrix of MCMC values for xi. The vector xi corresponds to uncorrelated random effects.
#' @return alphab A 1xB matrix of MCMC values for the shape parameter associated with beta.
#' @return alphae A 1xB matrix of MCMC values for the shape parameter associated with eta.
#' @return alphax A 1xB matrix of MCMC values for the shape parameter associated with xi.
#' @return kappab A 1xB matrix of MCMC values for the second shape parameter associated with beta.
#' @return kappae A 1xB matrix of MCMC values for the second shape parameter associated with eta.
#' @return kappax A 1xB matrix of MCMC values for the second shape parameter associated with xi
#' @examples
#' library(CM)
#' set.seed(23)
#'
#' #simulate pseudo data
#' n = 4*50
#' X = cbind(1,sin(pi*seq(0,2,by = ((2 - 0)/(n - 1)))))
#' X=X[1:n,]
#' p = dim(X)[2]
#' output<-qr(X)
#' Q = qr.Q(output)
#'
#' #Moran's I basis functions (Hughes and Haran, 2013)
#' r = n/2
#' Psi = diag(n) - Q%*%t(Q)
#' output2<-eigen(Psi)
#' Psi = output2$vectors[,1:r]
#'
#' #large scale parameters
#' beta = c(0.01,-2)
#' mu = X%*%beta
#' plot(mu)
#'
#' #small scale parameters
#' eta = sqrt(var(X%*%beta)/2)*rnorm(r)
#' rand.effects = Psi%*%eta
#' plot(rand.effects)
#'
#' #fine scale parameters
#' xi = sqrt(var(Psi%*%eta)/2)*rnorm(n)
#'
#' #simulate multinomial data with 5 categories
#' nu = mu+rand.effects+xi
#' plot(nu)
#'
#' # organize pi according to Section 2.1
#' pi1 = exp(nu[1:(n/4)])/(1+exp(nu[1:(n/4)]))
#' p2 = exp(nu[(n/4 +1):(2*n/4)])/(1+exp(nu[(n/4 +1):(2*n/4)]))
#' pi2 = (1-pi1)*p2
#' p3 = exp(nu[(2*n/4 +1):(3*n/4)])/(1+exp(nu[(2*n/4 +1):(3*n/4)]))
#' pi3 = (1-pi1 - pi2)*p3
#' p4 = exp(nu[(3*n/4 +1):(4*n/4)])/(1+exp(nu[(3*n/4 +1):(4*n/4)]))
#' pi4 = (1-pi1 - pi2 - pi3)*p4
#' pi5 = (1-pi1 - pi2 - pi3-pi4)
#'
#' #generate sample sizes
#' m = rpois(n/4,50)
#' data = matrix(0,n/4,5)
#' for (j in 1:(n/4)){
#'   data[j,]=rmultinom(1, m[j], prob = c(pi1[j],pi2[j],pi3[j],pi4[j],pi5[j]))
#' }
#'
#' #organize vectors according to Section 2.1 of Bradley (2019)
#' y = matrix(data[,1:4],n,1)
#' nn = matrix(cbind(m,m-data[,1],m-rowSums(data[,1:2]),m-rowSums(data[,1:3])),n,1)
#' indics = nn>0
#' y = y[indics==1]
#' nn = nn[indics==1]
#' XX = X
#' PsiPsi = Psi
#' X = X[indics==1,]
#' Psi = Psi[indics==1,]
#'
#' #number of MCMC reps
#' B = 2000
#'
#' #Implement Pseudo-Code 2 with updates for shape parameters
#' output<-GibbsMultinomialMLB(100,B,y,nn,X,Psi)
#'
#' #check some trace plots
#' plot((as.mcmc(output$beta[1,1000:B])))
#' plot((as.mcmc(output$beta[2,1000:B])))
#' plot(as.mcmc(output$alphab[1,1000:B]))
#' plot(as.mcmc(output$alphae[1,1000:B]))
#' plot(as.mcmc(output$alphax[1,1000:B]))
#' plot(as.mcmc(output$kappab[1,1000:B]))
#' plot(as.mcmc(output$kappae[1,1000:B]))
#' plot(as.mcmc(output$kappax[1,1000:B]))
#' plot((as.mcmc(output$eta[5,1000:B])))
#' plot((as.mcmc(output$eta[15,1000:B])))
#' plot((as.mcmc(output$eta[20,1000:B])))
#' plot((as.mcmc(output$xi[10,1000:B])))
#' plot((as.mcmc(output$xi[20,1000:B])))
#' plot((as.mcmc(output$xi[30,1000:B])))
#'
#' #hows our estimate of beta? true beta is (0.01,-2)
#' apply(output$beta[,1000:B], 1, quantile,probs=c(0.025,0.5,0.975))
#'
#' #hows our estimate of the proportions?
#' #get mcmc reps of nu
#' nu.gibbs = XX%*%(output$beta[,1000:B])+PsiPsi%*%(output$eta[,1000:B])
#' nu.gibbs[indics==1,]=output$xi[,1000:B] +XX[indics==1,]%*%(output$beta[,1000:B])+PsiPsi[indics==1,]%*%(output$eta[,1000:B])
#' nu.gibbs[indics==0,]=XX[indics==0,]%*%(output$beta[,1000:B])+PsiPsi[indics==0,]%*%(output$eta[,1000:B])
#'
#' #reorganize mcmc reps of nu into pi according to Section 2.1
#' pi1m = exp(nu.gibbs[1:(n/4),])/(1+exp(nu.gibbs[1:(n/4),]))
#' p2m = exp(nu.gibbs[(n/4 +1):(2*n/4),])/(1+exp(nu.gibbs[(n/4 +1):(2*n/4),]))
#' pi2m = (1-pi1m)*p2m
#' p3m = exp(nu.gibbs[(2*n/4 +1):(3*n/4),])/(1+exp(nu.gibbs[(2*n/4 +1):(3*n/4),]))
#' pi3m = (1-pi1m - pi2m)*p3m
#' p4m = exp(nu.gibbs[(3*n/4 +1):(4*n/4),])/(1+exp(nu.gibbs[(3*n/4 +1):(4*n/4),]))
#' pi4m = (1-pi1m - pi2m - pi3)*p4m
#'
#' pit = c(pi1,pi2,pi3,pi4)
#' pihat = c(rowMeans(pi1m),rowMeans(pi2m),rowMeans(pi3m),rowMeans(pi4m))
#' pibounds = apply(rbind(pi1m,pi2m,pi3m,pi4m), 1, quantile, probs = c(0.025,0.975))
#' pilow = as.matrix(pibounds[1,],200,1)
#' piup = as.matrix(pibounds[2,],200,1)
#' pilow[pilow<0] = 0
#' piup[piup<0] = 0
#' pihat[pihat<0] = 0
#'
#' #figure
#' plot(pit,pihat,ylab=c("Estimated Proportion"),xlab=c("True Proportion"),main=c("Estimated Versus Truth"))
#' abline(0,1)
#'
#' #version of figure from manuscript
#' plot(1:n,pit, type = "l", col = "black",xlab = "Arbitrary Ordering of Proportions",ylab = "Proportions")
#' lines(1:n,pihat, type = "l", col = "red",lty=1)
#' legend("topright", c("Truth","Estimated"),col=c("black", "red"), lty=c(1,1))
#' @export
GibbsMultinomialMLB<-function(report,B,y,nn,X,Psi,nslice=2){

  p = dim(X)[2]
  r = dim(Psi)[2]
  n = length(y)

  #run the Gibbs sampler from Section 2.3
  beta.gibbs = matrix(0,p,B)
  eta.gibbs = matrix(0,r,B)
  xi.gibbs = matrix(0,n,B)
  alphab = matrix(1,1,B)
  kappab = matrix(1.1,1,B)
  alphax = matrix(1,1,B)
  kappax = matrix(1.1,1,B)

  alphae = matrix(1,1,B)
  kappae = matrix(1.1,1,B)

  simalphb=1
  simkappb = 1.1
  simalphe=1
  simkappe = 1.1
  simalphx=1
  simkappx = 1.1

  eps = 0.05
  rho = 0.9

  for (b in 2:B){

    #the adaptive rejection algorithm 'ars' requires several initial values. It is very difficult to select each initial
    #values within a Gibbs sampler that keeps values within the parameter space. Thus, occasionally an error will occur for this reason.
    # to ignore the error we use the "tryCatch" function, and set the update using ars to be the previous value in the chain


    #update beta
    alpha = c(t(rho*y+eps),t(rho*y+eps),alphab[b-1]*matrix(1,1,p))
    kappa = c(t(nn),t(nn),kappab[b-1]*matrix(1,1,p))
    H = rbind(X[,],0.000001*X[,],diag(p))
    W = logitbetasim(alpha,kappa) + rbind(-Psi[,]%*%eta.gibbs[,b-1] - xi.gibbs[,b-1],  matrix(0,length(y)+p,1))
    beta.gibbs[,b] = solve(t(H)%*%H)%*%t(H)%*%W


    #update shapes beta using adaptive rejection
    #update shapes beta using slice sampler
    f<-function(x){
      like= -p*(lgamma(x)) - p*(lgamma(kappab[b-1] -x))-x -x*sum(beta.gibbs[,b])
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }
      if (abs(kappab[b-1] -x)<0.01){
        like = -Inf;
      }
      if (kappab[b-1]<x){
        like = -Inf;
      }

      return(like)

    }
    alphab_temp<-MfU.Sample.Run(alphab[b-1], f, nsmp = nslice)
    alphab[b] = abs(alphab_temp[nslice])

    f<-function(x){

      like=p*(lgamma(x)) - p*(lgamma(x -alphab[b]))-x -x*sum(log(1+exp(beta.gibbs[,b])))
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }
      if (abs(x -alphab[b])<0.01){
        like = -Inf;
      }
      if (x<alphab[b]){
        like = -Inf;
      }

      return(like)
    }

    kappab_temp<-MfU.Sample.Run(kappab[b-1], f, nsmp = nslice)
    kappab[b] = abs(kappab_temp[nslice])

    #update eta
    alpha = c(t(rho*y+eps),t(rho*y+eps),alphae[b-1]*matrix(1,1,r))
    kappa = c(t(nn),t(nn),kappae[b-1]*matrix(1,1,r))
    H = rbind(Psi[,],0.000001*Psi[,],diag(r))
    W = logitbetasim(alpha,kappa)  + rbind(-X[,]%*%beta.gibbs[,b]-xi.gibbs[,b-1],matrix(0,length(y)+r,1))
    eta.gibbs[,b] = solve(t(H)%*%H)%*%t(H)%*%W

    #update shapes eta using adaptive rejection
    f<-function(x){

      like=(-r*(lgamma(x)) - r*(lgamma(kappae[b-1] -x))-x -x*sum(eta.gibbs[,b]))
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }
      if (abs(kappae[b-1] -x)<0.01){
        like = -Inf;
      }
      if (kappae[b-1]<x){
        like = -Inf;
      }

      return(like)

    }
    alphae_temp<-MfU.Sample.Run(alphae[b-1], f, nsmp = nslice)
    alphae[b] = abs(alphae_temp[nslice])

    f<-function(x){

      like =(r*(lgamma(x)) - r*(lgamma(x -alphae[b]))-x -x*sum(log(1+exp(eta.gibbs[,b]))))
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      if (abs(x -alphae[b])<0.01){
        like = -Inf;
      }


      if (x<alphae[b]){
        like = -Inf;
      }

      return(like)

    }
    kappae_temp<-MfU.Sample.Run(kappae[b-1], f, nsmp = nslice)
    kappae[b] = abs(kappae_temp[nslice])

    #update xi
    alpha = c(t(rho*y+eps),t(rho*y+eps),alphax[b-1]*matrix(1,1,length(y)))
    kappa = c(t(nn),t(nn),kappax[b-1]*matrix(1,1,length(y)))
    H = rbind(diag(length(y)),0.000001*diag(length(y)),diag(length(y)))
    W = logitbetasim(alpha,kappa)  + rbind(-X[,]%*%beta.gibbs[,b]-Psi[,]%*%eta.gibbs[,b],matrix(0,2*length(y),1))
    xi.gibbs[,b] = solve(t(H)%*%H)%*%t(H)%*%W

    if (sum(is.nan(xi.gibbs[,b]))>0){
      xi.gibbs[,b] = xi.gibbs[,b-1]
    }

    if (sum(is.infinite(xi.gibbs[,b]))>0){
      xi.gibbs[,b] = xi.gibbs[,b-1]
    }

    #update shape xi using adaptive rejection
    f<-function(x){
      like=(-n*(lgamma(x)) - n*(lgamma(kappax[b-1] -x))-x -x*sum(xi.gibbs[,b]))
      if (x<0.01){
        like = -Inf;
      }

      if (abs(kappax[b-1] -x)<0.01){
        like = -Inf;
      }


      if (x>20000){
        like = -Inf;
      }

      if (kappax[b-1]<x){
        like = -Inf;
      }

      return(like)

    }
    alphax_temp<-MfU.Sample.Run(alphax[b-1], f, nsmp = nslice)
    alphax[b] =abs(alphax_temp[nslice])

    f<-function(x){
      like = (n*(lgamma(x)) - n*(lgamma(x - alphax[b]))-x -x*sum(log(1+exp(xi.gibbs[,b]))))
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      if (abs(x - alphax[b])<0.01){
        like = -Inf;
      }

      if (x< alphax[b]){
        like = -Inf;
      }

      return(like)

    }

    kappax_temp<-MfU.Sample.Run(kappax[b-1], f, nsmp = nslice)
    kappax[b] = abs(kappax_temp[nslice])

    if (b>(report-1)){
    if ((max(b,report) %% min(b,report)) ==0){
      print(c('iteration',b))
    }
    }

  }

  output<-list(beta=beta.gibbs,eta=eta.gibbs,xi=xi.gibbs,alphab=alphab,kappab=kappab,alphae=alphae,kappae=kappae,alphax=alphax,kappax=kappax)
  return(output)
}
