#' Binomial LCM
#'
#' This code implements a version of the collapsed Gibbs sampler from Bradley et al. (2018). The model assumes that data follow a binomial distribution with a logit link to a mixed effects model. The priors and random effects are assumed to follow a multivariate logit-beta distribution.
#' @param Niter The number of iterations of the collapsed Gibbs sampler
#' @param X An nxp matrix of covariates. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param G An nxr matrix of basis functions. Each column represents a basis function. Each row corresponds to one of the n replicates.
#' @param data An n dimensional vector consisting of count-valued observations.
#' @param data An n dimensional vector consisting of totals associated with binomial count-valued observations.
#' @param report Option to print the iteration number of the MCMC.
#' @param rho A parameter between zero and one to allow for zero counts
#' @param eps A second parameter between zero and one to allow for zero counts
#' @param cs A three dimensional parameter vector to allow prior to have non-zero means.
#' @param nslice The burnin for the slice sampler
#' @param select.h Option to select pseudo-optimal values for rho and eps. The MCMC is run for only 50 iterations, and the corresponding predictions are compared to a holdout datasets. The values of rho and eps are chosen to minimize the Hellinger distance between the holdout and the naive predictions using the function optim.
#' @import actuar MfUSampler stats devtools roxygen2 coda MASS
#' @return beta A pxNiter matrix of MCMC values for beta. The vector beta corresponds to the covariates X.
#' @return eta A rxNiter matrix of MCMC values for eta. The vector eta corresponds to the basis functions G.
#' @return xi A nxNiter matrix of MCMC values for xi. The vector xi corresponds to uncorrelated random effects.
#' @return pi A nxNiter matrix of MCMC values for the proportions of the binomial data, modeled with, exp(Xbeta+Geta+xi)/(1+exp(Xbeta+Geta+xi)).
#' @return pi.smooth A nxNiter matrix of MCMC values for the smoothed proportions associated with binomial data, modeled with, exp(Xbeta+Geta)/(1+exp(Xbeta+Geta)).
#' @return alphab A 1xNiter matrix of MCMC values for the shape parameter associated with beta.
#' @return alphae A 1xNiter matrix of MCMC values for the shape parameter associated with eta.
#' @return alphax A 1xNiter matrix of MCMC values for the shape parameter associated with xi.
#' @return kappab A 1xNiter matrix of MCMC values for the second shape parameter associated with beta.
#' @return kappae A 1xNiter matrix of MCMC values for the second shape parameter associated with eta.
#' @return kappax A 1xNiter matrix of MCMC values for the second shape parameter associated with xi
#' @examples
#' library(CM)
#'
#' set.seed(2)
#'
#' #simulate pseudo data
#' n = 1000
#' X = cbind(1,sin(pi*seq(0,2,by = ((1 - 0)/(n - 1)))))
#' X=X[1:n,]
#' p = dim(X)[2]
#' output<-qr(X)
#' Q = qr.Q(output)
#'
#' #Moran's I basis functions (Hughes and Haran, 2013)
#' r = 100
#' Psi = MoransI.Basis(Q,r,diag(n))
#'
#' #large scale parameters
#' beta = c(0.01,-2)
#' mu = X%*%beta
#' plot(mu)
#'
#' #small scale parameters
#' eta = sqrt(as.numeric(var(X%*%beta)/2))*rnorm(r)
#' rand.effects = Psi%*%eta
#' plot(rand.effects)
#'
#' #fine scale parameters
#' xi = sqrt(as.numeric(var(Psi%*%eta)/2))*rnorm(n)
#'
#' #simulate binomial data
#' nu = mu+rand.effects+xi
#' plot(nu)
#'
#' #inverse logit
#' pi1 = exp(nu[1:n])/(1+exp(nu[1:n]))
#'
#' #generate sample sizes
#' m = rpois(n,50)
#' data = matrix(0,n,1)
#' for (j in 1:n){
#'   data[j]=rmultinom(1, m[j], prob = c(pi1[j],1-pi1[j]))[1]
#' }
#'
#' plot(data)
#'
#' #number of MCMC reps
#' B = 2000
#' burnin=1000
#'
#' #Pseudo-Code 2 with updates for shape parameters
#' output<-GibbsBinomialMLB(B,data,m,X,Psi)
#'
#' #check some trace plots
#' plot((as.mcmc(output$beta[1,burnin:B])))
#' plot((as.mcmc(output$beta[2,burnin:B])))
#' plot(as.mcmc(output$alphab[1,burnin:B]))
#' plot(as.mcmc(output$alphae[1,burnin:B]))
#' plot(as.mcmc(output$alphax[1,burnin:B]))
#' plot(as.mcmc(output$kappab[1,burnin:B]))
#' plot(as.mcmc(output$kappae[1,burnin:B]))
#' plot((as.mcmc(output$eta[5,burnin:B])))
#' plot((as.mcmc(output$xi[10,burnin:B])))
#'
#' #hows our estimate of beta? true beta is (0.01,-2)
#' apply(output$beta[,burnin:B], 1, quantile,probs=c(0.025,0.5,0.975))
#'
#'#hows our estimate of the proportions?
#'pihat = rowMeans(output$pi)
#'pibounds = apply(output$pi, 1, quantile, probs = c(0.025,0.975))
#'pilow = as.matrix(pibounds[1,],n,1)
#'piup = as.matrix(pibounds[2,],n,1)
#'
#' #figure
#' plot(pi1,pihat,ylab=c("Estimate Proportion"),xlab=c("True Proportion"),main=c("Estimated Versus Truth"))
#' abline(0,1)
#'
#figure
#' plot(1:n,pi1, type = "l", col = "black",xlab = "Arbitrary Ordering of Proportions",ylab = "Proportions")
#' lines(1:n,pihat, type = "l", col = "red",lty=1)
#' legend("topright", c("Truth","Estimated"),col=c("black", "red"), lty=c(1,1))
#'
#' @export
GibbsBinomialMLB<-function(B,data,nn,X,G,report=100,rho=0.999,cs=c(0,0,0),eps=0.01,nslice=2,select.h=0){

  n = length(data)
  p = dim(X)[2]
  r = dim(G)[2]
  cbeta=cs[1]
  ceta=cs[2]
  cxi=cs[3]

  if (select.h==1){

    ind=sample(n,n)
    train=data[ind[(floor(0.1*n)+1):(floor(0.5*n))]]
    valid=data[ind[1:floor(0.1*n)]]

    Hbeta = rbind(X[ind[(floor(0.1*n)+1):(floor(0.5*n))],],0.000001*X[ind[(floor(0.1*n)+1):(floor(0.5*n))],],diag(p))
    HpHinvBeta = solve(t(Hbeta)%*%Hbeta)

    Heta = rbind(G[ind[(floor(0.1*n)+1):(floor(0.5*n))],],0.000001*G[ind[(floor(0.1*n)+1):(floor(0.5*n))],],diag(r))
    HpHinvEta = solve(t(Heta)%*%Heta)
    H = rbind(diag(length(train)),0.000001*diag(length(train)),diag(length(train)))

    GibbsBinomialMLB2<-function(B,data,nn,X,G,report=100,rho=0.999,eps=0.01,cb,ce,cx,nslice){

      p = dim(X)[2]
      r = dim(G)[2]
      #run the Gibbs sampler from Section 2.3
      n = length(data)
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

      for (b in 2:B){

        #the adaptive rejection algorithm 'ars' requires several initial values. It is very difficult to select each initial
        #values within a Gibbs sampler that keeps values within the parameter space. Thus, occasionally an error will occur for this reason.
        # to ignore the error we use the "tryCatch" function, and set the update using ars to be the previous value in the chain


        #update beta
        alpha = c(t(rho*data+eps),t(rho*data+eps),alphab[b-1]*matrix(1,1,p))
        kappa = c(t(nn),t(nn),kappab[b-1]*matrix(1,1,p))
        W = logitbetasim(alpha,kappa) + rbind(-G%*%eta.gibbs[,b-1] - xi.gibbs[,b-1],  matrix(-cb,length(data)+p,1))
        beta.gibbs[,b] = HpHinvBeta%*%t(Hbeta)%*%W


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
        alpha = c(t(rho*data+eps),t(rho*data+eps),alphae[b-1]*matrix(1,1,r))
        kappa = c(t(nn),t(nn),kappae[b-1]*matrix(1,1,r))
        W = logitbetasim(alpha,kappa)  + rbind(-X%*%beta.gibbs[,b]-xi.gibbs[,b-1],matrix(-ce,length(data)+r,1))
        eta.gibbs[,b] = HpHinvEta%*%t(Heta)%*%W

        #update shapes eta using slice sampler
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
        alpha = c(t(rho*data+eps),t(rho*data+eps),alphax[b-1]*matrix(1,1,length(data)))
        kappa = c(t(nn),t(nn),kappax[b-1]*matrix(1,1,length(data)))
        W = logitbetasim(alpha,kappa)  + rbind(-X%*%beta.gibbs[,b]-G%*%eta.gibbs[,b],matrix(-cx,2*length(data),1))
        xi.gibbs[,b] = (1/3)*t(H)%*%W

        if (sum(is.nan(xi.gibbs[,b]))>0){
          xi.gibbs[,b] = xi.gibbs[,b-1]
        }

        if (sum(is.infinite(xi.gibbs[,b]))>0){
          xi.gibbs[,b] = xi.gibbs[,b-1]
        }

        #update shape xi using slice sampler
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
            print(paste("MCMC Replicate: ",b))
          }
        }

      }

      nu.gibbs = X%*%(beta.gibbs)+G%*%(eta.gibbs) + xi.gibbs

      #reorganize mcmc reps of nu into pi according to Section 2.1
      pi.gibbs = exp(nu.gibbs)/(1+exp(nu.gibbs))

      output<-list(pi=pi.gibbs,beta=beta.gibbs,eta=eta.gibbs,xi=xi.gibbs,alphab=alphab,kappab=kappab,alphae=alphae,kappae=kappae,alphax=alphax,kappax=kappax)
      return(output)
    }


    SimpleClassifier<-function(data,phat,holdout=NULL,phathold=NULL){

      cutoff<-function(lambda){
        estpos=phat>lambda
        helliginer= sum((sqrt(data)-sqrt(estpos))^2)
        return(helliginer)
      }

      initB=0.1
      estpar<-optim(initB,cutoff,method="SANN")
      estpos.final=phat>estpar$par[1]

      helliginer=NULL
      estpos.hold=NULL
      if(length(holdout)>0){

        estpos.hold=phathold>estpar$par[1]
        helliginer= sum((sqrt(holdout)-sqrt(estpos.hold))^2)
      }

      output<-list(helliginer,estpar$par[1],estpos.final,estpos.hold)

      return(output)
    }


    fz<-function(z){

      hellingerDist=Inf

      if (z[1]<0.99){
        if(0.01<z[1]){
          if(z[2]<0.99){
            if (0.01<z[2]){
              if (z[1]+z[2]<0.99){

                output<-GibbsBinomialMLB2(50,train,nn[ind[(floor(0.1*n)+1):(floor(0.5*n))]],as.matrix(X[ind[(floor(0.1*n)+1):(floor(0.5*n))],]),G[ind[(floor(0.1*n)+1):(floor(0.5*n))],],report=1e15,rho=z[1],eps=z[2],z[3],z[4],z[5],2)

                nuest=X[ind[(floor(0.1*n)+1):(floor(0.5*n))],]%*%output$beta + G[ind[(floor(0.1*n)+1):(floor(0.5*n))],]%*%output$eta +output$xi
                nuesthold=X[ind[1:floor(0.1*n)],]%*%output$beta + G[ind[1:floor(0.1*n)],]%*%output$eta
                for (bb in 1:50){
                  nuesthold[,bb]=nuesthold[,bb]+logitbetasim(matrix(output$alphax[bb],length(valid),1),matrix(output$kappax[bb],length(valid),1))
                }
                piest = exp(nuest)/(1+exp(nuest))
                phat = rowMeans(piest[,25:50])

                piesthold = exp(nuesthold)/(1+exp(nuesthold))
                phathold = rowMeans(piesthold[,25:50])

                results<-SimpleClassifier(train,phat,valid,phathold)


                hellingerDist = results[[1]]
              }
            }
          }
        }
      }

      return(hellingerDist)
    }

    initB=c(0.5,0.15,-3.5,-3.5,-3.5)
    print("Selecting hyperparameters")
    estpar<-optim(initB,fz,method = "Nelder-Mead",control=list(maxit=300))
    rho=estpar$par[1]
    eps=estpar$par[2]
    cbeta=estpar$par[3]
    ceta=estpar$par[4]
    cxi=estpar$par[5]
  }

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

  Hbeta = rbind(X,0.000001*X,diag(p))
  HpHinvBeta = solve(t(Hbeta)%*%Hbeta)

  Heta = rbind(G,0.000001*G,diag(r))
  HpHinvEta = solve(t(Heta)%*%Heta)
  H = rbind(diag(length(data)),0.000001*diag(length(data)),diag(length(data)))

  print("Collapsed Gibbs Sampler Starting")
  for (b in 2:B){

    #the adaptive rejection algorithm 'ars' requires several initial values. It is very difficult to select each initial
    #values within a Gibbs sampler that keeps values within the parameter space. Thus, occasionally an error will occur for this reason.
    # to ignore the error we use the "tryCatch" function, and set the update using ars to be the previous value in the chain


    #update beta
    alpha = c(t(rho*data+eps),t(rho*data+eps),alphab[b-1]*matrix(1,1,p))
    kappa = c(t(nn),t(nn),kappab[b-1]*matrix(1,1,p))
    W = logitbetasim(alpha,kappa) + rbind(-G%*%eta.gibbs[,b-1] - xi.gibbs[,b-1],  matrix(-cbeta,length(data)+p,1))
    beta.gibbs[,b] = HpHinvBeta%*%t(Hbeta)%*%W


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
    alpha = c(t(rho*data+eps),t(rho*data+eps),alphae[b-1]*matrix(1,1,r))
    kappa = c(t(nn),t(nn),kappae[b-1]*matrix(1,1,r))
    W = logitbetasim(alpha,kappa)  + rbind(-X%*%beta.gibbs[,b]-xi.gibbs[,b-1],matrix(-ceta,length(data)+r,1))
    eta.gibbs[,b] = HpHinvEta%*%t(Heta)%*%W

    #update shapes eta using slice sampler
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
    alpha = c(t(rho*data+eps),t(rho*data+eps),alphax[b-1]*matrix(1,1,length(data)))
    kappa = c(t(nn),t(nn),kappax[b-1]*matrix(1,1,length(data)))
    W = logitbetasim(alpha,kappa)  + rbind(-X%*%beta.gibbs[,b]-G%*%eta.gibbs[,b],matrix(-cxi,2*length(data),1))
    xi.gibbs[,b] = (1/(2+0.000001^2))*t(H)%*%W

    if (sum(is.nan(xi.gibbs[,b]))>0){
      xi.gibbs[,b] = xi.gibbs[,b-1]
    }

    if (sum(is.infinite(xi.gibbs[,b]))>0){
      xi.gibbs[,b] = xi.gibbs[,b-1]
    }

    #update shape xi using slice sampler
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
      print(paste("MCMC Replicate: ",b))
    }
    }

  }

  nu.gibbs = X%*%(beta.gibbs)+G%*%(eta.gibbs) + xi.gibbs

  #reorganize mcmc reps of nu into pi according to Section 2.1
 pi.gibbs = exp(nu.gibbs)/(1+exp(nu.gibbs))

 nu.smooth.gibbs=X%*%(beta.gibbs)+G%*%(eta.gibbs)
 pi.smooth.gibbs= exp(nu.smooth.gibbs)/(1+exp(nu.smooth.gibbs))

 output<-list(pi=pi.gibbs,pi.smooth=pi.smooth.gibbs,beta=beta.gibbs,eta=eta.gibbs,xi=xi.gibbs,alphab=alphab,kappab=kappab,alphae=alphae,kappae=kappae,alphax=alphax,kappax=kappax)
 return(output)
}
