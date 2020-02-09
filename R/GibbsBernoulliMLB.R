#' Bernoulli LCM
#'
#' This code implements a version of the collapsed Gibbs sampler from Bradley et al. (2018). The model assumes that data follow a Bernoulli distribution with a logit link to a mixed effects model. The priors and random effects are assumed to follow a multivariate logit-beta distribution.
#' @param B The number of iterations of the collapsed Gibbs sampler
#' @param X An nxp matrix of covariates. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param G An nxr matrix of basis functions. Each column represents a basis function. Each row corresponds to one of the n replicates.
#' @param data An n dimensional vector consisting of totals associated with binomial count-valued observations.
#' @param itmax The values of rho and eps in "GibbsBinomialMLB" are chosen to minimize the Hellinger distance between the data and the estimated proportion. itmax defines the maximum number of iterations in "optim" to obtain these values for rho and eps.
#' @param report Option to print the iteration number of the MCMC.
#' @param nslice The burnin for the slice sampler
#' @import actuar MfUSampler stats coda MASS
#' @return beta A pxB matrix of MCMC values for beta. The vector beta corresponds to the covariates X.
#' @return eta A rxB matrix of MCMC values for eta. The vector eta corresponds to the basis functions G.
#' @return xi A nxB matrix of MCMC values for xi. The vector xi corresponds to uncorrelated random effects.
#' @return pi A nxB matrix of MCMC values for the proportions of the binomial data, modeled with, exp(Xbeta+Geta+xi)/(1+exp(Xbeta+Geta+xi)).
#' @return pi.smooth A nxB matrix of MCMC values for the smoothed proportions associated with binomial data, modeled with, exp(Xbeta+Geta)/(1+exp(Xbeta+Geta)).
#' @return alphab A 1xB matrix of MCMC values for the shape parameter associated with beta.
#' @return alphae A 1xB matrix of MCMC values for the shape parameter associated with eta.
#' @return alphax A 1xB matrix of MCMC values for the shape parameter associated with xi.
#' @return kappab A 1xB matrix of MCMC values for the second shape parameter associated with beta.
#' @return kappae A 1xB matrix of MCMC values for the second shape parameter associated with eta.
#' @return kappax A 1xB matrix of MCMC values for the second shape parameter associated with xi
#' @return hypparams The selected hyperparameters.
#' @examples
#'
#' library(CM)
#' 
#' set.seed(3000)
#' 
#' #simulate pseudo data
#' n = 400
#' #define a test function
#' #A non-linear test function
#' lambda <- function(t) (0.1 + sin(4*pi*t))
#' 
#' #define some 1-d locations
#' points = seq(0,1,length.out=n+2)
#' points=points[1:n+2]
#' 
#' #get the logit true proportion at these locations
#' truemean<-matrix(0,n,1)
#' for (j in 1:length(points)){
#'   truemean[j,1] = lambda(points[j])
#' }
#' 
#' #simulate logit proportion
#' plot(truemean)
#' 
#' #inverse logit
#' pi1 = exp(truemean)/(1+exp(truemean))
#' 
#' plot(pi1)
#' 
#' #generate sample sizes
#' m = matrix(1,n,1)
#' data = matrix(0,n,1)
#' for (j in 1:n){
#'   data[j]=rmultinom(1, m[j], prob = c(pi1[j],1-pi1[j]))[1]
#' }
#' 
#' plot(data)
#' 
#' #covariate intercept-only
#' X = matrix(1,n,1)
#' 
#' ##compute the basis function matrix
#' #compute thin-plate splines
#' r = 10
#' knots = seq(0,1,length.out=r)
#' 
#' # G
#' Bisquare<-function(locs,knots,tol,smoother){
#'   r = dim(knots)[1]
#'   n = dim(locs)[1]
#'   disKnots = matrix(0,r,r)
#'   for (j in 1:r){
#'     disKnots[r,]=sqrt((knots[j,]-knots)^2)
#'   }
#'   disKnots<-disKnots[disKnots!=0]
#'   r_l<-smoother*max(disKnots)
#'   
#'   psi = matrix(0,n,r)
#'   for (i in 1:n){
#'     for (j in 1:r){
#'       if (sum(sqrt((locs[i,]-knots[j,])^2))<tol){
#'         psi[i,j]= (1-(sum(sqrt((locs[i,]-knots[j,])^2))/r_l)^2)^2
#'       }
#'     }
#'   }
#'   for (j in 1:r){
#'     psi[,j] = psi[,j]/max(psi[,j])
#'   }
#'   return(psi)
#' }
#' G = Bisquare(as.matrix(points,n,1),as.matrix(knots,r,1),tol=0.5,smoother=0.5)
#' 
#' #number of MCMC reps
#' B = 10000
#' burnin=5000
#' 
#' #Pseudo-Code 2 with updates for shape parameters
#' ind=sample(n,n)
#' holdout=data[ind[1:floor(0.1*n)]]
#' train=data[ind[(floor(0.1*n)+1):n]]
#' Xtrain=as.matrix(X[ind[(floor(0.1*n)+1):n],])
#' Gtrain=G[ind[(floor(0.1*n)+1):n],]
#' Xhold=as.matrix(X[ind[1:floor(0.1*n)],])
#' Ghold=G[ind[1:floor(0.1*n)],]
#' 
#' output<-GibbsBernoulliMLB(B,train,Xtrain,Gtrain,loess.smooth(points, data)$y,itmax = 100)
#' 
#' #check some trace plots
#' plot((as.mcmc(output$beta[1,burnin:B])))
#' plot(as.mcmc(output$alphab[1,burnin:B]))
#' plot(as.mcmc(output$alphae[1,burnin:B]))
#' plot(as.mcmc(output$alphax[1,burnin:B]))
#' plot(as.mcmc(output$kappab[1,burnin:B]))
#' plot(as.mcmc(output$kappae[1,burnin:B]))
#' plot(as.mcmc(output$kappax[1,burnin:B]))
#' plot((as.mcmc(output$eta[5,burnin:B])))
#' plot((as.mcmc(output$xi[10,burnin:B])))
#' 
#' nu.gibbsf=X%*%(output$beta)+G%*%(output$eta)
#' ximean=mean(output$xi[,burnin:B])
#' for (b in 1:B){
#'   nu.gibbsf[,b]=nu.gibbsf[,b]+ximean
#' }
#' pi.gibbsf= exp(nu.gibbsf)/(1+exp(nu.gibbsf))
#' 
#' #hows our estimate of the proportions?
#' pihatf = apply(pi.gibbsf[,burnin:B], 1, mean)
#' piboundsf = apply(pi.gibbsf[,burnin:B], 1, quantile, probs = c(0.025,0.975))
#' 
#' plot(points,pi1,ylim=c(0,1))
#' lines(sort(points),pihatf,col="red")
#' lines(sort(points),piboundsf[1,],col="blue")
#' lines(sort(points),piboundsf[2,],col="blue")
#' @export
GibbsBernoulliMLBtt<-function(B,data,X,G,ini,report=100,itmax=20,nslice=2){
  n = length(data)
  p = dim(X)[2]
  r = dim(G)[2]
  nn=matrix(1,n,1)
  
  #dataset needs to be at least 10
  ind=sample(n,n)
  train=data[ind[(floor(0.1*n)+1):(floor(n))]]
  valid=data[ind[1:floor(0.1*n)]]
  
  Hbeta = rbind(as.matrix(X[ind[(floor(0.1*n)+1):(floor(n))],]), as.matrix(X[ind[(floor(0.1*n)+1):(floor(n))],]),diag(p))
  HpHinvBeta = solve(t(Hbeta)%*%Hbeta)
  
  H = rbind(diag(length(train)), diag(length(train)),diag(length(train)))
  
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
    
    sv=1
    
    alphae = matrix(1,1,B)
    kappae = matrix(1.1,1,B)
    invV = diag(r)
    Gamma = diag(r)
    
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
        like= -p*(lgamma(x)) - p*(lgamma(kappab[b-1] -x))-x +x*sum(beta.gibbs[,b])
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
      Heta = rbind(G, G,invV)
      HpHinvEta = solve(t(Heta)%*%Heta)
      eta.gibbs[,b] = HpHinvEta%*%t(Heta)%*%W
      
      #update shapes eta using slice sampler
      f<-function(x){
        
        like=(-r*(lgamma(x)) - r*(lgamma(kappae[b-1] -x))-x +x*sum(invV%*%eta.gibbs[,b]))
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
        
        like =(r*(lgamma(x)) - r*(lgamma(x -alphae[b]))-x -x*sum(log(1+exp(invV%*%eta.gibbs[,b]))))
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
        like=(-n*(lgamma(x)) - n*(lgamma(kappax[b-1] -x))-x +x*sum(xi.gibbs[,b]))
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
      
      #update V
      etatemp=sv*eta.gibbs[,b]
      for (jj in 2:r){
        alph_post = c(alphae[b],1000*matrix(1,jj-1,1))
        kappa_post = c(kappae[b],2000*matrix(1,jj-1,1))
        musv = rbind(etatemp[jj],matrix(0,jj-1,1))
        Hv = t(cbind(etatemp[1:(jj-1)],diag(jj-1)))
        hh= logitbetasim(alph_post,kappa_post);
        hh = ginv(t(Hv)%*%Hv)%*%t(Hv)%*%musv+ginv(t(Hv)%*%Hv)%*%t(Hv)%*%hh
        Gamma[jj,1:(jj-1)]= hh
      }
      
      
      f<-function(x){
        if (x<0){
          like=-Inf
        }else{
        like =r*log(x)-x +x*alphae[b]*sum(Gamma%*%eta.gibbs[,b])-kappae[b]*sum(log(1+exp(x*Gamma%*%eta.gibbs[,b])))
        }
        if (x<0.01){
          like = -Inf;
        }
        
        if (x>20000){
          like = -Inf;
        }
        
        return(like)
        
      }
      sv2<-MfU.Sample.Run(sv, f, nsmp = nslice)
      sv =abs(sv2[nslice])
      
      
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
  
  
  SimpleClassifier<-function(ini,phat,holdout=NULL,phathold=NULL){
    
    cutoff<-function(lambda){
      estpos=phat>lambda
      helliginer= sum((sqrt(ini)-sqrt(estpos))^2)
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
              
              output<-GibbsBinomialMLB2(100,train,matrix(1,length(train),1),as.matrix(X[ind[(floor(0.1*n)+1):(floor(n))],]),G[ind[(floor(0.1*n)+1):(floor(n))],],report=1e15,rho=z[1],eps=z[2],z[3],z[4],z[5],2)
              
              nuest=X[ind[(floor(0.1*n)+1):(floor(n))],]%*%output$beta + G[ind[(floor(0.1*n)+1):(floor(n))],]%*%output$eta+mean(output$xi)
              nuesthold=X[ind[1:floor(0.1*n)],]%*%output$beta + G[ind[1:floor(0.1*n)],]%*%output$eta+mean(output$xi)
              
              piest = exp(nuest)/(1+exp(nuest))
              phat = rowMeans(piest[,50:100])
              
              piesthold = exp(nuesthold)/(1+exp(nuesthold))
              phathold = rowMeans(piesthold[,50:100])
              
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
  estpar<-optim(initB,fz,method = "Nelder-Mead",control=list(maxit=itmax))
  rho=estpar$par[1]
  eps=estpar$par[2]
  cbeta=estpar$par[3]
  ceta=estpar$par[4]
  cxi=estpar$par[5]
  
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
  
  sv=1
  
  simalphb=1
  simkappb = 1.1
  simalphe=1
  simkappe = 1.1
  simalphx=1
  simkappx = 1.1
  Hbeta = rbind(as.matrix(X), as.matrix(X),diag(p))
  HpHinvBeta = solve(t(Hbeta)%*%Hbeta)
  
  Heta = rbind(G, G,diag(r))
  HpHinvEta = solve(t(Heta)%*%Heta)
  H = rbind(diag(length(data)), diag(length(data)),diag(length(data)))
  invV = diag(r)
  Gamma = diag(r)
  
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
      like= -p*(lgamma(x)) - p*(lgamma(kappab[b-1] -x))-x +x*sum(beta.gibbs[,b])
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
    Heta = rbind(G, G,invV)
    HpHinvEta = solve(t(Heta)%*%Heta)
    eta.gibbs[,b] = HpHinvEta%*%t(Heta)%*%W
    
    #update shapes eta using slice sampler
    f<-function(x){
      
      like=(-r*(lgamma(x)) - r*(lgamma(kappae[b-1] -x))-x +x*sum(invV%*%eta.gibbs[,b]))
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
      
      like =(r*(lgamma(x)) - r*(lgamma(x -alphae[b]))-x -x*sum(log(1+exp(invV%*%eta.gibbs[,b]))))
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
    xi.gibbs[,b] = (1/3)*t(H)%*%W
    
    if (sum(is.nan(xi.gibbs[,b]))>0){
      xi.gibbs[,b] = xi.gibbs[,b-1]
    }
    
    if (sum(is.infinite(xi.gibbs[,b]))>0){
      xi.gibbs[,b] = xi.gibbs[,b-1]
    }
    
    #update shape xi using slice sampler
    f<-function(x){
      like=(-n*(lgamma(x)) - n*(lgamma(kappax[b-1] -x))-x +x*sum(xi.gibbs[,b]))
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
    
    etatemp=sv*eta.gibbs[,b]
    for (jj in 2:r){
      alph_post = c(alphae[b],1000*matrix(1,jj-1,1))
      kappa_post = c(kappae[b],2000*matrix(1,jj-1,1))
      musv = rbind(etatemp[jj],matrix(0,jj-1,1))
      Hv = t(cbind(etatemp[1:(jj-1)],diag(jj-1)))
      hh= logitbetasim(alph_post,kappa_post);
      hh = ginv(t(Hv)%*%Hv)%*%t(Hv)%*%musv+ginv(t(Hv)%*%Hv)%*%t(Hv)%*%hh
      Gamma[jj,1:(jj-1)]= hh
    }
    
    
    f<-function(x){
      if (x<0){
        like=-Inf
      }else{
        like =r*log(x)-x +x*alphae[b]*sum(Gamma%*%eta.gibbs[,b])-kappae[b]*sum(log(1+exp(x*Gamma%*%eta.gibbs[,b])))
      }
      if (x<0.01){
        like = -Inf;
      }
      
      if (x>20000){
        like = -Inf;
      }
      
      return(like)
      
    }
    sv2<-MfU.Sample.Run(sv, f, nsmp = nslice)
    sv =abs(sv2[nslice])
    
    invV = sv*Gamma;
    
    
    if (b>(report-1)){
      if ((max(b,report) %% min(b,report)) ==0){
        print(paste("MCMC Replicate: ",b))
      }
    }
    
  }
  
  nu.gibbs = X%*%(beta.gibbs)+G%*%(eta.gibbs) + xi.gibbs
  
  #reorganize mcmc reps of nu into pi according to Section 2.1
  pi.gibbs = exp(nu.gibbs)/(1+exp(nu.gibbs))
  
  
  nu.smooth.gibbs=X%*%beta.gibbs+G%*%eta.gibbs
  for (b in 1:B){
    nu.smooth.gibbs[,b]=nu.smooth.gibbs[,b]+mean(xi.gibbs[,b])
  }
  pi.smooth.gibbs= exp(nu.smooth.gibbs)/(1+exp(nu.smooth.gibbs))
  
  
  output<-list(pi=pi.gibbs,pi.smooth=pi.smooth.gibbs,beta=beta.gibbs,eta=eta.gibbs,xi=xi.gibbs,alphab=alphab,kappab=kappab,alphae=alphae,kappae=kappae,alphax=alphax,kappax=kappax,hypparams=estpar)
  return(output)
}
