#' Poisson LCM
#'
#' This code implements a version of the collapsed Gibbs sampler from Bradley et al. (2018). The model assumes that data follow a poisson distribution with a log link to a mixed effects model. The priors and random effects are assumed to follow a multivariate log-gamma distribution.
#' @param Niter The number of iterations of the collapsed Gibbs sampler
#' @param X An nxp matrix of covariates. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param G An nxr matrix of basis functions. Each column represents a basis function. Each row corresponds to one of the n replicates.
#' @param data An n dimensional vector consisting of count-valued observations.
#' @param sigbeta Prior variance on the regression parameters.
#' @param printevery Option to print the iteration number of the MCMC.
#' @param updatekappa Updating kappa is not needed if X contains an intercept. The default is FALSE.
#' @param jointupdate Block update beta and eta together. Default is TRUE.
#' @import actuar MfUSampler stats coda MASS
#' @return betas A pxNiter matrix of MCMC values for beta. The vector beta corresponds to the covariates X.
#' @return etas A rxNiter matrix of MCMC values for eta. The vector eta corresponds to the basis functions G.
#' @return deltas A nxNiter matrix of MCMC values for delta. The vector delta corresponds to uncorrelated random effects.
#' @return lambda_rep A nxNiter matrix of MCMC values for the mean of the Poisson dataset, modeled with, exp(Xbeta+Geta+delta).
#' @return alpha_b A 1xNiter matrix of MCMC values for the shape parameter associated with beta.
#' @return alpha_eta A 1xNiter matrix of MCMC values for the shape parameter associated with eta.
#' @return alpha_delta A 1xNiter matrix of MCMC values for the shape parameter associated with delta.
#' @examples
#' #load the necessary packages
#' library(CM)
#'
#' set.seed(123)
#' #define a test function
#' #A non-linear test function
#' lambda <- function(t) exp(1.1 + sin(2*pi*t))
#'
#' #define some 1-d locations
#' points = seq(0,1,length.out=1001)
#' points=points[2:1001]
#' m = dim(as.matrix(points))[1]
#'
#' #get the true mean at these locations
#' truemean<-matrix(0,m,1)
#' for (j in 1:length(points)){
#'   truemean[j,1] = lambda(points[j])
#' }
#'
#' #simulate the data
#' data = matrix(0,m,1)
#' for (i in 1:m){
#'   data[i] = rpois(1,truemean[i])
#' }
#'
#' #see how many zeros there are
#' sum(data==0)
#'
#' #plot the data
#' plot(data,xlab="time",ylab="Poisson counts",main="Counts vs. time")
#'
#' #covariate intercept-only
#' X = matrix(1,m,1)
#' p <- dim(X)[2]
#'
#'
#' ##compute the basis function matrix
#' #compute thin-plate splines
#' r = 8
#' knots = seq(0,1,length.out=r)
#'
#' #orthogonalize G
#' G = THINSp(as.matrix(points,m,1),as.matrix(knots,r,1))
#' outG<-qr(G)
#' G<-qr.Q(outG)
#'
#' #orthogonalize X
#' outX<-qr(X)
#' X<-qr.Q(outX)
#'
#'
#' #Run the MCMC algorithm
#' output<-GibbsPoissonMLG(Niter=2000,X,G,data)
#'
#' #trace plots (without burnin)
#' plot(as.mcmc(output$betas[1000:2000]))
#' plot(as.mcmc(output$etas[1,1000:2000]))
#' plot(as.mcmc(output$etas[8,1000:2000]))
#' plot(as.mcmc(output$deltas[10,1000:2000]))
#' plot(as.mcmc(output$alpha_eta[1,1000:2000]))
#' #alpha_delta does not mix well, one could thin, however, we choose not to. See MacEachern and Berliner (1994) for a discussion on thinning.
#' plot(as.mcmc(output$alpha_delta[1,1000:2000]))
#'
#' #estimates (remove a burnin)
#' lambda_est = apply(output$lambda_rep[,1000:2000],1,mean)
#' lambda_lower= apply(output$lambda_rep[,1000:2000],1,quantile,0.025)
#' lambda_upper= apply(output$lambda_rep[,1000:2000],1,quantile,0.975)
#'
#' #plot estimates and truth
#' plot(1:m,truemean,ylim = c(0,max(lambda_upper)+1))
#' lines(1:m,lambda_est,col="red")
#' lines(1:m,lambda_lower,col="blue")
#' lines(1:m,lambda_upper,col="blue")
#'
#' #smooth estimates (remove a burnin)
#' lambda_est = apply(output$lambda_rep_smooth[,1000:2000],1,mean)
#' lambda_lower= apply(output$lambda_rep_smooth[,1000:2000],1,quantile,0.025)
#' lambda_upper= apply(output$lambda_rep_smooth[,1000:2000],1,quantile,0.975)
#'
#' #plot smooth estimates and truth
#' plot(1:m,truemean,ylim = c(0,max(lambda_upper)+1))
#' lines(1:m,lambda_est,col="red")
#' lines(1:m,lambda_lower,col="blue")
#' lines(1:m,lambda_upper,col="blue")
#' covmat = matrix(0,1000,1000)
#' for (jj in 1:1000){
#'   covmat = covmat+(output$lambda_rep[,1000+jj] - lambda_est)%*%t((output$lambda_rep[,1000+jj] - lambda_est))/1001
#'   print(jj)
#' }
#' vars = 1/sqrt(diag(covmat))
#' corrmat= diag(vars)%*%covmat%*% diag(vars)
#' image(corrmat)
#'
#' @export
GibbsPoissonMLG<-function(Niter=2000,X,G,data, sigbeta=10,printevery=100,updatekappa=FALSE,jointupdate=TRUE){

  p=dim(X)[2]
  if (jointupdate==TRUE){
    G = cbind(X,G)
  }

  m = length(data)
  r = dim(G)[2]
  sigbeta=sigbeta*diag(p)


  ### Gibbs sampler
  #initializations
  kappa_k=matrix(1,m,Niter)
  betas=matrix(0,p,Niter)
  alpha_b=matrix(1,1,Niter)
  kappa_b=matrix(1,1,Niter)
  etas=matrix(1,r,Niter)
  alpha_eta=matrix(1,1,Niter)
  kappa_eta=matrix(1,1,Niter)
  deltas=matrix(1,m,Niter)
  alpha_delta=matrix(1,1,Niter)
  kappa_delta=matrix(1,1,Niter)
  kappak = matrix(1,Niter,1)
  kappak2 = matrix(1,Niter,1)
  kappavi=10000
  mgam2e=1

  (m22 <- matrix(0.5, r, r))
  m22[upper.tri(m22)]=0
  diag(m22) = 1

  variances1 = m22
  variances2 = diag(m)
  halceta=1
  halcdelta=1
  halck=1
  halcbeta=1
  alphvi=10000

  sigdelta=1*matrix(1,Niter,1)
  mgam2d=1

  sum(data==0)
  zetaadj = data
  if (sum(data==0)>0){zetaadj =data +0.5#this value (i.e., 0.5) may not be best
  }

  mgam2b=1

  mgam2v=1
  varinclude=1
  variances1=diag(r)

  ve=1e-5


  #function to pass through slice sampler update of shape parameters
  f1<-function(a,b,betag,w1,w2,rho,bound,bound2){

    si = (length(betag));

    like = (w1+sum(betag))*a - (w2+sum(exp(betag)))*b - (rho+si)*(lgamma(a)) +(rho+si)*a*log(b);

    if (a<bound){
      like = -Inf;
    }

    if (a>bound2){
      like = -Inf;
    }

    if (b<0){
      like = -Inf;
    }

    return(like)

  }




  boundb=0.5
  bounde=0.5
  boundd=0.5

  print("Collapsed Gibbs Sampler Starting")
  for (t in 2:Niter)
  {

    if (jointupdate==FALSE){
    #update alpha_b
    Hbetas=rbind(X,sigbeta)
    fb<-function(x){
      return(f1(x,kappa_b[t-1],Hbetas*betas[,t-1],1,-1e-15,1,boundb,20000))}
    alpha_b_temp<-MfU.Sample.Run(alpha_b[t-1], fb, nsmp = 11)
    alpha_b[t] = abs(alpha_b_temp[11])
    if (updatekappa){
     sc = 1/(sum(exp(Hbetas*betas[,t-1]))+mgam2b);
      kappa_b[t] = rgamma(1,p*alpha_b[t]+alpha_b[t]+1,scale=sc);
     mgam2b = rgamma(1,1,scale=1/(kappa_b[t]+1));
    }
    }

    #update alpha_eta
    Hetas=rbind(G,variances1)
    fe<-function(x){
      return(f1(x,kappa_eta[t-1],Hetas%*%etas[,t-1],1,-1e-15,1,bounde,20000))}
    alpha_e_temp<-MfU.Sample.Run(alpha_eta[t-1], fe, nsmp = 11)
    alpha_eta[t] = abs(alpha_e_temp[11])
    if (updatekappa){
      sc = 1/(sum(exp(Hetas%*%etas[,t-1]))+mgam2e);
      kappa_eta[t] = rgamma(1,r*alpha_eta[t]+alpha_eta[t]+1,scale=sc);
      mgam2e = rgamma(1,1,scale=1/(kappa_eta[t]+1));
    }

    fd<-function(x){
      return(f1(abs(x),kappa_delta[t-1],sigdelta[t-1]*deltas[,t-1],1,-1e-15,1,boundd,20000))}
    alpha_d_temp<-MfU.Sample.Run(alpha_delta[t-1], fd, nsmp = 11)
    alpha_delta[t] = abs(alpha_d_temp[11])
    if (updatekappa){
      sc = 1/(sum(exp(sigdelta[t-1]*deltas[,t-1]))+mgam2d);
      kappa_delta[t] = rgamma(1,m*alpha_delta[t]+alpha_delta[t]+1,scale=sc);
      mgam2d = rgamma(1,1,scale=1/(kappa_delta[t]+1));
    }


    if (jointupdate==FALSE){
    ###update betas #####
    Hbetas=rbind(X,sigbeta)
    alphabetas=rbind(zetaadj+alpha_b[t],alpha_b[t]*matrix(1,p,1))
    kappabetas=rbind(exp(G%*%etas[,t-1]+deltas[,t-1])+kappa_b[t],kappa_b[t]*matrix(1,p,1))
    betas[,t]=solve(t(Hbetas)%*%Hbetas)%*%t(Hbetas)%*%log(rgamma((m+p),alphabetas,rate=kappabetas))
    }

    ###update etas #####
    if (jointupdate==FALSE){
      Hetas=rbind(G,variances1)
      alphaetas=rbind(zetaadj+alpha_eta[t],alpha_eta[t]*matrix(1,r,1))
      kappaetas=rbind(exp(X%*%betas[,t]+deltas[,t-1])+kappa_eta[t],kappa_eta[t]*matrix(1,r,1))
      etas[,t]=solve(t(Hetas)%*%Hetas)%*%t(Hetas)%*%log(rgamma((m+r),alphaetas,rate=kappaetas))
    }
    if (jointupdate==TRUE){
    Hetas=rbind(G,variances1)
    alphaetas=rbind(zetaadj+alpha_eta[t],alpha_eta[t]*matrix(1,r,1))
    kappaetas=rbind(exp(matrix(deltas[,t-1], ncol=1))+kappa_eta[t],kappa_eta[t]*matrix(1,r,1))
      #rbind(exp(X%*%betas[,t]+deltas[,t-1])+kappa_eta[t],kappa_eta[t]*matrix(1,r,1))
    etas[,t]=solve(t(Hetas)%*%Hetas)%*%t(Hetas)%*%log(rgamma((m+r),alphaetas,rate=kappaetas))
    }
    #### update deltas ####
    if (jointupdate==FALSE){
      alphadeltas=zetaadj+alpha_delta[t]*matrix(1,m,1)
      kappadeltas=exp(X%*%betas[,t]+G%*%etas[,t])+kappa_delta[t]*matrix(1,m,1)
      deltas[,t]=log(rgamma(m,alphadeltas,rate=kappadeltas))
    }
    if (jointupdate==TRUE){
      alphadeltas=zetaadj+alpha_delta[t]*matrix(1,m,1)
      kappadeltas=exp(G%*%etas[,t])+kappa_delta[t]*matrix(1,m,1)
      deltas[,t]=log(rgamma(m,alphadeltas,rate=kappadeltas))
    }

    #update var of eta
    fve<-function(x){
      variances2=variances1
      diag(variances2)=x
      Hetas=rbind(G,variances2)
      vv=f1(alpha_eta[t],kappa_eta[t],Hetas%*%etas[,t-1],1,-1e-15,1,bounde,20000)
      return(vv)
    }
    ve_temp<-MfU.Sample.Run(ve, fve, nsmp = 11)
    ve = abs(ve_temp[11])


    ###########update IV###########
    variances1b = matrix(0,1,r)
    variances1b[1] =ve
    for(k in 2:r){      ##### since diags has diag(V), k should start from 2 to Z+1
      Vk=rep(0,k-1)

      H_iv=rbind(diag(etas[1:(k-1),t],k-1,k-1),1000*diag(k-1))
      shape_iv=c(alpha_eta[t]*matrix(1,length(etas[1:(k-1),t]),1),alphvi*matrix(1,k-1,1))
      K_iv=c(kappa_eta[t-1]*matrix(1,length(etas[1:(k-1),t]),1),kappavi*matrix(1,k-1,1))
      mu_iv = c(-etas[k,t],matrix(0,length(shape_iv)-1,1))
      Vk=ginv(t(H_iv)%*%H_iv)%*%t(H_iv)%*%mu_iv + ginv(t(H_iv)%*%H_iv)%*%t(H_iv)%*%(log(rgamma(length(shape_iv),shape_iv,rate=1))-log(K_iv))

      Vtemp = matrix(0,1,r)
      Vtemp[1:length(Vk)] = Vk
      Vtemp[k] = ve
      variances1b = rbind(variances1b,Vtemp)
    }
    variances1=variances1b




    #update alpha_vi
    variances11=variances1
    variances11[variances11==ve]=0
    variances11=variances11[abs(variances11)>0]
    fv<-function(x){
      return(f1(x,kappavi,variances11,1,-1e-15,1,0.5,20000))}
    alpha_d_temp<-MfU.Sample.Run(alphvi, fv, nsmp = 11)
    alphvi = abs(alpha_d_temp[11])
    sc = 1/(sum(exp(variances11))+mgam2v);
    kappavi = rgamma(1,length(variances11)*alphvi+alphvi+1,scale=sc);


    if ((t%%printevery)==0){
    print(paste("MCMC Replicate: ",t))
    }

  }

  lambda_rep = matrix(0,dim(G)[1],Niter)
  for (j in 1:Niter){
    lambda_rep[,j] = exp(X%*%betas[,j]+G%*%etas[,j]+deltas[,j])
  }

  lambda_rep_smooth = matrix(0,dim(G)[1],Niter)
  deltmean = mean(deltas[,(Niter/2):Niter])
  for (j in 1:Niter){
    lambda_rep_smooth[,j] = exp(X%*%betas[,j]+G%*%etas[,j]+deltmean)
  }

  if (jointupdate==FALSE){

  output = list(betas,etas,deltas,lambda_rep,lambda_rep_smooth,alpha_b,alpha_eta,alpha_delta)

  names(output) <- c("betas","etas","deltas","lambda_rep","lambda_rep_smooth","alpha_b","alpha_eta","alpha_delta")

  return(output)}

  if (jointupdate==TRUE){

    betas = etas[1:p,]
    etas = etas[(p+1):r,]

    output = list(betas,etas,deltas,lambda_rep,lambda_rep_smooth,alpha_eta,alpha_delta)

    names(output) <- c("betas","etas","deltas","lambda_rep","lambda_rep_smooth","alpha_eta","alpha_delta")

    return(output)}

}
