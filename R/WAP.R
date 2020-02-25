#' WAP
#'
#' This code implements a version of the collapsed Gibbs sampler from Xu et al. (2019). The model assumes that data follow Weibull and Poisson distributions with a log link to a mixed effects model. The priors and random effects are assumed to follow a multivariate log-gamma distribution.
#' @param iter The number of iterations of the collapsed Gibbs sampler
#' @param X1 An nxp matrix of covariates for the Weibull response. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param psi1 An nxr matrix of basis functions for the Weibull response. Each column represents a basis function. Each row corresponds to one of the n replicates.
#' @param Y1 An n dimensional vector consisting of totals associated with Weibull observations.
#' @param X2 An nxp matrix of covariates for the Poisson response. Each column represents a covariate. Each row corresponds to one of the n replicates.
#' @param psi2 An nxr matrix of basis functions for the Poisson response. Each column represents a basis function. Each row corresponds to one of the n replicates.
#' @param Y2 An n dimensional vector consisting of totals associated with Poisson observations.
#' @param num_A Hyperparameter for the adaptive rejection algorithm
#' @param burn_in The size of the burn-in for the MCMC.
#' @param printevery Option to print the iteration number of the MCMC.
#' @param output_type Can be specified to be the entire chain 'chain', or the posterior means 'mean'.
#' @param alpha_rho The shape parameter of the prior for the Weibull distribution's shape parameter.
#' @param kappa_rho The rate parameter of the prior for the Weibull distribution's shape parameter.
#' @param alpha_rho The shape parameter of the prior for the elements of the covariance parameter.
#' @param alpha_rho The rate parameter of the prior for the elements of the covariance parameter.
#' @param prior_initials The initial values of the shape parameters
#' @param nslice The burnin for the slice sampler
#' @import actuar MfUSampler stats coda MASS Matrix
#' @return ans A list of updated parameter values from WAP (when output_type = 'chain'). A list of posterior means from WAP (when output_type='mean').
#' @examples
#'
#' library(CM)
#'
#' set.seed(123000)
#'n = 100
#'locs = seq(-2*pi,2*pi,length.out=n)
#'
#'a1=-0.5
#'b1=-1
#'a2=1
#'b2=0.5
#'
#'q1 = a1+b1*sin(locs)
#'q2 = a2 + b2*q1
#'
#'Y1 = q1+sqrt(var(a1+b1*sin(locs))/(n*5))*rnorm(n)
#'Y2 = q2+sqrt(var(a1 + b2*q1)/(n*5))*rnorm(n)
#'
#'X1=matrix(1,n,1)
#'X2 = matrix(1,n,1)
#'
#'r = 10
#'knots = seq(-2*pi,2*pi,length.out=r)
#'psi1 = THINSp(as.matrix(locs),as.matrix(knots))
#'psi2 = psi1
#'
#'Z1=rweibull(n,1,exp(-q1))
#'Z2=rpois(n,exp(q2))
#'
#'#plot data
#'plot(Z2)
#'plot(Z1)
#'
#'# example here is for the chain output
#'ans=WAP(X1,X2,Z1,Z2,psi1,psi2,num_A=100,iter=2000,burn_in = 1000,nslice=50)
#'iter = 400
#'burn_in = 100
#'
#'#estimate versus truth
#'q1_mcmc=X1%*%ans$beta1+psi1%*%ans$eta+psi1%*%ans$eta1+ans$gamma1;
#'q1hat = apply(q1_mcmc, 1, mean)
#'q1bounds = apply(q1_mcmc, 1, quantile, probs = c(0.025,0.975))
#'
#'plot(locs,q1,ylim=range(c(q1bounds)))
#'lines(sort(locs),q1hat,col="red")
#'lines(sort(locs),q1bounds[1,],col="blue")
#'lines(sort(locs),q1bounds[2,],col="blue")
#'
#'q2_mcmc=X2%*%ans$beta2+psi2%*%ans$eta+psi2%*%ans$eta2+ans$gamma2;
#'q2hat = apply(q2_mcmc, 1, mean)
#'q2bounds = apply(q2_mcmc, 1, quantile, probs = c(0.025,0.975))
#'
#'plot(locs,q2,ylim=range(c(q1bounds)))
#'lines(sort(locs),q2hat,col="red")
#'lines(sort(locs),q2bounds[1,],col="blue")
#'lines(sort(locs),q2bounds[2,],col="blue")
#' @export
WAP<-function(X1=X,X2=X,Y1=Y1,Y2=Y2,psi1,psi2,num_A=10,
                   iter=5000, burn_in=3000,output_type='chain',
                   alpha_rho=0.001,kappa_rho=1000,alpha_iv=1000,kappa_iv=0.001,
                   prior_initials=rep(1,17),nslice=2,printevery=100){
  # choose the output type of the MCMC, 'chain' for the whole MCMC chain after burnin
  # 'mean' for mean estimate of MCMC.
  # N is the number of observations.
  N1 <- length(Y1)
  N2 <- length(Y2)
  rank1<-dim(psi1)[2]
  rank2<-dim(psi2)[2]
  r=max(rank1,rank2)
  p <-dim(X2)[2]
  p1<-dim(X1)[2];p2<-dim(X2)[2]
  l_eta=1
  logl_eta1=1
  logl_eta2=1
  logl_beta1=1
  logl_beta2=1
  logl_gamma1=1
  logl_gamma2=1

  if(output_type=='chain'){
    # chain of hyperparameters in beta
    alphas_beta1=rep(0,iter);kappas_beta1=rep(0,iter);
    alphas_beta2=rep(0,iter);kappas_beta2=rep(0,iter);
    # eta chain
    etas=matrix(0,r,iter)
    alphas_eta=rep(0,iter);
    kappas_eta=rep(0,iter);
    # eta1 chain
    eta1s=matrix(0,rank1,iter)
    kappas_eta1=rep(0,iter);
    alphas_eta1=rep(0,iter);
    # eta2 chain
    eta2s=matrix(0,rank2,iter)
    kappas_eta2=rep(0,iter);
    alphas_eta2=rep(0,iter);
    # beta chain
    beta1s=matrix(0,p1,iter);
    beta2s=matrix(0,p2,iter);
    # Weibull shape chain rho
    rhos=matrix(0,N1,iter);
    tt_rhos=matrix(0,num_A,iter);
    # Inverse V chain
    IVs=list();
    # Inverse V1 chain
    IV1s=list();
    # Inverse V2 chain
    IV2s=list();
    # gamma chain
    gamma1s = matrix(0,N1,iter);
    gamma2s = matrix(0,N2,iter);
    # gamma1 prior chain
    alphas_gamma1=rep(0,iter);
    kappas_gamma1=rep(0,iter);
    # gamma2 prior chain
    alphas_gamma2=rep(0,iter);
    kappas_gamma2=rep(0,iter);
  }else{
    if(output_type=='mean'){
      # parameters for distributions
      mean_etas=rep(0,r);
      mean_eta1s=rep(0,rank1);
      mean_eta2s=rep(0,rank2);
      mean_beta1s=rep(0,p1);
      mean_beta2s=rep(0,p2);
      mean_gamma1s=rep(0,N1);
      mean_gamma2s=rep(0,N2);
      mean_rhos=matrix(0,N1,1)
    }
  }
  # inital of parameters
  eta=rep(0,r)
  eta1=rep(0,rank1)
  eta2=rep(0,rank2)
  beta1=rep(0,p1)
  beta2=rep(0,p2)
  gamma1=rep(0,N1)
  gamma2=rep(0,N2)
  # inital of hyperparameters
  alpha_beta1=prior_initials[1];
  kappa_beta1=prior_initials[2];
  alpha_beta2=prior_initials[3];
  kappa_beta2=prior_initials[4];
  kappa_eta=prior_initials[5]
  alpha_eta=prior_initials[6]
  kappa_eta1=prior_initials[7]
  alpha_eta1=prior_initials[8]
  kappa_eta2=prior_initials[9]
  alpha_eta2=prior_initials[10]
  beta1=rep(prior_initials[11],p1)
  beta2=rep(prior_initials[12],p2)
  tt_rho=rep(prior_initials[13],num_A)
  rho=matrix(0,N1,1)
  rho=rep(tt_rho,each=N1/num_A)
  IV=matrix(runif(r*r,0,1),r,r)
  IV[upper.tri(IV)]=0;diag(IV)=1
  IV1=matrix(runif(rank1*rank1,0,1),rank1,rank1)
  IV1[upper.tri(IV1)]=0;diag(IV1)=1
  IV2=matrix(runif(rank2*rank2,0,1),rank2,rank2);
  IV2[upper.tri(IV2)]=0;diag(IV2)=1
  alpha_gamma1=prior_initials[14];
  kappa_gamma1=prior_initials[15];
  alpha_gamma2=prior_initials[16];
  kappa_gamma2=prior_initials[17];

  # likelihood of M-H's proposal function
  lh<-function(rho,shape.w,rho.t){
    sum(log(rho)-log(shape.w)+log(rho.t)*(rho-1)+(-rho.t^rho/shape.w))+(alpha_rho-1)*log(rho)-rho/kappa_rho
  }


  # initial settings
   ## the variance is very small
  zeta_pois=0.5;zeta_weib=0 ## parameter for vi

  indic_weib= Y1>0
  indic_pois= Y2>0
  zfudge_weib = rep(zeta_weib,N1)
  zfudge_pois= rep(zeta_pois,N2)
  zfudge_weib[indic_weib==TRUE] =0
  zfudge_pois[indic_pois==TRUE] =0
  pp=1;r1=1000;r2=-1e-15
  if(rank1>=rank2){psi22=cbind(psi2,matrix(0,max(N1,N2),(rank1-rank2)));psi11=psi1}else
  {psi11=cbind(psi1,matrix(0,max(N1,N2),(rank2-rank1)));psi22=psi2}

  H_gamma1=sparseMatrix(i=c(1:(2*N1)),j=rep(1:N1,2),x=1)
  HH_gamma1=solve(t(H_gamma1)%*%H_gamma1)%*%t(H_gamma1)
  H_gamma2=sparseMatrix(i=c(1:(2*N2)),j=rep(1:N2,2),x=1)
  HH_gamma2=solve(t(H_gamma2)%*%H_gamma2)%*%t(H_gamma2)


  for (i in 2:iter)
  {
    H_eta=rbind(psi11,psi22,IV)
    K_eta_1=-rho*log(Y1)-(X1%*%beta1+psi1%*%eta1+gamma1)
    K_eta_2=-(X2%*%beta2+psi2%*%eta2+gamma2)
    K_eta_3=rep(-log(kappa_eta),r)
    K_eta=c(K_eta_1,K_eta_2,K_eta_3)
    shape_eta=c(rep(1,N1),Y2+zfudge_pois,rep(alpha_eta,r))
    m_eta=dim(H_eta)[1]
    b_eta=matrix(K_eta+log(rgamma(m_eta,shape=shape_eta,scale = 1)),nrow=m_eta,ncol=1)
    eta=solve(t(H_eta)%*%H_eta,tol=1e-50)%*%t(H_eta)%*%b_eta

    ############# beta 1 ################
    H_beta1=rbind(X1,diag(rep(1,p1)))
    K_beta1=c(-rho*log(Y1+zfudge_weib)-(psi11%*%eta+psi1%*%eta1+gamma1),rep(-log(kappa_beta1),p1))
    shape_beta1=c(rep(1,N1),rep(alpha_beta1,p1))
    m_beta1=dim(H_beta1)[1]
    b_beta1=matrix(K_beta1+log(rgamma(m_beta1,shape=shape_beta1,scale = 1)),nrow=m_beta1,ncol=1)
    beta1=as.numeric(ginv(t(H_beta1)%*%H_beta1,tol=1e-50)%*%t(H_beta1)%*%b_beta1)

    ###############gamma1 ################
    K_gamma1=c(-rho*log(Y1+zfudge_weib)-(psi11%*%eta+psi1%*%eta1+X1%*%beta1),rep(-log(kappa_gamma1),N1))
    shape_gamma1=c(rep(1,N1),rep(alpha_gamma1,N1))
    m_gamma1=dim(H_gamma1)[1]
    b_gamma1=matrix(K_gamma1+log(rgamma(m_gamma1,shape=shape_gamma1,scale = 1)),nrow=m_gamma1,ncol=1)
    gamma1=as.matrix(HH_gamma1%*%b_gamma1)

    ###########update eta1 #########
    H_eta1=rbind(psi1,IV1)
    K_eta1_1=-rho*log(Y1+zfudge_weib)-(X1%*%beta1+psi11%*%eta+gamma1)
    K_eta1_2=rep(-log(kappa_eta1),rank1)
    K_eta1=c(K_eta1_1,K_eta1_2)
    shape_eta1=c(rep(1,N1),rep(alpha_eta1,rank1))
    m_eta1=dim(H_eta1)[1]
    b_eta1=matrix(K_eta1+log(rgamma(m_eta1,shape=shape_eta1,scale = 1)),nrow=m_eta1,ncol=1)
    eta1=ginv(t(H_eta1)%*%H_eta1,tol=1e-50)%*%t(H_eta1)%*%b_eta1

    ###########update beta 2##########
    H_beta2=rbind(X2,diag(rep(1,p2)))
    K_beta2=c(-(psi22%*%eta+psi2%*%eta2+gamma2),rep(-log(kappa_beta2),p2))
    shape_beta2=c(Y2+zfudge_pois,rep(alpha_beta2,p2))
    m_beta2=dim(H_beta2)[1]
    b_beta2=matrix(K_beta2+log(rgamma(m_beta2,shape=shape_beta2,scale = 1)),nrow=m_beta2,ncol=1)
    beta2=as.numeric(ginv(t(H_beta2)%*%H_beta2,tol=1e-50)%*%t(H_beta2)%*%b_beta2)

    ###########update gamma 2##########
    K_gamma2=c(-(psi22%*%eta+psi2%*%eta2+X2%*%beta2),rep(-log(kappa_gamma2),N2))
    shape_gamma2=c(Y2+zfudge_pois,rep(alpha_gamma2,N2))
    m_gamma2=dim(H_gamma2)[1]
    b_gamma2=matrix(K_gamma2+log(rgamma(m_gamma2,shape=shape_gamma2,scale = 1)),nrow=m_gamma2,ncol=1)
    gamma2=as.numeric(HH_gamma2%*%b_gamma2)

    ###########update eta2 #########
    H_eta2=rbind(psi2,IV2)
    K_eta2_1=-(X2%*%beta2+psi22%*%eta+gamma2)
    K_eta2_2=rep(-log(kappa_eta2),rank2)
    K_eta2=c(K_eta2_1,K_eta2_2)
    shape_eta2=c(Y2+zfudge_pois,rep(alpha_eta2,rank2))
    m_eta2=dim(H_eta2)[1]
    b_eta2=matrix(K_eta2+log(rgamma(m_eta2,shape=shape_eta2,scale = 1)),nrow=m_eta2,ncol=1)
    eta2=ginv(t(H_eta2)%*%H_eta2,tol=1e-50)%*%t(H_eta2)%*%b_eta2

    ############### IV ###############
    New_IV=matrix(0,nrow = r,ncol=r)
    ##### diags included the diag of IV, so it has length of N
    for(j in 1:r){
      H_iv_j=rbind(eta[j],1)
      shape_iv=c(alpha_eta,alpha_iv)
      for(ii in j:r){
        Ik_star=log(kappa_eta)-(IV[ii,-c(j)]%*%eta[-c(j)])
        K_iv=c(Ik_star,-log(kappa_iv))
        m_iv=dim(H_iv_j)[1]
        b_iv=matrix(K_iv+log(rgamma(m_iv,shape = shape_iv,scale =1)),nrow=m_iv,ncol=1)
        New_IV[ii,j]=ginv(t(H_iv_j)%*%H_iv_j,tol=1e-50)%*%t(H_iv_j)%*%b_iv
      }
    }
      IV=New_IV


    ############### IV 1#############
    New_IV=matrix(0,nrow = rank1,ncol=rank1)
    ##### diags included the diag of V, so it has length of N
    for(j in 1:rank1){
      H_iv1_j=rbind(eta1[j],1)
      shape_iv1=c(alpha_eta1,alpha_iv)
      for(ii in j:rank1){
        Ik1_star=log(kappa_eta1)-(IV1[ii,-c(j)]%*%eta1[-c(j)])
        K_iv1=c(Ik1_star,-log(kappa_iv))
        m_iv1=dim(H_iv1_j)[1]
        b_iv1=matrix(K_iv1+log(rgamma(m_iv1,shape = shape_iv1,scale = 1)),nrow=m_iv1,ncol=1)
        New_IV[ii,j]=ginv(t(H_iv1_j)%*%H_iv1_j,tol=1e-50)%*%t(H_iv1_j)%*%b_iv1
      }
    }
      IV1=New_IV

    ############### IV 2#############
    New_IV=matrix(0,nrow = rank2,ncol=rank2)
    ##### diags included the diag of V, so it has length of N
    for(j in 1:rank2){
      H_iv2_j=rbind(eta2[j],1)
      shape_iv2=c(alpha_eta2,alpha_iv)
      for(ii in j:rank2){
        Ik2_star=log(kappa_eta2)-(IV2[ii,-c(j)]%*%eta2[-c(j)])
        K_iv2=c(Ik2_star,-log(kappa_iv))
        m_iv2=dim(H_iv2_j)[1]
        b_iv2=matrix(K_iv2+log(rgamma(m_iv2,shape = shape_iv2,scale = 1)),nrow=m_iv2,ncol=1)
        New_IV[ii,j]=ginv(t(H_iv2_j)%*%H_iv2_j,tol=1e-50)%*%t(H_iv2_j)%*%b_iv2
      }
    }
      IV2=New_IV

    ############## update rho #####
    ### use adaptive rejection algorithem
    q1it = X1%*%beta1+psi11%*%eta+psi1%*%eta1+gamma1
    weibull_scale= exp(-q1it)## use adaptive rejection algorithem

    for(j in 1:num_A){
      a=(j-1)*N1/num_A+1;b=j*N1/num_A
      new_rho=rgamma(1,shape=tt_rho[j],scale=1);if(new_rho==0){new_rho=0.01}
      P_new=lh(rho = new_rho,shape.w= weibull_scale[a:b],rho.t=Y1[a:b])
      P=lh(rho = tt_rho[j],shape.w= weibull_scale[a:b],rho.t=Y1[a:b])
      g1_new=(dgamma(new_rho,shape=tt_rho[j],scale=1,log=T))
      g1=(dgamma(tt_rho[j],shape =new_rho,scale=1,log=T))
      accept_temp=min(log(1),P_new-P+g1-g1_new)
      u=log(runif(1,0,1))
      if (is.na(accept_temp)==0){
        accept=accept_temp
      if(u<accept){tt_rho[j]=new_rho}
      }
    }
    rho=rep(tt_rho,each=N1/num_A)

    ######################priors for eta #######################
    eta_shap=alpha_eta*(pp+r)+1
    eta_r2=r2-sum(exp(ginv(IV)%*%eta))
    kappa_temp=rgamma(1,shape = eta_shap,rate=-eta_r2 )#+0.0001
    if (kappa_temp>0.01){
    kappa_eta =kappa_temp
    }

    f<-function(x){
      like =x* (p+r)*log(kappa_eta)-x-(p+r)*lgamma(x)+x*sum(ginv(IV)%*%eta)
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

    }

    l_eta_temp<-MfU.Sample.Run(alpha_eta, f, nsmp = nslice)
    alpha_eta = abs(l_eta_temp[nslice])

    ######################priors for eta1 #######################
    eta1_shap=(alpha_eta1*(pp+rank1)+1)
    eta1_r2=r2-sum(exp(ginv(IV1)%*%eta1))
    kappa_temp=rgamma(1,shape = eta1_shap,rate= -eta1_r2)#+0.0001
    if (kappa_temp>0.01){
      kappa_eta1 =kappa_temp
    }

    f<-function(x){
      like = x*(pp+r1)*log(kappa_eta1)-x-(pp+r1)*lgamma(x)+x*sum(ginv(IV1)%*%eta1)
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }



      return(like)

    }

    logl_eta1_temp<-MfU.Sample.Run(alpha_eta1, f, nsmp = nslice)
    alpha_eta1 = abs(logl_eta1_temp[nslice])

    ######################priors for eta2 #######################
    eta2_shap=(alpha_eta2*(pp+rank2)+1)
    eta2_r2=r2-sum(exp(ginv(IV2)%*%eta2))
    kappa_temp=rgamma(1,shape =eta2_shap,rate  = -eta2_r2)#+0.0001
    if (kappa_temp>0.01){
      kappa_eta2 =kappa_temp
    }

    f<-function(x){
      like = (x*(pp+r1)*log(kappa_eta2)-x-(pp+r1)*lgamma(x)+x*sum(ginv(IV2)%*%eta2))
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

    }

    logl_eta2_temp<-MfU.Sample.Run(alpha_eta2, f, nsmp = nslice)
    alpha_eta2 = abs(logl_eta2_temp[nslice])


    ##################### priors for beta1 ###################
    beta1_shap=(alpha_beta1*(pp+p1)+1)
    beta1_r2=r2-sum(exp(beta1))
    kappa_temp=rgamma(1,shape =beta1_shap,rate  =-beta1_r2)#+0.0001
    if (kappa_temp>0.01){
      kappa_beta1 =kappa_temp
    }

    f<-function(x){
      like = (x*(pp+r1)*log(kappa_beta1)-x-(pp+r1)*lgamma(x)+x*sum(beta1))
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

    }

    logl_beta1_temp<-MfU.Sample.Run(alpha_beta1, f, nsmp = nslice)
    alpha_beta1 = abs(logl_beta1_temp[nslice])

    ##################### priors for beta2 ###################
    beta2_shap=(alpha_beta2*(pp+p2)+1)
    beta2_r2=r2-sum(exp(beta2))
    kappa_temp=rgamma(1,shape = beta2_shap,rate  = -beta2_r2)#+0.0001
    if (kappa_temp>0.01){
      kappa_beta2 =kappa_temp
    }

    f<-function(x){
      like = x*(pp+r1)*log(kappa_beta2)-x-(pp+r1)*lgamma(x)+x*sum(beta2)
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

    }

    logl_beta2_temp<-MfU.Sample.Run(alpha_beta2, f, nsmp = nslice)
    alpha_beta2 = abs(logl_beta2_temp[nslice])

    ##################### priors for gamma1 ###################
    gamma1_shap=(alpha_gamma1*(pp+N1)+1)
    gamma1_r2=r2-sum(exp(gamma1))
    kappa_temp=rgamma(1,shape = gamma1_shap,rate  = -gamma1_r2)#+0.0001
    if (kappa_temp>0.01){
      kappa_gamma1 =kappa_temp
    }


    f<-function(x){
      like = x*(pp+r1)*log(kappa_gamma1)-x-(pp+r1)*lgamma(x)+x*sum(gamma1)
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

    }

    logl_gamma1_temp<-MfU.Sample.Run(alpha_gamma1, f, nsmp = nslice)
    alpha_gamma1 = abs(logl_gamma1_temp[nslice])

    ##################### priors for gamma2 ###################
    gamma2_shap=(alpha_gamma2*(pp+N2)+1)
    gamma2_r2=r2-sum(exp(gamma2))
    kappa_temp=rgamma(1,shape =gamma2_shap ,rate  = -gamma2_r2)#+0.0001
    if (kappa_temp>0.01){
      kappa_gamma2 =kappa_temp
    }


    f<-function(x){
      like = x*(pp+r1)*log(kappa_gamma2)-x-(pp+r1)*lgamma(x)+x*sum(gamma2)
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

    }

    logl_gamma2_temp<-MfU.Sample.Run(alpha_gamma2, f, nsmp = nslice)
    alpha_gamma2 = abs(logl_gamma2_temp[nslice])
    
    if ((i%%printevery)==0){
      print(paste("MCMC Replicate: ",i))
    }

    # output type
    if(output_type=='chain'){
      # parameters
      etas[,i]=eta
      eta1s[,i]=eta1
      eta2s[,i]=eta2
      beta1s[,i]=beta1
      beta2s[,i]=beta2
      gamma1s[,i]=gamma1
      gamma2s[,i]=gamma2
      tt_rhos[,i]=tt_rho
      rhos[,i]=rho
      # hyperparameters
      IVs[[i]]=IV
      IV1s[[i]]=IV1
      IV2s[[i]]=IV2
      alphas_beta1[i]=alpha_beta1
      kappas_beta1[i]=kappa_beta1
      alphas_beta2[i]=alpha_beta2
      kappas_beta2[i]=kappa_beta2
      alphas_eta[i]=alpha_eta
      kappas_eta[i]=kappa_eta
      alphas_eta1[i]=alpha_eta1
      kappas_eta1[i]=kappa_eta1
      alphas_eta2[i]=alpha_eta2
      kappas_eta2[i]=kappa_eta2
      alphas_gamma1[i]=alpha_gamma1
      kappas_gamma1[i]=kappa_gamma1
      alphas_gamma2[i]=alpha_gamma2
      kappas_gamma2[i]=kappa_gamma2
    }else{
      if(output_type=='mean'){
        if(i>burn_in){
          mean_etas=mean_etas+eta
          mean_eta1s=mean_eta1s+eta1
          mean_eta2s=mean_eta2s+eta2
          mean_beta1s=mean_beta1s+beta1;
          mean_beta2s=mean_beta2s+beta2;
          mean_gamma1s=mean_gamma1s+gamma1
          mean_gamma2s=mean_gamma2s+gamma2
          mean_rhos=mean_rhos+rho
        }
      }
    }
  }


  # return answer
  if(output_type=='chain'){
    eta.mcmc=etas[,burn_in:iter]
    eta1.mcmc=eta1s[,burn_in:iter]
    eta2.mcmc=eta2s[,burn_in:iter]
    beta1.mcmc=beta1s[,burn_in:iter]
    beta2.mcmc=beta2s[,burn_in:iter]
    gamma1.mcmc=gamma1s[,burn_in:iter]
    gamma2.mcmc=gamma2s[,burn_in:iter]
    rho.mcmc=rhos[,burn_in:iter]
    ans<-list(eta=eta.mcmc,eta1=eta1.mcmc,eta2=eta2.mcmc,
              beta1=beta1.mcmc,beta2=beta2.mcmc,gamma1=gamma1.mcmc,
              gamma2=gamma2.mcmc,rho=rho.mcmc,psi11=psi11,psi22=psi11)
  }else{
    if(output_type=='mean'){
    mean_etas=mean_etas/(iter-burn_in)
    mean_eta1s=mean_eta1s/(iter-burn_in)
    mean_eta2s=mean_eta2s/(iter-burn_in)
    mean_beta1s=mean_beta1s/(iter-burn_in)
    mean_beta2s=mean_beta2s/(iter-burn_in)
    mean_gamma1s=mean_gamma1s/(iter-burn_in)
    mean_gamma2s=mean_gamma2s/(iter-burn_in)
    mean_rhos=mean_rhos/(iter-burn_in)
    ans=list(eta=mean_etas,eta1=mean_eta1s,eta2=mean_eta2s,
             beta1=mean_beta1s,beta2=mean_beta2s,gamma1=mean_gamma1s,
             gamma2=mean_gamma2s,rho=mean_rhos,psi11=psi11,psi22=psi11)
    }
  }

  return(ans)
}