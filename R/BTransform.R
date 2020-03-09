#' Bayesian Transformations
#'
#' This code posterior replicates from the transformation model described in Bradley et al. (2020).
#' @param B The number of iterations of the Gibbs sampler
#' @param data_continuous Normal distributed data vector
#' @param data_poisson Poisson distributed data vector
#' @param data_multinomial Binomial distributed data vector
#' @param nn The sample size assocated with the binomial observations. Must be same length of data_multinomial
#' @param nslice The burnin for the slice sampler
#' @param report Option to print the iteration number of the MCMC.
#' @import actuar MfUSampler stats coda MASS progress Matrix LearnBayes
#' @return ans A list of updated parameter values.
#' @examples
#'
#'library(CM)
#'
#'Xu = matrix(runif(10*1000),1000,9)
#'Xt = Xu
#'Xt[,1] = sin(pi*Xu[,1]*Xu[,2])
#'Xt[,2] = (X[,3] - 0.5)^2
#'beta = as.matrix(c(10,20,10,5,0,0,0,0,0),10,1)
#'f<-Xt%*%beta
#'
#'Xmat = as.matrix(cbind(Xt[,1],Xt[,3:9]))
#'Xmat1 = cbind(Xmat[1:350,],matrix(0,350,18))
#'Xmat2 = cbind(matrix(0,350,9),Xmat[351:700,],matrix(0,350,9))
#'Xmat3 = cbind(matrix(0,300,18),Xmat[701:1000,])
#'
#'Xmat = rbind(Xmat1,Xmat2,Xmat3)
#'
#'X1=matrix(1,1000,1)
#'X2=matrix(0,1000,1)
#'X3=matrix(0,1000,1)
#'X2[351:700]=1
#'X3[701:1000]=1
#'
#'Xmat = cbind(X1,X2,X3,Xmat)
#'
#'N=1000
#'r = 500
#'A = matrix(0,N,N)
#'A[1:(N-1),2:N] = diag(N-1)
#'A = A+t(A)
#'MI = (diag(N)-Xmat%*%ginv(t(Xmat)%*%Xmat)%*%t(Xmat))%*%A%*%(diag(N)-Xmat%*%ginv(t(Xmat)%*%Xmat)%*%t(Xmat))
#'output<-eigen(MI)
#'Psi=Re(output$vectors[,1:r])
#'
#'y = f[1:350] + rnorm(350)
#'y2 = rpois(350,(f[351:700]))
#'y3= rbinom(300, 31, (f[701:1000])/31)
#'
#'outCM<-BTransform(2000,y,y2,y3,matrix(31,300,1))
#'
#'response = rbind(outCM$xi1,outCM$xi2,outCM$xi3)
#'## Fit the preferred model
#'
#'etadiscrip = matrix(0,r,2000)
#'betadiscrip = matrix(0,29,2000)
#'xidiscrip=matrix(0,1000,2000)
#'
#'GibbsNormalGAU<-function(B,X,G,Z,etas2=NA,betas2=NA,xi=NA,taus2=NA,Hphib2=NA,sig2xi=NA,printevery = 100){
#'
#'  outLGP=GibbsNormalGAU(2,Xmat,Psi,(response[,1]),matrix(0,r,1),matrix(0,29,1),matrix(0,1000,1),1,diag(r),1,1)
#'  betadiscrip[,1]=outLGP[[1]][,2]
#'  etadiscrip[,1]=outLGP[[2]][,2]
#'  xidiscrip[,1]=outLGP[[3]][,2]
#'
#'  for (j in 2:2000){
#'    outLGP=GibbsNormalGAU(2,Xmat,Psi,response[,j],etadiscrip[,j-1],betadiscrip[j-1],xidiscrip[,j-1],outLGP[[4]][2],outLGP[[5]][2]*diag(r),outLGP[[6]][2],outLGP[[7]][2])
#'    betadiscrip[,j]=outLGP[[1]][,2]
#'    etadiscrip[,j]=outLGP[[2]][,2]
#'    xidiscrip[,j]=outLGP[[3]][,2]
#'    if (j%%100==0){
#'      print(cbind("MCMC Replicate:", j))
#'    }
#'  }
#'
#'
#'  yhat2=Xmat%*%betadiscrip[,500:1000]+Psi%*%etadiscrip[,500:1000]+xidiscrip[,500:1000]
#'
#'  #countinuous predictions with truth
#'  f_est=apply(yhat2[1:350,],1,mean)
#'  f_lower= apply(yhat2[1:350,],1,quantile,0.025)
#'  f_upper= apply(yhat2[1:350,],1,quantile,0.975)
#'
#'  #plot estimates and truth
#'  plot(f_est,f[1:350],ylim = c(0,max(f[1:350])+1))
#'  abline(0,1)
#'
#'  #Poisson count predictions with truth
#'f_est=apply(exp(yhat2[351:700,]),1,median)
#'f_lower= apply(exp(yhat2[351:700,]),1,quantile,0.025)
#'  f_upper= apply(exp(yhat2[351:700,]),1,quantile,0.975)
#'
#'  #plot estimates and truth
#'  plot(f_est,f[351:700,],ylim = c(0,max(f[351:700,])+1))
#'  abline(0,1)
#'
#'  #Binomial count predictions with truth
#'  f_est=31*apply(exp(yhat2[701:1000,])/(1+exp(yhat2[701:1000,])),1,median)
#'  f_lower= apply(31*exp(yhat2[701:1000,])/(1+exp(yhat2[701:1000,])),1,quantile,0.025)
#'  f_upper= apply(31*exp(yhat2[701:1000,])/(1+exp(yhat2[701:1000,])),1,quantile,0.975)
#'
#'   #plot estimates and truth
#' plot(f_est,f[701:1000],ylim = c(0,max(f_est)+1))
#' abline(0,1)
#' @export
BTransform<-function(B,data_continuous,data_poisson,data_multinomial,nn,report=100,nslice=2){


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





 n1 = length(data_continuous)
 n2 = length(data_poisson)
 n3 = length(data_multinomial)


logit<-function(WW){
return(log(WW/(1-WW)))
}


  #run the Gibbs sampler from Section 2.3
  xi1.gibbs = matrix(0,n1,B)
  alphax1 = matrix(1,1,B)
  kappax1 = matrix(2,1,B)
  vvv=matrix(1,n1,1)

  xi2.gibbs = matrix(0,n2,B)
  alphax2 = matrix(1,1,B)
  kappax2 = matrix(2,1,B)

alpha_eps=matrix(1,1,B)
kappa_eps=matrix(2,1,B)


  xi3.gibbs = matrix(0,n3,B)
  alphax3 = matrix(1,1,B)
  kappax3 = matrix(2,1,B)


  print("Collapsed Gibbs Sampler Starting")
  for (b in 2:B){

if (is.null(data_continuous)==FALSE){
    #update xi
    alpha = alpha_eps[b-1]*matrix(1,n1,1)+alphax1[,b-1]
    kappa = kappa_eps[b-1]*matrix(1,n1,1)+kappax1[,b-1]
    xi1.gibbs[,b] = sqrt(vvv)*rnorm(n1)  + data_continuous

    if (sum(is.nan(xi1.gibbs[,b]))>0){
      xi1.gibbs[,b] = xi1.gibbs[,b-1]
    }

    if (sum(is.infinite(xi1.gibbs[,b]))>0){
      xi1.gibbs[,b] = xi1.gibbs[,b-1]
    }

    #update shape xi using slice sampler
    for (iv in 1:n1){
vvv[iv] = rigamma(1,0.5+n1/2,0.5+(xi1.gibbs[iv,b]^2)/2)
}

}



if (is.null(data_poisson)==FALSE){
    #### update xi2 ####
    alphadeltas=data_poisson+alphax2[,b]
    kappadeltas=1*matrix(1,n2,1)
    xi2.gibbs[,b]=log(rgamma(n2,alphadeltas,rate=kappadeltas))

	#update xi2
#	for (iii in 1:n2){
f<-function(x){
      like= n2*x*log(0.5)- (n2)*(lgamma(x)) +x*sum(xi2.gibbs[,b]) -x
      if (x<0.01){
        like = -Inf;
      }

      if (x>20000){
        like = -Inf;
      }

      return(like)

      }

    alpha_e_temp<-MfU.Sample.Run(alphax2[b-1], f, nsmp = nslice)
    alphax2[b] = abs(alpha_e_temp[nslice])
	#}
}


if (is.null(data_multinomial)==FALSE){
    #update xi3
    alpha = data_multinomial+alphax3[,b-1]
    kappa = nn+2*alphax3[,b-1]
    xi3.gibbs[,b] = logitbetasim(alpha,kappa)

#for (iii in 1:n3){
#update shape xi using adaptive rejection
    f<-function(x){
      like=(n3*(lgamma(x)) - n3*(lgamma(x))-x +x*sum(xi3.gibbs[,b]))
      if (x<0.01){
        like = -Inf;
      }


      if (x>20000){
        like = -Inf;
      }

      return(like)

    }
    alphax_temp<-MfU.Sample.Run(alphax3[b-1], f, nsmp = nslice)
    alphax3[b] =abs(alphax_temp[nslice])

#}
}

    if (b>(report-1)){
    if ((max(b,report) %% min(b,report)) ==0){
      print(paste("MCMC Replicate: ",b))
    }
    }

  }


 output<-list(xi1=xi1.gibbs,xi2=xi2.gibbs,xi3=xi3.gibbs,a1=alphax1,a2=alphax2,a3=alphax3)
 return(output)
}
