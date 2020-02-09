# CM
This package provides functionality to implement latent conjugate multivariate models for non-Gaussian data.

For examples of how to use the package, please see the help files. Please note that CM is under active development.

# References

Bradley, JR, Holan, SH, and Wikle, CK. (2020). Bayesian Hierarchical Models with Conjugate Full-Conditional Distributions for Dependent Data from the Natural Exponential Family. *Journal of the American Statistical Association.*

Bradley, JR, Wikle, CK, and Holan, SH. (2019). Spatio-Temporal Models for Big Multinomial Data using the Conditional Multivariate Logit-Beta Distribution. *Journal of Time Series Analysis.* 40: 363 - 382.

Bradley, JR, Holan, SH, and Wikle, CK. (2018). Computationally Efficient Distribution Theory for Bayesian Inference of High-Dimensional Dependent Count-Valued Data (with discussion). *Bayesian Analysis.* 13: 253 - 302. *Rejoinder:* 302 - 310.

Hu, G, Bradley, JR. (2018). A Bayesian Spatio-Temporal Model for Analyzing Earthquake Magnitudes. *Stat.* 7(1): e179.

# Installation
```
#install.packages("devtools",dependencies=TRUE) #install devtools from CRAN
#library(devtools) #load devtools
#install_github("hadley/devtools")#install devtools from Github
library(devtools)
install_github("JonathanBradley28/CM")
library(CM)
```
# Examples 
The help files contain small simulation examples. For example, if you type help(GibbsPoissonMLG) or ??GibbsPoissonMLG, the help file will contain the following code.
```
set.seed(123)
#define a test function
#A non-linear test function
lambda <- function(t) exp(1.1 + sin(2 * pi * t))

#define some 1-d locations
points = seq(0,1,length.out=1001)
points=points[2:1001]
m = dim(as.matrix(points))[1]

#get the true mean at these locations
truemean<-matrix(0,m,1)
for (j in 1:length(points)){
  truemean[j,1] = lambda(points[j])
}

#simulate the data
data = matrix(0,m,1)
for (i in 1:m){
  data[i] = rpois(1,truemean[i])
}

#see how many zeros there are
sum(data==0)

#plot the data
plot(data,xlab="time",ylab="Poisson counts",main="Counts vs. time")

#covariate intercept-only
X = matrix(1,m,1)
p <- dim(X)[2]

##compute the basis function matrix
#compute thin-plate splines
r = 8
knots = seq(0,1,length.out=r)
#orthogonalize G
G = THINSp(as.matrix(points,m,1),as.matrix(knots,r,1))
outG<-qr(G)
G<-qr.Q(outG)
#orthogonalize X
outX<-qr(X)
X<-qr.Q(outX)

#Run the MCMC algorithm
output<-GibbsPoissonMLG(Niter=2000,X,G,data)
```

Peform MCMC diagnostics using the R package "code." Note that alpha_delta does not mix well. Of course, one could thin the MCMC to reduce autocorrelation, however, we choose not to. See MacEachern and Berliner (1994) for a discussion on thinning Markov chains.

```
#trace plots (without burnin)
plot(as.mcmc(output$betas[1000:2000]))
plot(as.mcmc(output$etas[1,1000:2000]))
plot(as.mcmc(output$etas[8,1000:2000]))
plot(as.mcmc(output$deltas[10,1000:2000]))
plot(as.mcmc(output$alpha_eta[1,1000:2000]))
plot(as.mcmc(output$alpha_delta[1,1000:2000]))


#estimates (remove a burnin)
lambda_est = apply(output$lambda_rep[,1000:2000],1,mean)
lambda_lower= apply(output$lambda_rep[,1000:2000],1,quantile,0.025)
lambda_upper= apply(output$lambda_rep[,1000:2000],1,quantile,0.975)


#plot estimates and truth
plot(1:m,truemean,ylim = c(0,max(lambda_upper)+1))
lines(1:m,lambda_est,col="red")
lines(1:m,lambda_lower,col="blue")
lines(1:m,lambda_upper,col="blue")

#smooth estimates (remove a burnin)
lambda_est = apply(output$lambda_rep_smooth[,1000:2000],1,mean)
lambda_lower= apply(output$lambda_rep_smooth[,1000:2000],1,quantile,0.025)
lambda_upper= apply(output$lambda_rep_smooth[,1000:2000],1,quantile,0.975)

#plot smooth estimates and truth
plot(1:m,truemean,ylim = c(0,max(lambda_upper)+1))
lines(1:m,lambda_est,col="red")
lines(1:m,lambda_lower,col="blue")
lines(1:m,lambda_upper,col="blue")


covmat = matrix(0,1000,1000)
for (jj in 1:1000){
  covmat = covmat+(output$lambda_rep[,1000+jj] - lambda_est)%*%t((output$lambda_rep[,1000+jj] - lambda_est))/1000
  print(jj)
}
vars = 1/sqrt(diag(covmat))
corrmat= diag(vars)%*%covmat%*% diag(vars)
image(corrmat)

```

To access examples for Binomial, Bernoulli, and Multinomial data, copy and paste from the help files

```
help(GibbsBinomialMLB)
help(GibbsBernoulliMLB)
help(GibbsMultinomialMLB)
help(GibbsPoissonMLG)
```
