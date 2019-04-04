#install.packages("devtools") #uncomment this line if you have not installed "devtools"

library(devtools)

#the author token will be removed once the repository  goes public.

install_github("JonathanBradley28/CM",auth_token="7fcfcc652f9d3b1441a8050c13d44f50abb65eaa")

library(CM)

###################
#The help file contains simulation examples. For example, if you type help(GibbsPoisson) or ??GibbsPoisson, the help file will contain the following code

set.seed(123)

#define a test function

#A non-linear test function

lambda <- function(t) exp(1.1 + sin(2 &ast; pi &ast; t))




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





#plot data and truth

plot(1:m,truemean,ylim = c(0,max(lambda_upper)+1))

lines(1:m,lambda_est,col="red")

lines(1:m,lambda_lower,col="blue")

lines(1:m,lambda_upper,col="blue")

########################
#To access examples for Binomial, Bernoulli, and Multinomial data, copy and paste from the help files

help(GibbsBinomialMLB)

help(GibbsBernoulliMLB)

help(GibbsMultinomialMLB)

help(GibbsPoissonMLB)
