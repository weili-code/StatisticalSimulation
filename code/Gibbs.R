################
# Example (Gibbs sampler: Bivariate distribution)
################
# initialize constants and parameters
N <- 5000               #length of chain
burn <- 1000            #burn-in length
X <- matrix(0, N, 2)    #the chain, a bivariate sample

rho <- -.75             #correlation
mu1 <- 0
mu2 <- 2
sigma1 <- 1
sigma2 <- .5
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2

#generate the chain
X[1, ] <- c(mu1, mu2)            #initialize
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]

# compare sample statistics to parameters
colMeans(x)
cov(x)
cor(x)

mu<-c(mu1,mu2)
Sigma<-matrix(c(sigma1^2, rho*sigma1*sigma2,rho*sigma1*sigma2,sigma2^2),2,2)
library(MASS)
z<-mvrnorm(n = N, mu, Sigma, tol = 1e-6, empirical = FALSE)

par(mfrow=c(1,2)) 
plot(x, main="sample from Gibbs", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(x[,2]))
plot(z, main="sample from R", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(x[,2]))



################
# Example (Slice sampling -- 2D slice sampler) R.C. Exercise 7.11
################

# sample from f(x)=1/2 e^{-sqrt(x)}

T=5000
f=function(x){
  1/2*exp(-sqrt(x))}

X=c(runif(1)); U=c(runif(1))
for (t in 1:T){
  U=c(U,runif(1,0,f(X[t])))
  X=c(X,runif(1,0,(log(2*U[t+1]))^2))
}

par(mfrow=c(1,2))
hist(X,prob=TRUE,col="wheat2",xlab="",main="")
acf(X)

################
# Example (Slice sampling -- generalized)
################

# Example of a 3D slice sampler
# exp -x^2/2 (1+sin^2 3x) (1+cos^4 5x)
x=rep(0,5000)
for (i in 2:5000){
  omega1=(1+sin(3*x[i-1])^2)*runif(1)
  omega2=(1+cos(5*x[i-1])^4)*runif(1)
  omega3= runif(1)*exp(-x[i-1]^2/2)
  temp=sqrt(-2*log( omega3))
  repeat{
    y=-temp+2*temp*runif(1) # proposal for x
    if ((sin(3*y)^2>omega1-1)&(cos(5*y)^4>omega2-1)) break  
  }
  x[i]=y
}

par(mfrow=c(1,1)) 
# histogram of the generated samples
hist_plot=hist(x,breaks=75,col="grey",proba=T,xlab="x",ylab="",main="",
               freq=FALSE)
labs=seq(-3,3,.01)

# To overlay the histogram with a (crude) normalized density
dense=(1+sin(3*labs)^2)*(1+cos(5*labs)^4)*exp(-labs^2/2)
dense_norm=dense*max(hist_plot$density)/max(dense)
# hist_plot$density is the estimated density value 
# (relative frequency) from histogram (returned by freq=FALSE)
lines(labs,dense_norm,col="sienna4",lwd=2)
