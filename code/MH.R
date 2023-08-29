
#############
# Example: simulate a beta distribution (R.C. example 6.1)
#############
a=2.7; b=6.3
Nsim=500
X=rep(runif(1),Nsim) # initialize the chain
for (i in 2:Nsim){
  Y=runif(1)
  r=dbeta(Y,a,b)/dbeta(X[i-1],a,b)
  X[i]=X[i-1] + (Y-X[i-1])*(runif(1)<r)
}

par(mfrow=c(1,1)) 
plot(X,type="l") #trace plot

Z<-rbeta(Nsim,2.6,6.3)
par(mfrow=c(1,2)) #plots
hist(X,freq=F,main="samples from MH")
hist(Z,freq=F,main="samples from R")

#############
# Example: simulate a N(0,1) distribution using Metropolis (R.C. example 6.4)
############
Nsim=5000
X=rep(runif(1),Nsim) # initialize the chain
delta=1 #try 1, 10
for (i in 2:Nsim){
  e=runif(1,min=-delta,max=delta)
  Y=X[i-1]+e
  r=exp((1/2)*(X[i-1]^2-Y^2))
  X[i]=X[i-1] + (Y-X[i-1])*(runif(1)<r)
}

par(mfrow=c(1,1)) 
plot(X,type="l") #trace plot

Z<-rnorm(Nsim)
par(mfrow=c(1,2)) #plots
hist(X,freq=F,main="samples from MH")
hist(Z,freq=F,main="samples from R")

# plot auto-correlation function 
acf(X,lag.max=30)

# thinning
X.thin<-X[seq(25,Nsim,by=10)]
acf(X.thin,lag.max=30)

par(mfrow=c(1,2)) #plots
hist(X.thin,freq=F,main="thinned amples from MH")
hist(Z,freq=F,main="samples from R")

#################
# Example: simulate a bivariate normal using Metropolis (joint update)
#################
# object bivariate density function up to a multiplicative constant
p.star <- function(x, mu, Sigma){
  exp(- 0.5 * (x - mu) %*% solve(Sigma) %*% (x - mu))
}
mu <- c(1, 2)
Sigma <- c(1, 1, 1, 4)
Sigma <- matrix(Sigma, 2, 2)
Sigma
#use a bivariate-rectangular proposal distribution centered at the preceding value xi−1
#with half-extent d = 2 in the direction of the coordinate x1 and d = 2 in the direction of x2.

N <- 10000 # draws
x.current <- c(0, 0) # x_0
X<- matrix(0, N, 2) # to hold sampled values
d<-2

rbvunif <- function(mu, delta){
  u1 <- runif(1, mu[1] - d, mu[1] + d)
  u2 <- runif(1, mu[2] - d, mu[2] + d)
  c(u1, u2)
}

for (i in 1:N){
  x.proposed <- rbvunif(mu=x.current, delta=delta) # proposal
  r <- p.star(x.proposed, mu, Sigma) / p.star(x.current, mu, Sigma) 
  if (r >= 1 || runif(1) <= r) {
    X[i, ] <- x.proposed
    x.current <- x.proposed
  }
  else {
    X[i, ] <- x.current
  }
}

# compare sample statistics to parameters
colMeans(X)
cov(X)

library(MASS)
Z<-mvrnorm(n = N, mu, Sigma, tol = 1e-6, empirical = FALSE)

par(mfrow=c(1,2)) 
plot(X, main="sample from Metropolis", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(X[,2]))
plot(Z, main="sample from R", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(X[,2]))

# burn-in and/or thinning may be used as well.

par(mfrow=c(1,2)) 
plot(X[,1],type="l")
plot(X[,2],type="l") 

# ACF can be plotted for each coordinate

#################
# Example: simulate a bivariate normal using Metropolis (blockwise update)
#################
# object bivariate density function up to a multiplicative constant
p.star <- function(x, mu, Sigma){
  exp(- 0.5 * (x - mu) %*% solve(Sigma) %*% (x - mu))
}
mu <- c(1, 2)
Sigma <- c(1, 1, 1, 4)
Sigma <- matrix(Sigma, 2, 2)
#use a bivariate-rectangular proposal distribution centered at the preceding value xi−1
#with half-extent d1 = 2 in the direction of the coordinate x1 and d2 = 4 in the direction of x2.

N <- 10000 # draws
x.current <- c(0, 0) # x_0
X<- matrix(0, N, 2) # to hold sampled values
d1<-2
d2<-4

for (i in 1:N){
  #update first component x1
  x1.proposed <- x.current[1]+runif(1,min=-d1,max=d1)
  x.proposed<-c(x1.proposed,x.current[2])
  r <- p.star(x.proposed, mu, Sigma) / p.star(x.current, mu, Sigma) 
  if (r >= 1 || runif(1) <= r) {
    X[i,1] <- x1.proposed
    x.current[1]<- x1.proposed
  }
  else {
    X[i,1] <- x.current[1]
  }
  
  #update second component x2
  x2.proposed <- x.current[2]+runif(1,min=-d2,max=d2)
  x.proposed<-c(x.current[1],x2.proposed)
  r <- p.star(x.proposed, mu, Sigma) / p.star(x.current, mu, Sigma) 
  if (r >= 1 || runif(1) <= r) {
    X[i,2] <- x2.proposed
    x.current[2]<- x2.proposed
  }
  else {
    X[i,2] <- x.current[2]
  }
  
}

# compare sample statistics to parameters
colMeans(X)
cov(X)

library(MASS)
Z<-mvrnorm(n = N, mu, Sigma, tol = 1e-6, empirical = FALSE)

par(mfrow=c(1,2)) 
plot(X, main="sample from Metropolis", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(X[,2]))
plot(Z, main="sample from R", cex=.5, xlab=bquote(X[1]), ylab=bquote(X[2]), ylim=range(X[,2]))