### This script contains some examples for factorization ###


###############################
### generate multivariate normal using Cholesky
###############################

mvrnormR <- function(nsim,mu, Sigma) { 
  # generate multivariate normal using cholesky decomposition
  # input: nsim=number of rows, mu=p-dim vector, Sigma=p by p 
  # variance matrix
  # output: a matrix of multivariate normals (nsim by p)
  p<-ncol(Sigma)
  out<-matrix(rnorm(nsim * p), ncol = p, byrow=FALSE) %*% chol(Sigma)
  out<-out+matrix(mu, nrow=nsim, ncol=p, byrow=TRUE)
  out
}

# ALWAYS SET-UP A SEED NUMBER before generating random numbers!
set.seed(1)
mu<-c(1,2,3)
Sigma<-diag(c(0.3,0.3,0.1))
output<-mvrnormR(10,mu,Sigma)

###############################
### solving OLS using QR
###############################

n <- 50
x <- seq(1,500,len=n)
X <- cbind(1,x,x^2,x^3)
colnames(X) <- c("Intercept","x","x2","x3")
beta <- matrix(c(1,1,1,1),4,1)

# ALWAYS SET-UP A SEED NUMBER before generating random numbers!
set.seed(1)

y <- X%*%beta+rnorm(n,sd=1)

# try to solve OLS
solve(t(X)%*%X) %*% t(X) %*%y
solve(t(X)%*%X, t(X) %*%y )

# why?
options(digits=4)
t(X)%*%X
log(t(X)%*%X)

# use QR
QR <- qr(X)
Q <- qr.Q( QR )
R <- qr.R( QR )
betahat <- backsolve(R, t(Q)%*%y )
betahat

QR <- qr(X)
betahat2 <- solve.qr(QR, y)
betahat2

###############################
### ### Computing PCA ### ##### 
###############################

library(stats) # use: cov()

# get some data:
x <- c(2.5,.5,2.2,1.9,3.1,2.3,2,1,1.5,1.1)
y <- c(2.4,0.7,2.9,2.2,3.0,2.7,1.6,1.1,1.6,.9)
plot(x,y,xlim=c(-1,4),ylim=c(-1,4)); abline(h=0,v=0,lty=3)


x1 <- x - mean(x)
y1 <- y - mean(y)
plot(x1,y1); abline(h=0,v=0,lty=3)

M <- matrix(c(x1,y1),ncol=2) # make a matrix of the given data
M

# compute sample covariances/covariances
# This gives t(M)%*%m/(nrow(M)-1) since columns of m are demeanded
cov.M <- cov(M)
cov.M  

cov.eig <- eigen(cov.M)
cov.eig

# Eigenvectors of cov.m (orthogonal)
# cov.eig$vectors give eigenvectors
cov.eig$vectors[,1] %*% cov.eig$vectors[,2]

# let's plot these eigenvectors onto the data to present the new basis
plot(x1,y1); abline(h=0,v=0,lty=3)
abline(a=0,b=(cov.eig$vectors[1,1]/cov.eig$vectors[2,1]),col="red")
abline(a=0,b=(cov.eig$vectors[1,2]/cov.eig$vectors[2,2]),col="green")

# construct principal components
V_sel=as.matrix(cov.eig$vectors[,c(1,2)],ncol=2) # select both PCs
Z=M%*%V_sel # PCs matrix
cov(Z) # note the covariance

# plot PC scores (for ann data points)
plot(Z[,1],Z[,2],ylim=c(-2,2));abline(h=0,v=0,lty=3)

# recover initial dataset cbind(x,y) from z
M2=Z%*%t(V_sel)
M2[,1]=M2[,1]+mean(x)
M2[,2]=M2[,2]+mean(y)
# check
M2-cbind(x,y)

### prcomp() and svd()
## prcomp() uses SVD to perform pricipal components analysis

# (1) get singular value decomposition of M (demeaned data matrix)
svd.M <- svd(M)
svd.M
# note that the V vectors from svd.mare same to that from cov.eig
# except for the signs

pca.M <- prcomp(M, center=FALSE, scale=FALSE)
pca.M

# When variables have different units (for columns), then it may 
# sense to scale them.
# so that all have unit variance (scale=TRUE)  (it uses correlation 
# matrix). This is necessary if the data has different  units.

svd.M <- svd(scale(M))
# ?scale
# scale(x, center = TRUE, scale = TRUE)
svd.M$v 
pca.M <- prcomp(m, center=TRUE, scale=TRUE)
pca.M
# The V matrix are different from above obtained if without scaling.

# pca.m$rotation gives the eigenvectors
# We check 
plot(x1,y1); abline(h=0,v=0,lty=3)
abline(a=0,b=(pca.m$rotation[1,1]/pca.m$rotation[2,1]),col="red")
abline(a=0,b=(pca.m$rotation[1,2]/pca.m$rotation[2,2]),col="green")

summary(pca.M)

