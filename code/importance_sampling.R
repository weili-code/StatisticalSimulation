#########
#example 3.5 of R.C. importance sampling estimators for P(Z>4.5)
########
pnorm(-4.5,log=T)

Nsim=10^5
y<-rexp(Nsim)+4.5
weit<-dnorm(y)/dexp(y-4.5)
plot(cumsum(weit)/1:Nsim,type="l")
abline(a=pnorm(-4.5),b=0,col="red")


#########
#example 
########

# X pdf proportional to exp(-x^2/2)*((sin(6*x))^2+3*((cos(x))^2)*(sin(4*x))^2+1)
# goal: estimate E(X), and draw from X
# Use: g~N(0,1)

fp <- function (x) {exp(-x^2/2)*((sin(6*x))^2+3*((cos(x))^2)*(sin(4*x))^2+1) }
plot(fp,-3,3)

N<-50000

# use self-normalized importance sampling
w<-function(x){ fp(x)/dnorm(x) }
x<-rnorm(N,0,1)
est<-cumsum(w(x)*x)/cumsum(w(x))
plot(est,type="l")
tail(est,n=10)

# obtain samples from self-normalized importance sampling
w_hat<-w(x)/sum(w(x))
x_sim<-sample(x, size=2500,replace=T,prob=w_hat)
hist(x_sim, breaks=100,freq=T)

# The approximate pdf of X
hist(x_sim, breaks=100,freq=F)
lines(density(x_sim))

# compare the pdf with the its kernel 
plot(density(x_sim))
curve(fp, add=TRUE)


#########
#example Voss E.6
########

# X~N(0,1), Find prob(X in [3,4])
# using four different importance pdfs
# Y~N(1,1), N(2,1), N(3.5,1) and exp(1)+3
# compare the means and variances of these four IS estimators

set.seed(95)

GetISEstimate <- function(N, gen.Y, psi) {
  # psi is the pdf of Y
  Y <- gen.Y(N)
  phi <- function(x) dnorm(x, 0, 1) # pdf of X
  weighted.samples <- (Y >= 3 & Y <= 4) * phi(Y) / psi(Y)
  return(list(p=mean(weighted.samples),
              var=var(weighted.samples)))
}

# try Y~N(1,1)
gen.Y1 <- function(N) { rnorm(N, 1, 1) }
psi1 <- function(x) { dnorm(x, 1, 1) }
GetISEstimate(10000, gen.Y1, psi1)

PrintISVariance <- function(name, gen.Y, psi) {
  est <- GetISEstimate(10000, gen.Y, psi)
  cat(name, "  -->  Var = ", est$var, "\n", sep="")
  return(est$var)
}

gen.Y2 <- function(N) { rnorm(N, 2, 1) }
psi2 <- function(x) { dnorm(x, 2, 1) }
gen.Y3 <- function(N) { rnorm(N, 3.5, 1) }
psi3 <- function(x) { dnorm(x, 3.5, 1) }
gen.Y4 <- function(N) { rexp(N) + 3 }
psi4 <- function(x) { dexp(x-3) }

PrintAllVariances <- function() {
  v1 <- PrintISVariance("N(1,1)  ", gen.Y1, psi1)
  v2 <- PrintISVariance("N(2,1)  ", gen.Y2, psi2)
  v3 <- PrintISVariance("N(3.5,1)", gen.Y3, psi3)
  v4 <- PrintISVariance("Exp(1)+3", gen.Y4, psi4)
  return(c(v1, v2, v3, v4))
}

res <- PrintAllVariances()
