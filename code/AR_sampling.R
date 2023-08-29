########
#example 2.1 of R.C. Sample from exp(1) using probability integral transform
########
set.seed(2021)
Nsim<-10^4 #number of random variables
U<-runif(Nsim)
X<--log(U) #transforms of uniforms
Y<-rexp(Nsim) #exponentials from R
par(mfrow=c(1,2)) #plots
hist(X,freq=F,main="probability integral transform")
hist(Y,freq=F,main="from R")



#########
#example 2.7 of R.C. Sample from beta(2.7,6.3) using accept-and-reject method
#########
set.seed(2021)
optimize(f=function(x){dbeta(x,2.7,6.3)},interval=c(0,1),maximum=T)$objective
Nsim<-2500
a<-2.7;b<-6.3
M<-2.67
u<-runif(Nsim,max=M) #uniform over (0,M)
y<-runif(Nsim) #generation from g
x<-y[u<dbeta(y,a,b)] #accepted subsample
z<-rbeta(Nsim,2.6,6.3)
par(mfrow=c(1,2)) #plots
hist(x,freq=F,main="from Accept-Reject")
hist(z,freq=F,main="from R")

output_num<-length(x)

#########
# Sample from beta(2.7,6.3) using the fundamental theorem
#########
set.seed(2021)
optimize(f=function(x){dbeta(x,2.7,6.3)},interval=c(0,1),maximum=T)$objective
Nsim<-2500
a<-2.7;b<-6.3
M<-2.67
y_accept<-vector()
num_accept<-0

while (num_accept<output_num) {
  y<-runif(n=1) #generation from g
  u<-runif(n=1, min=0, max=M*1) # g(y)=1 throughout
  if (u<dbeta(y,a,b)) {
    y_accept<-c(y_accept, y)
    num_accept<-num_accept+1
  }
  print(num_accept)
}

z<-rbeta(Nsim,2.6,6.3)
par(mfrow=c(1,2)) #plots
hist(y_accept,freq=F,main="from Accept-Reject")
hist(z,freq=F,main="from R")