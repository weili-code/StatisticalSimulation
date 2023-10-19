
########
# Example: simulated annealing
# example 5.2, 5.9
#######
#a simulated sample of 400 observations from this mixture with mu1 = 0 and mu2 = 2:5
# be careful, some parts of the code in the book are incorrect.
N <- 400

components <- sample(1:2,prob=c(0.25,0.75),size=N,replace=TRUE)
mu <- c(0,2.5)
da <- rnorm(N)+mu[components]
hist(da)

# loglikelihood function
like=function(mu){ sum(log( .25*dnorm(da-mu[1])+.75*dnorm(da-mu[2]) ))}

mu1<-seq(from=-2,to = 5, length.out=250)
mu2<-seq(from=-2,to = 5, length.out=250)
lli<-matrix(, nrow = 250, ncol = 250)
for (i in 1:250){
  for (j in 1:250){
    lli[i,j]=like(c(mu1[i],mu2[j]))
  }
}


# plot 
persp(mu1,mu2,lli,theta=15,phi=20,col="green4",
      ltheta=-120,shade=.75,border=NA,box=FALSE)
# theta, phi angles defining the viewing direction. 
# theta gives the azimuthal direction and phi the colatitude.

# movable 3d plot
library(rgl)
persp3d(mu1,mu2,lli,theta=155,phi=30,col="green4",
        ltheta=-120,shade=.75,border=NA,box=FALSE)

# define Simulated Annealing Functions
SA=function(x){
  # x is the start point
  temp=1 
  scale=1
  iter=1
  dif=1
  
  the=matrix(x,ncol=2) # matrix to store theta values
  curlike=hval=like(x) #current likelihood value
  
  while (dif>10^(-4)){
    prop=the[iter,]+rnorm(2)*scale
    # if the draw falls outside the domain, or fails the accept rule
    # we will not accept prop
    if ((max(-prop)>2)||(max(prop)>5)||
        (temp[iter]*log(runif(1))>like(prop)-curlike))
        prop=the[iter,]
    curlike=like(prop)
    hval=c(hval,curlike)
    the=rbind(the,prop)
    iter=iter+1
    temp=c(temp,(1/1.5)^iter)
    #temp=c(temp,1/(100*log(iter+1)))
    
    # rule to stop iteration
    dif=(iter<3000)+(max(hval)-max(hval[1:(iter/2)]))
  }
  
  list(theta=the,like=hval,ite=iter)
}

# start point (3,1)
outcome<-SA(c(3,1))
outcome$ite
sols<-outcome$the
tail(sols,10)

# plot heatmap
image(mu1,mu2,lli,xlab=expression(mu[1]),ylab=expression(mu[2]))
contour(mu1,mu2,lli,nle=100,add=T)

# plot the iterations
for (j in 1:(outcome$ite-1)){
  arrows(sols[j,1], sols[j,2], sols[j+1,1], sols[j+1,2], length = 0.1, 
         angle = 30, code = 2, col = "orange")
}
sols[outcome$ite,]
# it moves to the correct mode in just few jumps.
