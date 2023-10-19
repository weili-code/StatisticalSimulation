

#####################
### Example 5.6, 5.7 stochastic gradient descent (local search)
####################


h=function(x,y){(x*sin(20*y)+y*sin(20*x))^2*cosh(sin(10*x)
            *x)+(x*cos(10*y)-y*sin(10*x))^2*cosh(cos(20*y)*y)}

x=y=seq(-2,2,le=435) #defines a grid for persp
z=outer(x,y,Vectorize(h))

# plot 
par(bg="wheat",mar=c(1,1,1,1)) #bg stands for background
persp(x,y,z,theta=155,phi=30,col="green4",
        ltheta=-120,shade=.75,border=NA,box=FALSE)

# movable 3d plot
library(rgl)
persp3d(x,y,z,theta=155,phi=30,col="green4",
      ltheta=-120,shade=.75,border=NA,box=FALSE)


# goal: to find minimizer of h


SG=function(x,alpha,beta){
  theta=matrix(x,ncol=2)
  diff=iter=1
  hval=h(x[1],x[2])
  
  while (diff>10^-5){

     zeta=rnorm(2)
     zeta=zeta/sqrt(as.vector(t(zeta)%*%zeta))
     h_plus=theta[iter,]+beta(iter)*zeta
     h_minus=theta[iter,]-beta(iter)*zeta
     grad=alpha(iter)*zeta*(h(h_plus[1],h_plus[2])-
                                h(h_minus[1],h_minus[2]))/(2*beta(iter))

     scale=sqrt(t(grad)%*%grad)
     while (scale>1){
       zeta=rnorm(2)
       zeta=zeta/sqrt(as.vector(t(zeta)%*%zeta))
       h_plus=theta[iter,]+beta(iter)*zeta
       h_minus=theta[iter,]-beta(iter)*zeta
       grad=alpha(iter)*zeta*(h(h_plus[1],h_plus[2])-
                                 h(h_minus[1],h_minus[2]))/(2*beta(iter))
       scale=sqrt(t(grad)%*%grad)
       print(paste0("At iteration " , iter, ", zeta needs to be regenerated...", sep = ''))
     }   
     
     theta=rbind(theta,theta[iter,]-grad) 
     hval=c(hval,h(theta[1],theta[2]))
     diff=scale
     if (iter==500){diff=0}
     
     print(paste0("At iteration " , iter, sep = ''))
     iter=iter+1
     
  }
  
  list(theta=theta,hval=hval,ite=iter)
}


alpha=function(x){1/(x+1)}
beta=function(x){1/((x+1)^{0.1})}


start=c(.65,.8)
outcome=SG(start,alpha,beta)
outcome$ite
sols<-outcome$the
tail(sols,10)

# plot heatmap
image(x,y,z)
contour(x,y,z,nle=100,add=T)

# plot the iterations
for (j in 1:(outcome$ite-1)){
  arrows(sols[j,1], sols[j,2], sols[j+1,1], sols[j+1,2], length = 0.1, 
         angle = 30, code = 2, col = "red")
}
sols[outcome$ite,]