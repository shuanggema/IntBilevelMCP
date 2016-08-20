library(mvtnorm)
generate <- function(k,n,p)
{
 beta1 <- c(rep(0,45),c(0.7,1.1,1.3,1.1,0.9),rep(0,50),c(0.9,1.3,0.9,1.7,0.3),rep(0,p-105)) 
 beta2 <- c(rep(0,45),c(0.5,1.2,1.1,0.9,0.8),rep(0,50),c(1.1,0.4,1.3,0.5,1.2),rep(0,p-105)) 
 beta3 <- c(rep(0,45),c(1.1,0.8,0.9,1.2,1.2),rep(0,35),c(1.2,0.7,1.1,0.7,1.3),rep(0,p-90))
 beta4 <- c(rep(0,45),c(1.2,0.9,0.6,1.1,1.3),rep(0,35),c(0.8,1.5,1.9,0.6,0.8),rep(0,p-90))

  
# m = 4
 covar <- matrix(rep(0,p*p),nrow=p)
 
 for ( i in 1:p)
  for ( j in 1:p)
   covar[i,j]=0.2^(abs(i-j))

  for (i in 1:p)
  {
    covar[i,i]=1
  }


  x1 <- rmvnorm(n = n, mean = rep(0, p), covar) 
  x2 <- rmvnorm(n = n, mean = rep(0, p), covar)
  x3 <- rmvnorm(n = n, mean = rep(0, p), covar)
  x4 <- rmvnorm(n = n, mean = rep(0, p), covar)

  lp1 <- x1%*%beta1+rnorm(n = n, 0 , k)
  lp2 <- x2%*%beta2+rnorm(n = n, 0 , k)
  lp3 <- x3%*%beta3+rnorm(n = n, 0 , k)
  lp4 <- x4%*%beta4+rnorm(n = n, 0 , k)

  y1 <- rbinom(length(lp1),1,plogis(lp1))
  y2 <- rbinom(length(lp2),1,plogis(lp2))
  y3 <- rbinom(length(lp3),1,plogis(lp3))
  y4 <- rbinom(length(lp4),1,plogis(lp4))

  list(x1=x1,x2=x2,x3=x3,x4=x4,y1=y1,y2=y2,y3=y3,y4=y4)
}





