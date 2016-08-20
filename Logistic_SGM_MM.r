dyn.load("sparse_MCP.so")

standard <- function (x) {
   p <- ncol(x)
   n <- nrow(x)
   x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}

d_cal <- function(x,pos_s,pos_e)
{
  p <- length(pos_s)
  n=nrow(x)
  index <- seq(1,p)
  trans<-function(i,x,pos_s,pos_e)
  {
    n=nrow(x)
    x_t = matrix(x[,pos_s[i]:pos_e[i]],nrow=n)
    d=diag(t(x_t)%*%x_t)/n
    list(d=d)
  }
  D=unlist(lapply(index,trans,x,pos_s,pos_e))
 
  list(D=as.vector(D))
}

sgroup <- function(x,y,x0,G,n_i,ps,pe,lambda1,lambda2,gamma1,gamma2)
{
  n <- length(y)
  totalp <- ncol(x)
  p <- length(G)
  iter <- 20
  m <- ncol(x0)

  pos_e <- cumsum(G)
  pos_s <- pos_e - G + 1

  beta <- numeric(totalp)
  beta0 <- numeric(m)
  param <- c(n, p, totalp,1000)
  epsilon <- 1E-10
  h <- numeric(1)
  z_beta0 <- numeric(n)
  d=rep(n_i,p)/n
  dif <- 1000

  repeat
  {
    pi.x  <- 1/(1+ exp( - x0%*%beta0 - x%*%beta))
    for ( i in 1:n )
    { 
      if (pi.x[i] <0.001) pi.x[i]=0.001
      else if (pi.x[i]>0.999) pi.x[i]=0.999
    }
    z <- x0%*%beta0 + x%*%beta + (y-pi.x)/(pi.x*(1-pi.x))
    weight <- pi.x*(1-pi.x) 
    for ( i in 1:m)
      z_beta0[ps[i]:pe[i]] <- z[ps[i]:pe[i]] - mean(z[ps[i]:pe[i]])
    
    fit <- .C("sparse_MCP", y=as.double(z_beta0),x=as.double(t(x)),G=as.integer(G),d=as.double(d),
                          param=as.integer(param),tuning=as.double(c(lambda1,lambda2,gamma1,gamma2)),
                          epsilon=as.double(epsilon),beta=as.double(beta),diff=as.double(1))

    old.beta0 <- beta0
    old.beta <- beta
    r0 <- z - x%*%fit$beta
    inv <- solve(t(x0)%*%x0,diag(dim(x0)[2]))
    beta0 <- inv%*%t(x0)%*%r0
    beta <- fit$beta
    diff <- sum((beta - old.beta)^2)
    d2=tail(dif, n=1)-diff
    dif <- c(dif, diff)
    h <- h + 1
    if ( h == iter | diff < epsilon | d2 < 0)
    {
      break
    }  
  }

  beta[abs(beta)<1e-5]=0
  list(beta=beta,beta0=beta0,dif=dif)
}