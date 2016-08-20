dyn.load("bi_level.so")

Bilevel_Binary <- function(x,y,x0,ps,pe,group,lambda1,lambda2,gamma1,gamma2)
{
  n <- length(y)
  totalp <- dim(x)[2]
  p <- length(group)
  iter <- 20
  m <- dim(x0)[2]

  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1

  x.t = x 

  beta <- numeric(totalp)
  beta0 <- numeric(m)
  param <- c(n, p, totalp,m,500) 
  epsilon <- 1E-10
  h <- numeric(1)
  z_beta0 <- numeric(n)
  v <- rep(1,n)

  tuning <- c(lambda1,lambda2,gamma1,gamma2)

  repeat
  {
    pi.x  <- 1/(1+ exp( - x0%*%beta0 - x.t%*%beta))
    for ( i in 1:n )
    { 
      if (pi.x[i] <0.01) pi.x[i]=0.01
      else if (pi.x[i]>0.99) pi.x[i]=0.99
    }
    z <- x0%*%beta0 + x.t%*%beta + (y-pi.x)/(pi.x*(1-pi.x))
    weight <- pi.x*(1-pi.x) 
    for ( i in 1:m)
      z_beta0[ps[i]:pe[i]] <- z[ps[i]:pe[i]] - mean(z[ps[i]:pe[i]])
    
    fit <- .C("Bi_level_MCP_MCP", y=as.double(z_beta0),z=as.double(t(x0)),x=as.double(t(x.t)),w=as.double(weight),G=as.integer(group),param=as.integer(param),tuning=as.double(tuning),epsilon=as.double(epsilon),alpha=as.double(beta0),beta=as.double(beta))

    old.beta0 <- beta0
    old.beta <- beta
    r0 <- z - x.t%*%fit$beta
    inv <- solve(t(x0)%*%x0,diag(dim(x0)[2]))
    beta0 <- inv%*%t(x0)%*%r0
    beta <- fit$beta
    diff <- sum((beta - old.beta)^2)
    h <- h + 1
    if ( h == iter | diff < epsilon )
    {
      break
    }  
  }
  list(beta0=beta0,beta=beta,diff=diff)
}

