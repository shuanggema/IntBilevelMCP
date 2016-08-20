dyn.load("sparse_MCP.so")

standard <- function (x) {
   p <- ncol(x)
   n <- nrow(x)
   x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}

sgroup <- function(x,y,x0,G,lambda1,lambda2,gamma1,gamma2,orthonormalize=F)
{
  n <- length(y)
  totalp <- dim(x)[2]
  p <- length(G)
  iter <- 1000

  pos_e <- cumsum(G)
  pos_s <- pos_e - G + 1

  if ( orthonormalize == F)
  {
    fun <- function(i,x)
    {
      sigma = t(x[,pos_s[i]:pos_e[i]])%*%x[,pos_s[i]:pos_e[i]]/n
      R = chol(sigma)  #cholesky decomp always exists for X'X/n where dj < n
      x.t = x[,pos_s[i]:pos_e[i]]%*%solve(R)
      return(list(R=R,x.t=x.t))
    }
    result <- lapply(seq(1:p),fun,x) 
    x.t <- numeric()
    for ( j in 1:p )
      x.t <- cbind(x.t,result[[j]]$x.t)
  }
  else if ( orthonormalize == T)
    x.t = x 

  d <- rep(1,totalp)
  beta <- numeric(totalp)
  param <- c(n, p, totalp,iter)
  diff <- 1 
  epsilon <- 1E-10
  y_c = y-mean(y)
  fit <- .C("sparse_MCP", y=as.double(y_c),x=as.double(t(x.t)),G=as.integer(G),d=as.double(d),
                          param=as.integer(param),tuning=as.double(c(lambda1,lambda2,gamma1,gamma2)),
                          epsilon=as.double(epsilon),beta=as.double(beta),diff=as.double(diff))

  r0 <- y - x.t%*%fit$beta
  inv <- solve(t(x0)%*%x0,diag(dim(x0)[2]))
  beta0 <- inv%*%t(x0)%*%r0
                 
  beta <- numeric()
  b = fit$beta
  if ( orthonormalize == F)
  {
    for ( j in 1:p )
      beta <- c(beta,solve(result[[j]]$R,b[pos_s[j]:pos_e[j]]))
  }
  else if ( orthonormalize == T)
    beta <- b

  beta[abs(beta)<1e-5]=0
  list(beta=beta,beta0=beta0,diff=fit$diff)
}

