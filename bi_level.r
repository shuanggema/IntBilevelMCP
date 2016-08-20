dyn.load("bi_level.so")

standard <- function (x) 
{
   p <- ncol(x)
   n <- nrow(x)
   x.mean <- matrix(rep(apply(x,2,mean),n),n,p,byrow=T)
   x.std <- t(t((x-x.mean))/(apply(x,2,sd)*(n-1)^0.5)*n^0.5)
}

bi_level <- function(x,y,z,w,group,tuning)
{
  n <- dim(x)[1]
  p <- dim(x)[2]
  d <- length(group)
  m <- dim(z)[2]
  iter <- 1000
  param <- c(n,d,p,m,iter)
  x.std <- standard(x)
  alpha <- numeric(m)
  beta <- numeric(p)
  epsilon <- 1e-10

  fit <- .C("Bi_level_MCP_MCP", y=as.double(y),z=as.double(t(z)),x=as.double(t(x.std)),w=as.double(w),G=as.integer(group),param=as.integer(param),tuning=as.double(tuning),epsilon=as.double(epsilon),alpha=as.double(alpha),beta=as.double(beta))

  list(beta=fit$beta,alpha=fit$alpha)
}

