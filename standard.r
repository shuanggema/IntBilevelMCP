
standardize<-function(j,x)
{
  x.t   <- x[,j]
  x.non <- x.t[x.t!=0]
  N <- length(x.t)
  n <- length(x.non)
  x.mean <- mean(x.non)
  x.std <- (x.non-x.mean)/sd(x.non)*N^0.5/(n-1)^0.5
  x.t[x.t!=0] = x.std
  list(x.t=x.t)
}

