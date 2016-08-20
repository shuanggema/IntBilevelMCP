library(snow)
source("Logistic_SGM_MM.r")

lambda2.max <- function(x,y,ps,pe)
{
  n <- nrow(x)
  totalp <- ncol(x)
  m <- length(ps)

  z_beta0 <- numeric(n)
  pi.x  <- rep(0.5,n)
  z <- 0 + (y-pi.x)/(pi.x*(1-pi.x))
  weight <- pi.x*(1-pi.x) 
  for ( i in 1:m)
    z_beta0[ps[i]:pe[i]] <- z[ps[i]:pe[i]] - mean(z[ps[i]:pe[i]])
      
  lambda.max = max(abs(colSums(as.vector(z_beta0)*x)))/n
  list(lambda2.max=lambda.max)
}

lambda1.max <- function(x,y,d,ps,pe,group,pos_s,pos_e,lambda2)
{
  p <- length(group)
  totalp <- ncol(x)
  n <- nrow(x)
  m <- length(ps)
  
  z_beta0 <- numeric(n)
  pi.x  <- rep(0.5,n)
  z <- 0 + (y-pi.x)/(pi.x*(1-pi.x))
  weight <- pi.x*(1-pi.x) 
  for ( i in 1:m)
    z_beta0[ps[i]:pe[i]] <- z[ps[i]:pe[i]] - mean(z[ps[i]:pe[i]])

  fun <- function(j,pos_s,pos_e,x,y,group,d,lambda2)
  {
    t=x[,pos_s[j]:pos_e[j]]
    dt=d[pos_s[j]:pos_e[j]]
    n=nrow(t)
    z=colSums(t*as.vector(y)/n)
    tmp=abs(z)-lambda2
    tmp[tmp<=0]=0
    s_value=sign(z)*tmp/dt^0.5
    inner = sum(s_value^2)^0.5/group[j]^0.5
    #inner=sum( abs(colSums(t*as.vector(y)/n)) ^2)^0.5/group[j]^0.5
    list(inner=inner)
  }
  
  lambda.max=max(unlist(lapply(seq(1,p),fun,pos_s,pos_e,x,z_beta0,group,d,lambda2)))
  list(lambda1.max=lambda.max)
}

fun <- function(i,x,y,x0,m,fold.d,fold.g,n.fold,group,lambda1,lambda2,gamma)
{
  delete = unlist(fold.d[i])
  x.train = x[-delete,]
  y.train = y[-delete]
  x.valid = x[delete,]
  y.valid = y[delete]
  x0.train=x0[-delete,]
  x0.valid=x0[delete,]

  Ni_cal <- function(x,pos_s,pos_e)
  {
    p <- length(pos_s)
    n=nrow(x)
    index <- seq(1,p)
    trans<-function(i,x,pos_s,pos_e)
    {
      n=nrow(x)
      x_t = matrix(x[,pos_s[i]:pos_e[i]],nrow=n)
      n_i=diag(t(x_t)%*%x_t)
      list(n_i=n_i)
    }
    N_i=unlist(lapply(index,trans,x,pos_s,pos_e))
 
    list(N_i=as.vector(N_i))
  } 

  n.t <- numeric(m)
  for ( k in 1:m )
  {
    n.t[k] = 0
    for ( j in 1:n.fold)
    {
      if ( j != i )
      n.t[k] = n.t[k] + length(fold.g[[k]][[j]])
    }
  }
  pe.t <- cumsum(n.t)
  ps.t <- pe.t - n.t + 1
  
  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1
  n_i=Ni_cal(x.train,pos_s,pos_e)$N_i

  source("Logistic_SGM_MM.r")
  fit=sgroup(x.train,y.train,x0.train,group,n_i,ps.t,pe.t,lambda1,lambda2,gamma,gamma)

  lp <- x0.valid%*%fit$beta0+x.valid%*%fit$beta
  nl=-sum(y.valid*lp - log(1+exp(lp)))
  list(nl=nl)
}

parCv <- function(cl,lambda1.seq,lambda2.seq,m,x,y,x0,group,fold.d,fold.g,n.fold,gamma)
{
  step <- length(lambda1.seq)
  PE <- numeric(step)
  index <- seq(1:n.fold)

  for( i in 1:step)
  {
    lambda1=lambda1.seq[i]
    lambda2=lambda2.seq[i]
    v <- clusterApply(cl, index,fun,x,y,x0,m,fold.d,fold.g,n.fold,group,lambda1,lambda2,gamma)
    PE[i] <- sum(unlist(v))/n.fold
  }
  list(PE=PE)
}

cv.optim <- function(x,y,x0,n_i,group,ps,pe,n.fold,epsilon1,epsilon2,n.step,gamma,cl)
{
  m <- length(ps)
  fold.g=list(NULL,NULL,NULL,NULL)
  fold.d=list(NULL,NULL,NULL,NULL,NULL)
  for ( k in 1:m )
     fold.g[[k]] <- splitList(sample(ps[k]:pe[k]),n.fold)
  for ( k in 1:n.fold )
  {
     a <- numeric(0)
     for ( j in 1:m )
     {
        a <- c(a,fold.g[[j]][[k]])
     }
     fold.d[[k]]=a
  }

  pos_e <- cumsum(group)
  pos_s <- pos_e - group + 1

  d=rep(n_i,ncol(x))/nrow(x)
  
  l2.max=lambda2.max(x,y,ps,pe)$lambda2.max
  l2.min <- epsilon2*l2.max
  ss <- (log(l2.max)-log(l2.min))/(n.step-1)
  l2.seq <- numeric(n.step)
  n.l1 <- numeric(n.step)
  for ( i in 1:n.step )
  {
    l2.seq[i] <- exp(log(l2.max)-ss*(i-1))
    n.l1[i] <- 2*i-1
  }

  l1.max=lambda1.max(x,y,d,ps,pe,group,pos_s,pos_e,l2.max)$lambda1.max

  lambda1.seq <- numeric()
  lambda2.seq <- numeric()

  lambda1.seq <- c(lambda1.seq,l1.max)
  lambda2.seq <- c(lambda2.seq,l2.max)

  for ( i in 2:n.step )
  {
    l1.max=lambda1.max(x,y,d,ps,pe,group,pos_s,pos_e,l2.seq[i])$lambda1.max
    l1.min <- epsilon1*l1.max
    ss <- (log(l1.max)-log(l1.min))/(n.l1[i]-1)
    l1.seq <- numeric(n.l1[i])
    for ( j in 1:n.l1[i] )
    {
      l1.seq[j] <- exp(log(l1.max)-ss*(j-1))
       if ( i != n.step | j != n.l1[n.step] )
       {
         lambda1.seq <- c(lambda1.seq,l1.seq[j])
         lambda2.seq <- c(lambda2.seq,l2.seq[i])
       }
    }
  }

  cv = parCv(cl,lambda1.seq,lambda2.seq,m,x,y,x0,group,fold.d,fold.g,n.fold,gamma)

  list(lambda1.seq=lambda1.seq,lambda2.seq=lambda2.seq,PE=cv$PE)
}
