#include "stdio.h"
#include "math.h"
#include "stdlib.h"

double Sgn(double A)
{
  double value;
  if ( fabs(A) < 1E-12) 
     value = 0;
  if ( A > 0)
     value = 1;
  else if ( A < 0)
     value = -1;
  return value;
}

int Sum(int *G, int n)
{
  int i,value=0;
  for( i = 0; i < n; i++)
     value=value+G[i];
     
  return value;
}

double Norm(double *z,int n)
{
  int i;
  double sum=0,value;
  for ( i = 0; i < n; i++)
  {
    sum += z[i]*z[i];
  }
  value = pow(sum,0.5);
  return value;
}

double Inner_prod(double *X, double *Beta, int P)
{
  int j;
  double sum=0;
  
  for (j = 0; j < P; j++)
    sum = sum + X[j]*Beta[j];

  return sum;
}

double S1_func(double Z , double Lambda)  
{
  double value; 
  if (Z > Lambda) 
      value = Z - Lambda;
  else if (Z < -Lambda) 
      value = Z + Lambda;
  else 
      value = 0;
  return value ;
}

double G_func(double *Beta, double Lambda1, double Lambda2,
              double A, double B, int *G, int *Pos_s, int *Pos_e, int J)
{
  double g=0,sum=0,norm,inv=0;
  int j,s=Pos_s[J],e=Pos_e[J];

  for ( j = s ; j <= e; j++)
  {
    sum=sum+pow(Beta[j],2);
  }
  norm= pow(sum,0.5);

  if ( norm <= 1E-4 )
    inv=10000;
  else if ( norm > 1E-4 )
    inv=1/norm;

  if ( norm <= A*Lambda1*pow(G[J],0.5) )
    g=1+inv*(Lambda1*pow(G[J],0.5)-norm/A);
  else if ( norm > A*Lambda1*pow(G[J],0.5) )
    g=1;

  return g;
}

void Inner_loop(double *Z, double *Beta, double Lambda1, double Lambda2,
                double A, double B, double *D, 
                double Gg, int N, int P, int *G, int *Pos_s, int *Pos_e, int J)
{
  int j,s=Pos_s[J],e=Pos_e[J];
  double sc[G[J]], s_norm, criterion,cons;
  
  for ( j = s;j <= e;j++)
  {
    //let i-pos_s be the position for z[j]; 
    if ( fabs(Z[j-s]) > Lambda2*B*Gg*D[j] )
       sc[j-s]=Z[j-s]/D[j];
    else if ( fabs(Z[j-s]) <= Lambda2*B*Gg*D[j] )
       sc[j-s] = S1_func(Z[j-s],Lambda2)/(1-1/B/Gg/D[j])/D[j];
  }

  s_norm = Norm(sc,G[J]);     
  criterion = Lambda1*A*pow(G[J],0.5);
  
  if ( s_norm >= criterion)
  {
    for ( j = s; j <= e; j++ )
    {
      Beta[j] = sc[j-s];
    }
  }
  else 
  {
    cons = 1 - pow(G[J],0.5)*Lambda1/s_norm;
    if ( cons <= 0 )
      cons = 0;
  	  
    for ( j = s; j <= e; j++ )  
    {
      Beta[j] = cons*sc[j-s]*A/(A-1);
    }
  } 
}

void Update_beta(double *X, double *D, double *R, double *Beta, double Lambda1, double Lambda2,
                 double A, double B,
                 int *G, int N, int P, int *Pos_s, int *Pos_e, int J)
{
  int i,j,s=Pos_s[J],e=Pos_e[J],count=0;
  double Gg, diff[G[J]],dif, sum,old_old_beta[G[J]],old_beta[G[J]],z[G[J]];//,rk[N]

  for ( j = s; j <= e; j++)
    old_old_beta[j-s] = Beta[j];  

 /* for ( i = 0; i < N; i++)
  {
  	 rk[i]=0;
  	 for ( j = s; j <= e; j++)
  	 {
  	 	rk[i]=X[P*i+j]*Beta[j];
  	 }
  	 rk[i]=R[i]+rk[i];
  }*/

  /*for ( j = s;j <= e;j++)
  {
    //let i-pos_s be the position for z[j];
    z[j-s] = 0;
    for ( i = 0; i < N; i++)
    {
      //z[j-s] = z[j-s] + X[P*i+j]*rk[i]/N;
      z[j-s]=z[j-s]+X[P*i+j]*(R[i]+X[P*i+j]*Beta[j])/N;
    } 
  }*/
  
  for ( j = s;j <= e;j++)
  {
    //let i-pos_s be the position for z[j];
    z[j-s] = 0;
    for ( i = 0; i < N; i++)
    {
      z[j-s] = z[j-s] + X[P*i+j]*R[i]/N;
    }
    z[j-s] = z[j-s] + D[j]*Beta[j];
  }
  
  do
  {
    for ( j = s; j <= e; j++)
      old_beta[j-s] = Beta[j];  

    Gg=G_func(Beta, Lambda1,Lambda2, A, B, G, Pos_s, Pos_e, J);
    Inner_loop(z, Beta, Lambda1,Lambda2,A,B, D, Gg, N, P, G, Pos_s, Pos_e, J);
    
    for ( j = s; j <= e; j++)
      diff[j-s] = old_beta[j-s] - Beta[j];  
   
    dif =  Norm(diff,G[J]);  
    count = count+1;
  }
  while ( count <= 100 && dif > 1E-10); 
  //printf("iter=%d\n",count); 
  //printf("iter=%e\n",dif);

  //update r here;
  for ( i = 0; i < N; i++)
  {
    sum = 0;
    for ( j = s; j <= e; j++)
    {
      sum = sum + X[P*i+j]*(Beta[j] - old_old_beta[j-s]);
    }
    R[i] = R[i] - sum;
  }
}

double sparse_MCP(double *Y, double *X, int *G, double *D, 
                  int *Param, double *Tuning, double *Epsilon, double *Beta,double *Diff)
{
  int i,j,pos_s[Param[1]],pos_e[Param[1]],count=0;
  //d is the number of group, p is the number of covariates, Tuning[4]=(lambda1,lambda2,a,b);
  int n=Param[0],d=Param[1],p=Param[2],iter=Param[3];
    
  double **Xt,r[Param[0]],old_beta[Param[2]],lambda1=Tuning[0],lambda2=Tuning[1],a=Tuning[2],b=Tuning[3];
  
  double diff, diff_beta[Param[2]], y_t[Param[0]];

  Xt = malloc(sizeof(double *) * Param[0]); 
  for ( i = 0; i < Param[0]; i++)   
  { 
     Xt[i] = malloc(sizeof(double) * Param[2]); 
  }
    
  for ( i = 0; i < Param[0]; i++ )   
  {
    for ( j = 0; j < Param[2]; j++ )   
    {
       Xt[i][j] = *(X+Param[2]*i+j); 
    }
  }   

  for( j = 0; j < d; j++)
  {
    if ( j == 0)
       pos_s[j]=0;
    else
       pos_s[j]=Sum(G,j);
    pos_e[j]=pos_s[j]+G[j]-1;
  }
  
  // do initialization in outer loop;
  for( i = 0; i < n; i++)
  {
    y_t[i] = Inner_prod(Xt[i],Beta,p);
    r[i] = Y[i] - y_t[i];
  }

  do
  {
    for ( j = 0; j < p; j++)
      old_beta[j] = Beta[j];  
        
    for ( j = 0; j < d; j++)
    {
      Update_beta(X, D, r, Beta, lambda1,lambda2,a,b, G, n, p, pos_s, pos_e, j); 
    }    

    for ( j = 0; j < p; j++)
      diff_beta[j] = old_beta[j]-Beta[j];
          
    diff = Norm(diff_beta,p);     
    count=count+1;
  }
  while ( count <= iter && diff > Epsilon[0]); 
  printf("iter=%d\n",count); 
  
  Diff[0]=diff;
  for ( i = 0 ; i < n; i++)
  {  
    free(Xt[i]);
  }  
  free(Xt);
  return 0;
}

