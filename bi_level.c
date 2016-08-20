#include "stdio.h"
#include "math.h"
#include "stdlib.h"

int Sum(int *G, int N)
{
  int i,value=0;
  for( i = 0; i < N; i++)
     value=value+G[i];
     
  return value;
}

double Norm(double *D1, double *D2,int P, int M)
{
  int i;
  double sum=0,value;
  
  for ( i = 0; i < P; i++)
  {
    sum = sum+D1[i]*D1[i];
  }
  
  for ( i = 0; i < M; i++)
  {
    sum = sum+D2[i]*D2[i];
  }
  
  value = pow(sum,0.5);
  return value;
}

double Inner_prod2(double *X, double *Z, double *Beta, double *Alpha, int P, int M)
{
  int j;
  double sum=0;
  
  for (j = 0; j < P; j++)
    sum = sum + X[j]*Beta[j];
  
  for (j = 0; j < M; j++)
  	 sum = sum + Z[j]*Alpha[j];

  return sum;
}

double Inner_prod(double *X, double *Beta, int P)
{
  int j;
  double sum=0;
  
  for (j = 0; j < P; j++)
    sum = sum + X[j]*Beta[j];

  return sum;
}

double Penalty(double Lambda, double Gamma, double Theta)
{
  double value;
  
  if ( Theta <= Lambda*Gamma)
    value = Lambda*Theta - pow(Theta,2)/2/Gamma;
  else 
    value = pow(Lambda,2)/2;
  return value;
}

double Penalty_der(double Lambda, double Gamma, double Theta)
{
  double value = 0;
  
  if ( Theta <= Lambda*Gamma)
    value = Lambda - Theta/Gamma;
  else 
    value = 0;
  return value;
}

double Lambda_approx(double *Tuning, double *Beta, int *Pos_s, int *Pos_e, int M, int J, int K)
{
  int k;
  double sum=0, value=0;
  double lambda1=Tuning[0],lambda2=Tuning[1],a=Tuning[2],b=Tuning[3];
  
  for ( k = Pos_s[J]; k <= Pos_e[J]; k++)
    sum = sum + Penalty(lambda2,b,fabs(Beta[k]));
    
  sum = sum*M/(Pos_e[J]-Pos_s[J]+1);
  
  value = Penalty_der(lambda1,a,sum)*Penalty_der(lambda2,b,fabs(Beta[Pos_s[J]+K]));
  
  return value;
}

void Update_alpha(double *Z, double *R, double *W, double *Alpha, int M, int N)
{
  int i, j;
  double sum1, sum2;
  
  for ( j = 0; j < M; j++)
  {
  	  sum1 = 0;
  	  sum2 = 0;
  	  for ( i = 0; i < N; i++)
  	  {
  	  	 sum1 = sum1 + Z[M*i+j]*W[i]*R[i];
  	  	 sum2 = sum2 + Z[M*i+j]*W[i]*Z[M*i+j];
  	  } 
  	  Alpha[j] = sum1/sum2 + Alpha[j];
  }
  
}

double S_func(double Z, double C)
{
  double value;
  if ( Z > C)
    value = Z - C;
  else if ( Z < -C)
    value = Z + C;
  else 
    value = 0;
    
  return value;	
}

void Update_beta(double *X, double *R, double *W, double *Beta, double Lambda, 
                 int N,int P,int *Pos_s, int J, int K)
{
  int i, col;
  double sum1=0, sum2=0,z;
  
  col = Pos_s[J]+K;
  for ( i = 0; i < N; i++)  
  {
  	 sum1 = sum1 + X[P*i+col]*W[i]*R[i];
  	 sum2 = sum2 + X[P*i+col]*W[i]*X[P*i+col];
  }
  z = (sum1+sum2*Beta[col])/N;
  
  Beta[col] = S_func(z,Lambda)*N/sum2;
}

double Bi_level_MCP_MCP(double *Y, double *Z, double *X, double *W, int *G, int *Param, double *Tuning, double *Epsilon,
              double *Alpha, double *Beta)
{
  int i,j,k,pos_s[Param[1]],pos_e[Param[1]],count=0,col;
  //d is the number of group, p is the number of covariates;
  int n=Param[0],d=Param[1],p=Param[2],m=Param[3],iter=Param[4];
    
  double **Xt,**Zt,r[Param[0]],old_beta[Param[2]], old_alpha[Param[3]];
  
  double lambda, diff, diff_beta[Param[2]], diff_alpha[Param[3]], y_t[Param[0]];

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

  Zt = malloc(sizeof(double *) * Param[0]); 
  for ( i = 0; i < Param[0]; i++)   
  { 
     Zt[i] = malloc(sizeof(double) * Param[3]); 
  }
    
  for ( i = 0; i < Param[0]; i++ )   
  {
    for ( j = 0; j < Param[3]; j++ )   
    {
       Zt[i][j] = *(Z+Param[3]*i+j); 
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
    y_t[i] = Inner_prod2(Xt[i],Zt[i],Beta,Alpha,p,m);
    r[i] = Y[i] - y_t[i];
  }

  do
  {
    for ( j = 0; j < p; j++)
      old_beta[j] = Beta[j];  
    
    for ( j = 0; j < m; j++)
      old_alpha[j] = Alpha[j];   
      
    Update_alpha(Z, r, W, Alpha, m, n);
    //update residuals;
    for ( i = 0; i < n; i++)
      r[i] = r[i] + Inner_prod(Zt[i],old_alpha,m) - Inner_prod(Zt[i],Alpha,m);
    
    for ( j = 0; j < d; j++)
    {
      for ( k = 0; k < G[j]; k++)
      {
        lambda = Lambda_approx(Tuning, Beta, pos_s, pos_e ,m, j, k);
        Update_beta(X, r, W, Beta, lambda, n, p, pos_s, j, k); 
        //update residuals;
        col = pos_s[j]+k;
        for ( i = 0; i < n; i++)
          r[i] = r[i] + Xt[i][col]*old_beta[col] - Xt[i][col]*Beta[col];
      }     
    }    

    for ( j = 0; j < p; j++)
      diff_beta[j] = old_beta[j]-Beta[j];
      
    for ( j = 0; j < m; j++)
      diff_alpha[j] = old_alpha[j] - Alpha[j];
      
    diff = Norm(diff_beta,diff_alpha,p,m);        
    count=count+1;
  }
  while ( count <= iter && diff > Epsilon[0]); 
  printf("iter=%d\n",count); 
  Update_alpha(Z, r, W, Alpha, m, n);
  
  for ( i = 0 ; i < n; i++)
  {  
    free(Xt[i]);
    free(Zt[i]);
  }  
  free(Xt);
  free(Zt);
  return 0;
}

/*main()
{
  double X[200]={-2.02792237503735,-0.330644698053631,-1.92366321719462,-0.632964974024835,-0.809818885517832,-0.415076360499092,-2.64693571515404,-0.0627154301610861,-0.719457117128918,-0.271977153111523,-0.527705744654312,0.00915369739644391,-0.868843360382762,0.545132266966089,0.0867123067718595,0.0469080776837625,-1.53861149340959,-0.0483574941335827,0.291039430475395,-0.229123202960524,
-0.978698914985543,-0.330644698053631,0.62897333994519,-0.632964974024835,-0.185500061676002,-0.415076360499092,-0.60498961197707,-0.0627154301610861,-2.2664161590977,-0.271977153111523,-1.70026503461857,0.00915369739644391,-0.400959509551295,0.545132266966089,-1.76403820768562,0.0469080776837625,-1.45773445160678,-0.0483574941335827,-2.41736012629236,-0.229123202960524,
-0.86732205063256,-0.330644698053631,-2.05501598466391,-0.632964974024835,1.05098620238303,-0.415076360499092,1.12126574087383,-0.0627154301610861,0.721695269927193,-0.271977153111523,0.119458030283411,0.00915369739644391,1.18668837946454,0.545132266966089,2.09020698954258,0.0469080776837625,1.372253697752,-0.0483574941335827,1.57953002522718,-0.229123202960524,
1.76662937663,-0.330644698053631,0.222712512905671,-0.632964974024835,-1.75208922387492,-0.415076360499092,1.0794379528628,-0.0627154301610861,1.90110621736895,-0.271977153111523,1.32877564148231,0.00915369739644391,2.13511380872574,0.545132266966089,0.533934163059534,0.0469080776837625,1.4077013764245,-0.0483574941335827,-0.962738039730843,-0.229123202960524,
0.864641894940799,-0.330644698053631,0.556984195061232,-0.632964974024835,2.2524563496748,-0.415076360499092,0.30786204255213,-0.0627154301610861,-0.336056808868801,-0.271977153111523,-1.90092598687977,0.00915369739644391,-1.76080773893053,0.545132266966089,-1.47347390877732,0.0469080776837625,-1.13125992976465,-0.0483574941335827,0.726023702441664,-0.229123202960524,
0.24853441381693,0.608242600431901,0.514001830789287,1.32746643386832,-0.111206876197815,0.092980511080467,0.14867191816847,0.458082280937348,0.139825719559856,1.37667353806605,0.536132618877385,-0.00380131680436065,-0.0582383158651378,0.661784522643246,0.105331731417793,1.8523954708542,0.269530160120906,-0.842615001034159,0.156701001575793,0.851776098937761,
0.24853441381693,2.63538871318355,0.514001830789287,1.48027370431487,-0.111206876197815,1.66371751133851,0.14867191816847,2.27006945125226,0.139825719559856,2.18023058824919,0.536132618877385,1.52406474725898,-0.0582383158651378,1.11634669637674,0.105331731417793,1.40216058252916,0.269530160120906,2.81200143346513,0.156701001575793,1.21381003954883,
0.24853441381693,0.367934072334627,0.514001830789287,-1.51115503178626,-0.111206876197815,-0.802907965503185,0.14867191816847,-2.13742379749575,0.139825719559856,0.0988196955036095,0.536132618877385,-0.892554101628161,-0.0582383158651378,-1.68895308487999,0.105331731417793,-1.29949292852177,0.269530160120906,-0.745549143952234,0.156701001575793,1.39788120432411,
0.24853441381693,-0.77276029415103,0.514001830789287,1.01970579798543,-0.111206876197815,2.15583686267916,0.14867191816847,-0.0668383168891106,0.139825719559856,-1.5579831512077,0.536132618877385,-2.16065697579368,-0.0582383158651378,-1.49613899867882,0.105331731417793,-1.59764926714281,0.269530160120906,-0.899008495369878,0.156701001575793,-2.36274587560487,
0.24853441381693,-1.1855816015309,0.514001830789287,0.848533965741811,-0.111206876197815,-1.03424511709949,0.14867191816847,-0.210312466999323,0.139825719559856,-0.73785490505353,0.536132618877385,1.487179159985,-0.0582383158651378,-1.31870047029162,0.105331731417793,-0.591954246137595,0.269530160120906,-0.0830413224409467,0.156701001575793,0.0448945475967874},        
         Y[10]={-8.11776845685042,-0.69855405639037,0.572924710871367,-0.464404582462445,3.58222388781311,-0.537564109454422,
          3.28491501683857,-2.04947236162492,-2.63460735780985,2.21836997418644},
         Z[20]={1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1},
 Beta[20],Alpha[2]={0,0};
  int G[10]={2,2,2,2,2,2,2,2,2,2}, Param[5]={10,10,20,2,1000};
  double Tuning[4]={1,1,1000,1000},Epsilon[1]={1E-10}, W[10]={1,1,1,1,1,1,1,1,1,1};
  int j;
  for ( j = 0; j < 20; j++)
     Beta[j] = 0;
  Bi_level_MCP_MCP(Y, Z, X, W, G, Param, Tuning, Epsilon, Alpha, Beta);
  for ( j = 0; j < 20; j++)
     printf("beta[%d]=%4.10f\n",j, Beta[j]);
  for ( j = 0; j < 2; j++)
     printf("alpha[%d]=%4.10f\n",j, Alpha[j]);
}*/
