#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

int main(void){
	
   int     i, j, n, m, n_2, m_2;
   double  h, L, x, a, T, tau, max, h_2, tau_2,sum,max_2,sum_2;
   x   = 0;
   a   = 0.01;
   L   = 1;
   
   n   = 50;
   m   = 100;
   
   n_2 = n;
   m_2 = 50*m;
   h   = L/(n-1);
   h_2 = L/(n_2-1);
   T=1;
   tau = T/(m-1);
   tau_2 = T/(m_2-1);
   //printf ("%f", tau);
   double *u,*u_2;
   double *u_befor,*w;
   double *u_befor_befor;
   
   if (!(u = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
       
   if (!(u_2 = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    if (!(u_befor = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    if (!(u_befor_befor = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    if (!(w = (double*)malloc(n*sizeof(double)))) {
           printf("\n\nmem error\n\n");
           return (-1);
    }
    
    u_befor[0]=sin(tau);
    u_befor_befor[0]=0;
   
    for(j=0;j<n;j++){
    x=j*h;
    u[j]=0;
    //printf ("%f", u[j]);
	u_befor[j]=0;
    u_befor_befor[j]=0;
    w[j]=0;
	}
	

	
	for (j=2;j<m;j++){
		u[0]=sin(tau*j);
		u[n-1]=0;
		w[0]=2*a*pow((u_befor[1]-u_befor[0])/h,2);
    for (i=1;i<n-1;i++){
    	
    	x=i*h;
        w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
    	u[i]=(pow(tau,2)/(pow(h,2)))*((w[i]-w[i-1])*(u_befor[i+1]-u_befor[i])+w[i]*(u_befor[i+1]-2*u_befor[i]+u_befor[i-1]))+2*u_befor[i]-u_befor_befor[i];
    	//w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
   
    
    }
    for (i=0;i<n;i++){
    	u_befor_befor[i]=u_befor[i];
       	u_befor[i]=u[i];
    }
    
    //printf ("\n\n\n");
   for (i=0;i<n;i++){
   //	printf ("%f\n", u[i]);
   }
    
}
    u_befor[0]=sin(tau_2);
    u_befor_befor[0]=0;
   
    for(j=0;j<n_2;j++){
    x=j*h_2;
    u_2[j]=0;
    //printf ("%f", u[j]);
	u_befor[j]=0;
    u_befor_befor[j]=0;
    w[j]=0;
	}



	for (j=2;j<m_2;j++){
		u_2[0]=sin(tau_2*j);
		u_2[n_2-1]=0;
		w[0]=2*a*pow((u_befor[1]-u_befor[0])/h_2,2);
    for (i=1;i<n-1;i++){
    	
    	x=i*h_2;
        w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h_2,2);
    	u_2[i]=(pow(tau_2,2)/(pow(h_2,2)))*((w[i]-w[i-1])*(u_befor[i+1]-u_befor[i])+w[i]*(u_befor[i+1]-2*u_befor[i]+u_befor[i-1]))+2*u_befor[i]-u_befor_befor[i];
    	//w[i]=4*x*(1-x)+2*a*pow((u_befor[i]-u_befor[i-1])/h,2);
    }
    for (i=0;i<n_2;i++){
    	u_befor_befor[i]=u_befor[i];
       	u_befor[i]=u_2[i];
    }
    
  //  printf ("\n\n\n");
    for (i=0;i<n_2;i++){
   	//printf ("%f\n", u_2[i]);
   }
    
}

     printf ("\n\n\n");
   for (i=0;i<n;i++){
   	printf ("%f\n", u[i]);
   }
   
       printf ("\n\n\n");
   for (i=0;i<n;i++){
   	printf ("%f\n", u_2[i]);
   }
   
   
    max=fabs(u[0]-u_2[0]);
    for(i=0;i<n;i++){
   	if (fabs(u[i]-u_2[i])>=max)
   	max=fabs(u[i]-u_2[i]);
   }
   
    printf ("\n\n\n");
    printf ("%f\n", max);
    
    sum=0;
    for(i=0;i<n;i++){
    sum=sum+fabs(u[i]-u_2[i]);
    }
    
    printf ("\n\n\n");
    printf ("%f\n", sum);
    
    max_2=fabs(u_2[1]);
    for(i=0;i<n;i++){
   	if (fabs(u_2[i])>=max_2)
   	max_2=fabs(u_2[i]);
   }
   
   max=max/max_2;
   
   printf ("\n\n\n");
   printf ("%f\n", max);
   
    sum_2=0;
    for(i=0;i<n;i++){
    sum_2=sum_2+fabs(u_2[i]);
    }
    
    sum=sum/sum_2;
    
    printf ("\n\n\n");
    printf ("%f\n", sum);
   
   
   
   
  return 1;
}
