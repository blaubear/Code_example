#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 #include <iostream>
 


int  main(void){
    int i, j,N,k,n,qwer;
    double R,h,f0,sum;
    double c=4;
    N=7;
//	for(N=0;N<20;N++){
    sum=0;
    qwer=0;
    k=0;
    for(n=0;n<=N;n++){
    for(i=0;i<=n;i++){
    	k=k+1;
    	sum=sum+abs(n-i);
	}
}
c=sum/k;
if(c==7/3){
	R=c;
	qwer=N;
}
 //}
printf("%lf",R);

}




