#include<iostream>
#include<conio.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
using namespace std;
doule fun(double x);
double funct(double x,double a, double z);
double Trapez(int n, double a, double b, double (*func)(double));
int index_ik(int i,int k);
double romberg(int *n0, double a, double b, double (*func)(double),double eps);
float * conv_a(float limit, float upper_limit,int points, float (*func)(float));
int mmax=100;
doule fun(double x){return x;
}
//Making of the function which needs to be integrated
double funct(double x,double a, double z)
{ 
  return (exp(-(pow(x,2)/pow(a,2))))/(a*sqrt(3.14*((x*x)+(z*z))));
}

//function to calculate the trepez integral of the given function
double Trapez(int n, double a, double b, double (*func)(double))
{
	double h=(b-a)/n;
	double sum=(func(a)+func(b))/2;
	for(int i=1;i<n-1;i++)
	{sum+=func(a+(h*i));
	}
	return h*sum;
 
}

/* Funktion, die eindeutigen Index definiert von ik -> index_ik  fuer  i>=k und i,k<=mmax */ 
int index_ik(int i,int k)
   { return (i-k + k*(mmax+1)-k*(k-1)/2); }


//actual romberg integral function to run the process of integration until accuracy ac is acheived
double romberg(int *n0, double a, double b, double (*func)(double),double eps)
{ int k,m,n;  /* fuer Indizes */ 
  double *h;  /* Schrittweiten fuer j=0,...,mmax */  
  double *Tsum; /* Trapezsummen fuer j=0,...,mmax */  
  double *tildeT; /* Neville-Schema tilde T_{jk}  bei h=0 */
  double result;
 
  h=(double *)malloc((mmax+1)*sizeof(double));       /* Speicher fuer h in Schritt m */
  Tsum=(double *)malloc((mmax+1)*sizeof(double));    /* und T fuer diese h */
  tildeT=(double *)malloc(((mmax+1)*(mmax+2))/2*sizeof(double));  /* Speicher fuer Neville Schema */ 

  h[0]=(b-a)/(double)(*n0-1);
  n=*n0;
  Tsum[0]=Trapez(n,a,b,func);
  tildeT[index_ik(0,0)]=Tsum[0];

  for(m=1;m<=mmax;m++) 
    { h[m]=h[m-1]/2;         /* Trapezsumme fuer halbiertes h */ 
      n=2*n;
      Tsum[m]=Trapez(n,a,b,func);

      tildeT[index_ik(m,0)]=Tsum[m]; /* fuer  i=m und k=0 */ 

      /* generate tildeT i=m k=1,...,m */ 
      for(k=1;k<=m;k++)
	{ tildeT[index_ik(m,k)]=-h[m-k]*h[m-k]/(h[m]*h[m]-h[m-k]*h[m-k])*tildeT[index_ik(m,k-1)]
	                  + h[m]  *h[m]  /(h[m]*h[m]-h[m-k]*h[m-k])*tildeT[index_ik(m-1,k-1)];
         }

      printf("%5d   %15.6e  %15.6e \n",m,tildeT[index_ik(m,m)],Tsum[m]);  /* just to observe convergence */
      if(fabs(tildeT[index_ik(m,m)]-tildeT[index_ik(m-1,m-1)])<=eps) break; /* stop when accuracy reached */

    }
  *n0=n;

  if(m > mmax) {
    printf("Konvergenz nicht erreicht innerhalb mmax = %d Schritte!\n",mmax);
  }
  result=tildeT[index_ik((int) fmin(m,mmax),(int) fmin(m,mmax))];

  free(tildeT);
  free(Tsum);
  free(h);

  return result;

}
//function to make a convergent series for a
float * conv_a(float limit, float upper_limit,int points, float (*func)(float))
{   
   float interval=(upper_limit-limit)/points;
   float * series;
   int index=0;
   for(int i=limit+interval;i<=upper_limit;i+=interval)  //started from limit+interval to avoid the cases where the approahing limit value is undefined
   {
   	 series[index]=func(i); 
   }	
return series;
}


//the Main function
int main()
{ char choice;
  double a=2,b=5;
  int n=3;
  do
  {
  	cin>>limit>>upper_limit>>points;
  float * series=conv(limit,upper_limit,points,fun);
  
  
  
  cout<<"do you want to try once more because apparently you have nothing to do? press(y/n) ";
  
  cin>>choice;
  }while(choice=='y') ;
  return 0;
}
