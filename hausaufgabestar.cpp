#include<iostream>
#include<conio.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#define pi 3.1428
#include<vector>


using namespace std;


double funct(double x[]);//actual integral to be calculated
double Trapez(int n, double *paras, double (*func)(double *));//Trapez integral
int index_ik(int i,int k);//Index function for Tilde T
double romberg(int *n0, double * paras, double (*func)(double*),double eps);//Romberg integral method
vector<double> conv_type(double min, double max,int length, double func_type(double));//function to produce a convergent series for a
int mmax=20;//maximum number of loops allowed in romberg integral
double new_funct(double x[]);//convergent function type function example


//convergence function
vector<double> conv_type(double min, double max,int length, double func_type(double)) 
{
    
    
    double interval=(max-min)/length;
	double value=min;
    vector<double> series(length);
    for (int i=0;i<length;i++){series[i]=func_type(value);value+=interval;}
    

    return series;
}

//trial function to ensure that the program is working
double new_funct(double x)
{  return x;
}


//Making of the function which needs to be integrated
// x[] is the array of the values for x,a,z in the specified order
double funct(double x[])
{ 
	
   return ((exp(-pow((x[0]/x[1]),2)))/(x[1]*sqrt(pi)*sqrt((x[0]*x[0])+(x[2]*x[2]))));
  
}

//function to calculate the trepez integral of the given function
//paras is the pointer containing the values paras[0]=start point,  paras[1]=end point,   paras[2]=a, paras[3]=z
double Trapez(int n, double paras[], double (*func)(double []))
{
	double h=((paras[1])-(paras[0]))/n;
	double temp[3]; //to calculate the function at the start point(paras[0]=a)
	temp[0]=paras[0];
	temp[1]=paras[2];
	temp[2]=paras[3];
	double temp_b[3];//to calculate the function at the end point(paras[1]=b)
	temp_b[0]=paras[1];
	temp_b[1]=paras[2];
	temp_b[2]=paras[3];
	
	
	double sum=(func(temp)+func(temp_b))/2;
	double denominator = sum;
	for(int i=1;i<n-1;i++)	
	{
	temp[0]=(paras[0])+(h*i);
	sum+=func(temp);
	}
	
	return h*sum;
 
}

//Funktion, die eindeutigen Index definiert von ik -> index_ik  fuer  i>=k und i,k<=mmax  
int index_ik(int i,int k)
   { return (i-k + k*(mmax+1)-k*(k-1)/2); }


//actual romberg integral function to run the process of integration until accuracy ac is acheived*/
double romberg(int *n0, double  paras[], double (*func)(double*),double eps)
{ int k,m,n;  /* fuer Indizes */ 
  double *h;  /* Schrittweiten fuer j=0,...,mmax */  
  double *Tsum; /* Trapezsummen fuer j=0,...,mmax */  
  double *tildeT; /* Neville-Schema tilde T_{jk}  bei h=0 */
  double result;
 
  h=(double *)malloc((mmax+1)*sizeof(double));       /* Speicher fuer h in Schritt m */
  Tsum=(double *)malloc((mmax+1)*sizeof(double));    /* und T fuer diese h */
  tildeT=(double *)malloc(((mmax+1)*(mmax+2))/2*sizeof(double));  /* Speicher fuer Neville Schema */ 

  h[0]=((paras[1])-(paras[0]))/(double)(*n0-1);
  n=*n0;
  Tsum[0]=Trapez(n,paras,func);
  tildeT[index_ik(0,0)]=Tsum[0];

  for(m=1;m<=mmax;m++) 
    { h[m]=h[m-1]/2;         /* Trapezsumme fuer halbiertes h */ 
      n=2*n;
      Tsum[m]=Trapez(n,paras,func);

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


//the Main function
int main()
{ 
  double arguments[4]; //start point, end point , a, z , starting n to be given by user in array 'arguments'
  int n;
  int length;//the number of data points to be evaluated i.e how many different values for a(right now its not running for more than 9)
  double results[length];//the array of results of evaluated potentials
  cout<<"how many data points for respective a's do you want your integral to be evaluated??\n and how many stuetzstellen do you want to start with? \nwhat is the z coordinate";
  cin>>length>>n>>arguments[3];
  vector<double> conv_parameter(length);//vector to hold the return vector of the convergent series function called in next line
  conv_parameter=conv_type(0,4,length,&new_funct);//calling the function to produce a convergent series
 
//for loop to calculate thr integral
  for(int k=0;k<length-1;k++)
      {arguments[1]=conv_parameter.at(k+1);
       arguments[0]=-arguments[1];
       arguments[2]=arguments[1];
       results[k]=romberg(&n,arguments,&funct,0.0003);
      }
 
//for loop to show the result  
 for(int i=0;i<length-1;i++){cout<<"a = "<<(conv_parameter.at(i+1))<<" potential =   "<<results[i]<<endl;}
 
 return 0;

}
