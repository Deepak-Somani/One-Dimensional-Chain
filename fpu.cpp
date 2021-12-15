// Fermi-Pasta Ulam Chain
// By Deepak Somani
// Runge-Kutta 4th order application

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>

#define  N   128
 
double t1=1.1,t2=0.9;
#define seed            1
#define MAX   			375000
#define MAX2            125000
#define e  2.718
#define r               0.05
#define bins            20


using namespace std;

void runge4(double t, double x[], double dt);
void KE(double x[],double T[]);
void col(double x[],double J[]);
double f(double t, double x[], int i);
double K=1, M=1,d1,d2,a=0.1 ,c=0.0000000005;
unsigned long long int s=0,ncol=0;
double dt = 0.0005;
double bin[bins],binc[bins];
double energy;
double xm[N],xd[N];
double rel[N+1];


int main(){srand(0);
    double x[2*N+2];

    fstream file,op,fl,kine,flux,diff,dens;
	dens.open("bins.dat",ios::out);
    kine.open("kine.dat",ios::out);
    fl.open("ncol.dat",ios::out);
    flux.open("flux.dat",ios::out);
    diff.open("diff.dat",ios::out);
    op.open("output.dat",ios::out);
      for(int i=0;i<N;i++){
        x[i]=i;
        x[i+N]=(rand()%100)/100-0.5;
    }
    x[2*N+1]=x[2*N]=0;
    double T[N];
    double J[N]; //non collision part of heat flow
    double J2[N]; //collision part of heat flow
	double t = 0;
	int j;
	double binl = ((N+1)*1.0)/bins;
	for(int i = 0;i<bins; i++){
		bin[i]=0;
	}		
    for(int i=0;i<N;i++){
            T[i]=0;
            J[i]=0;
            J2[i]=0;
			xm[i]=0;
			xd[i]=0;
    }
	for (j=0; t<=MAX2 ;j++)	{
		t += dt;
		runge4(t,x,dt);
		for(int i =0;i<N;i++){
			xm[i]+=x[i]*dt;
		}
	}
	for(int i =0;i<N;i++){
		xm[i]/=t;
	}
	
	t=0;
    
 
    for (j=0; t<=MAX ;j++)	{
		t += dt;
		runge4(t,x,dt);
        col(x,J2);
		KE(x,T);
		for(int i=0;i<bins;i++)binc[i]=0;
		for(int i=0;i<N;i++){     //j[i] is the heat flow between i-1 th particle and i th particle
			xd[i]=(x[i]-xm[i])*(x[i]-xm[i])*dt;
			int b=(x[i]+1)/binl;
		    binc[b]++;
		}
		for(int i=1;i<N;i++){     //j[i] is the heat flow between i-1 th particle and i th particle
            J[i]+=-K*x[i+N]*((x[i]-x[i-1]-1)+a*(x[i]-x[i-1]-1)*(x[i]-x[i-1]-1)*(x[i]-x[i-1]-1))*dt;
		}
		for(int i=0;i<bins;i++)bin[i]+=binc[i]*dt;

		if(j%MAX==0){
            op<<t<<"  ";
            flux<<t<<" ";
            kine<<t<<"  ";
			fl<<t<<"   ";
			diff<<t<<"   ";
			dens<<t<<"   ";
            for(int i=0;i<N;i++){
                    kine<<T[i]<<" ";
                    op<<x[i]<<"  ";
                    flux<<(J[i]+J2[i])/t<<" ";
					diff<<xd[i]/t<<"  ";
            }
			for(int i=0;i<bins;i++)dens<<bin[i]/t<<"  ";
			dens<<endl;
			diff<<endl;
			fl<<ncol<<endl;
            op<<endl;
            kine<<endl;
            flux<<endl;

		}
	}
	file.open("steadystate.dat",ios::out);
    for(int i=0;i<2*N+2;i++){
        file<<x[i]<<endl;
    }
	
	return(0);
}

void runge4(double t, double x[], double dt){
    int i=0;
    double k[4][2*N+2];
    double w[3][2*N+2];
    for(i=0;i<2*N+2;i++)w[0][i]=x[i]+0.5*(k[0][i]=dt*f(t,x,i));
    for(i=0;i<2*N+2;i++)w[1][i]=x[i]+0.5*(k[1][i]=dt*f(t+dt/2,w[0],i));
    for(i=0;i<2*N+2;i++)w[2][i]=x[i]+(k[2][i]=dt*f(t+dt/2,w[1],i));
    for(i=0;i<2*N+2;i++)k[3][i]=dt*f(t+dt,w[2],i);

    for (i=0;i<2*N+2;i++) x[i]+=(k[0][i]+2.0*k[1][i]+2.0*k[2][i]+k[3][i])/6.0;

}

double f(double t, double x[], int i){
    if(i<N)return(x[i+N]);
    if(i==N){
		d1=x[i-N]+1;
	    d2=x[i-N+1]-x[i-N];
        return -K*(d1-d2)/M-K*a*((d1-1)*(d1-1)*(d1-1)-(d2-1)*(d2-1)*(d2-1))/M-x[2*N]*x[N];


    }
  if(i==2*N-1){
	    d1=x[i-N]-x[i-N-1];
	    d2=N-x[i-N];
        return -K*(d1-d2)/M-K*a*((d1-1)*(d1-1)*(d1-1)-(d2-1)*(d2-1)*(d2-1))/M-x[2*N+1]*x[2*N-1];


  }
    if(i==2*N){

        return (x[N]*x[N]-t1);
    }
    if(i==2*N+1){

        return (x[2*N-1]*x[2*N-1]-t2);
    }

	d1=x[i-N]-x[i-N-1];
	d2=x[i-N+1]-x[i-N];
    return -K*(d1-d2)/M-K*a*((d1-1)*(d1-1)*(d1-1)-(d2-1)*(d2-1)*(d2-1))/M ;
}

void KE(double x[],double T[]){

    for(int i=N;i<2*N;i++){

        T[i-N]=T[i-N]*(s*1.0/(s*1.0+1.0))+(M*x[i]*x[i])/(s*1.0+1.0);

    }
    s++;
}

void col(double x[],double J[]){
    double nrel[N+1];
	for(int i=0;i<N-1;i++)
	{
		nrel[i+1]=x[i+1+N]-x[i+N];
	}
	nrel[0]=x[N];
	nrel[N]=-x[2*N-1];
	
	for(int i=0;i<N+1;i++){
		if(nrel[i]*rel[i]<0)ncol++;
		rel[i]=nrel[i];
	}
	
}