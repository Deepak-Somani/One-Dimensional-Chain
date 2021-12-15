//=====================One-Dimension Non-Linear Chain=======================
// Runge-Kutta Fourth Order Application
// By Deepak Somani

#include<stdio.h>
#include<math.h>
#define ksp 10 // spring constant
#define m 1 // mass
#define h 0.05 // time step
#define alpha 1 // quadratic term coefficient
#define beta 0.5 // cubic term coefficient

// h is considered to be the step difference of time 'dt'.
// The equation for non linear one dimension system is taken as:
// dv[n]/dt = -(k/m)*[(2*x[n] - x[n+1] + x[n-1] - alpha*({x[n+1]-x[n]}^2-{x[n]-x[n-1]}^2) - beta*({x[n+1]-x[n]}^3-{x[n]-x[n-1]}^3)] - 1
// k is spring constant and m is mass of particle which is taken identical for all.
// and dx/dt = v. - 2
// (ksp) is taken as 10.
// (m) is taken as 1. 
// (alpha) is taken as 1. 

float velo(float x0, float x1, float x2){
	return -h*(ksp/m)*((2*x1-x0-x2)-alpha*(pow((x2-x1),2)-pow((x1-x0),2)) -beta*(pow((x2-x1),3)-pow((x1-x0),3)));
	// step-1 is used for runge kutta step.
}

float dist(float v){
	return h*v;
	// step-2 is used for runge kutta step.
}

int main(){
	int n;
	printf("Enter the value of the dimension of chain:\n");
	scanf("%d",&n);
	float x[99][n+2],v[99][n+2],k[4][n+2],l[4][n+2],mx[99],mv[99];
	// l is step function for displacement and k is step function for velocity.
	int i,j;
	//the very first and last points are considered to be hinged to the ground, so x=0 and v=0 for them.
		for(j=0;j<99;j++){
		x[j][0]=0;
		x[j][n+1]=0;
		v[j][0]=0;
		v[j][n+1]=0;
	}
	for(j=0;j<4;j++){
		l[j][0]=0;
		k[j][0]=0;
		l[j][n+1]=0;
		k[j][n+1]=0;
	}
	
	
	// Initiating the values of displacement and velocities of different particles...
	
	printf("Initialise the values of displacement and velocity of particles\n");
	
	for(i=1;i<=n;i++){
		scanf("%f",&x[0][i]);
		scanf("%f",&v[0][i]);
	}
	
	printf("\n");
	
	for(j=0;j<98;j++){
		for(i=1;i<=n;i++){
			l[0][i]=dist(v[j][i]);
		}
	
		for(i=1;i<=n;i++){
			k[0][i]=velo(x[j][i-1],x[j][i],x[j][i+1]);
		}
	
		for(i=1;i<=n;i++){
			l[1][i]=dist(v[j][i]+k[0][i]/2);
		}
	
		for(i=1;i<=n;i++){
			k[1][i]=velo(x[j][i-1]+(l[0][i-1]/2), x[j][i]+(l[0][i]/2), x[j][i+1]+(l[0][i+1]/2) );
		}
	
		for(i=1;i<=n;i++){
			l[2][i]=dist(v[j][i]+k[1][i]/2);
		}
	
		for(i=1;i<=n;i++){
			k[2][i]=velo(x[j][i-1]+(l[1][i-1]/2), x[j][i]+(l[1][i]/2), x[j][i+1]+(l[1][i+1]/2) );
		}
	
		for(i=1;i<=n;i++){
			l[3][i]=dist(v[j][i]+k[2][i]);
		}
	
		for(i=1;i<=n;i++){
			k[3][i]=velo(x[j][i-1]+(l[2][i-1]), x[j][i]+(l[2][i]), x[j][i+1]+(l[2][i+1]) );
		}
		
		for(i=1;i<=n;i++){
			x[j+1][i]=x[j][i]+((l[0][i]+2*l[1][i]+2*l[2][i]+l[3][i])/6);
			v[j+1][i]=v[j][i]+((k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6);
			//at t= t+h every time the displacement and velocity of each particle gets updated.
		}
	}
	
	
	printf("\n");
	printf("The displacement of particles from their equlibrium position are as given for every time step taken:\n");
	for(j=0;j<99;j++){
		for(i=1;i<=n;i++){
			printf("%f\t",x[j][i]);
		}
		printf("\n");
	}
	printf("\n\n");
	printf("The velocity of particles are as given for every time step taken:\n");
	for(j=0;j<99;j++){
		for(i=1;i<=n;i++){
			printf("%f\t",v[j][i]);
		}
		printf("\n");
	}
	printf("\n\n");
	for(j=0;j<99;j++){
		mx[j]=0;
		mv[j]=0;
	}
	for(j=0;j<99;j++){
		for(i=1;i<=n;i++){
			mx[j]=mx[j]+x[j][i];
			mv[j]=mv[j]+v[j][i];
		}
	}
	printf("The displacement and velocity of the center of mass of every particle:\n");
	printf("COM Displacement:\n");
	for(j=0;j<99;j++){
		printf("%f\n",mx[j]);
	}
	printf("\nCOM Velocity:\n");
	for(j=0;j<99;j++){
		printf("%f\n",mv[j]);
	}
	return 0;
} 
