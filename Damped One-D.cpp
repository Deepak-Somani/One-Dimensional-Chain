// Damped One-D chain
// Deepak Somani
// RK 4th order

#include<stdio.h>
#define h 0.01
// h is considered to be the step difference of time.
// The equation comes out as: dv[n]/dt = -(k/m)*(2*x[n] - x[n+1] + x[n-1]). - 1
// k is spring constant and m is mass of particle which is taken identical for all.
// and dx/dt = v. - 2
// (k/m) is taken as 10.

float velo(float x0, float x1, float x2, float v0, float v1, float v2){
	return -h*(10*(2*x1-x0-x2) + 1*(2*v1-v0-v2));
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
			k[0][i]=velo(x[j][i-1],x[j][i],x[j][i+1],v[j][i-1],v[j][i],v[j][i+1]);
		}
	
		for(i=1;i<=n;i++){
			l[1][i]=dist(v[j][i]+k[0][i]/2);
		}
	
		for(i=1;i<=n;i++){
			k[1][i]=velo(x[j][i-1]+(l[0][i-1]/2), x[j][i]+(l[0][i]/2), x[j][i+1]+(l[0][i+1]/2),v[j][i-1]+(k[0][i-1]/2), v[j][i]+(k[0][i]/2), v[j][i+1]+(k[0][i+1]/2) );
		}
	
		for(i=1;i<=n;i++){
			l[2][i]=dist(v[j][i]+k[1][i]/2);
		}
	
		for(i=1;i<=n;i++){
			k[2][i]=velo(x[j][i-1]+(l[1][i-1]/2), x[j][i]+(l[1][i]/2), x[j][i+1]+(l[1][i+1]/2),v[j][i-1]+(k[1][i-1]/2), v[j][i]+(k[1][i]/2), v[j][i+1]+(k[1][i+1]/2) );
		}
	
		for(i=1;i<=n;i++){
			l[3][i]=dist(v[j][i]+k[2][i]);
		}
	
		for(i=1;i<=n;i++){
			k[3][i]=velo(x[j][i-1]+(l[2][i-1]), x[j][i]+(l[2][i]), x[j][i+1]+(l[2][i+1]),v[j][i-1]+(k[2][i-1]), v[j][i]+(k[2][i]), v[j][i+1]+(k[2][i+1]) );
		}
		
		for(i=1;i<=n;i++){
			x[j+1][i]=x[j][i]+((l[0][i]+2*l[1][i]+2*l[2][i]+l[3][i])/6);
			v[j+1][i]=v[j][i]+((k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])/6);
			//at t= t+h every time the displacement and velocity of each particle gets updated.
		}
	}
	
	
	printf("\n");
	for(j=0;j<99;j++){
		for(i=1;i<=n;i++){
			printf("%f\t",x[j][i]);
		}
		printf("\n");
	}
	printf("\n\n");
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
	for(j=0;j<99;j++){
		printf("%f\n",mx[j]/n);
	}
	printf("\n\n");
	for(j=0;j<99;j++){
		printf("%f\n",mv[j]/n);
	}
	return 0;
} 
