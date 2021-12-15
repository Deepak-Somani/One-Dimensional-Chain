// Linear Harmonic motion
// By Deepak Somani
// Runge-Kutta 4th order

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<string.h>

// k=1,m=1;

#define n 32
// No. of particles

#define h 0.0005
// Time Step

#define r 1000
// Total no. of time steps

#define Tm 1.0
// Mean Temperature

#define Th 1.1
#define Tc 0.9
// Left and Right extreme constant Temperature

float rkmom(int i, float q1, float q0, float q2, float ns1, float ns2, float v1, float v2){
	if(i==1) h*(-2*q1 + q0 + q2 - ns1*v1);
	else if(i==n) h*(-2*q1 + q0 + q2 - ns2*v2);
	else return h*(-2*q1 + q0 + q2);
}
// Runge-Kutta time-step function for momentum

float rkpos(float p){
	return h*p;
}
// Runge-Kutta time-step function for position

float rknos(float p, int i){
	if(i==1) return h*(pow(p,2)/Th - 1);
	else return h*(pow(p,2)/Tc - 1);
}
// Runge-Kutta time-step function for noise terms

int main(){
	float q[n+1][r+1],p[n+1][r+1],P1[n+1],P2[n+1],P3[n+1],P4[n+1],Q1[n+1],Q2[n+1],Q3[n+1],Q4[n+1],T[n+1],zl[r+1],zr[r+1],Z1[5],Z2[5],j[n+1],J,kappa;	
	int i,x;
	
	for(x=0;x<r+1;x++){
		q[0][x]=-1;
		q[n+1][x]=n;
		p[0][x]=0;
		p[n+1][x]=0;
	}
	// Assigning constant position and 0 velocity to fixed particles at left and right end
	
	srand((unsigned int)time(NULL));
	
	for(i=1;i<=n;i++){
		q[i][0]=i-1;
		p[i][0]= -0.5 + ((float)rand()/(float)(RAND_MAX)) * 1.0;
		T[i] = 0.9 + sin(3.14*i/n)* 0.2;
	}
	// Assigning initial position and random initial velocity to each particle
	
	zl[0]=0;
	zr[0]=0;
	// Assigning initial conditions to noise terms
	
	for(i=1;i<=n;i++){
		j[i] = 0;
	}
	
	for(x=0;x<r;x++){
	
		for(i=1;i<=n;i++){
			P1[i] = rkmom(i,q[i][x],q[i-1][x],q[i+1][x],zl[x],zr[x],p[1][x],p[n][x]);
			Q1[i] = rkpos(p[i][x]);
			if(i==1) Z1[1] = rknos(p[i][x],i);
			if(i==n) Z2[1] = rknos(p[i][x],i);
		}
		
		for(i=1;i<=n;i++){
			P2[i] = rkmom(i,q[i][x] + Q1[i]/2,q[i-1][x] + Q1[i-1]/2,q[i+1][x] + Q1[i+1]/2,zl[x] + Z1[1]/2,zr[x] + Z2[1]/2,p[1][x] + P1[1]/2,p[n][x] + P1[n]/2);
			Q2[i] = rkpos(p[i][x] + P1[i]/2);
			if(i==1) Z1[2] = rknos(p[i][x] + P1[i]/2,i);
			if(i==n) Z2[2] = rknos(p[i][x] + P1[i]/2,i);
		}
		
		for(i=1;i<=n;i++){
			P3[i] = rkmom(i,q[i][x] + Q2[i]/2,q[i-1][x] + Q2[i-1]/2,q[i+1][x] + Q2[i+1]/2,zl[x] + Z1[2]/2,zr[x] + Z2[2]/2,p[1][x] + P2[1]/2,p[n][x] + P2[n]/2);
			Q3[i] = rkpos(p[i][x] + P2[i]/2);
			if(i==1) Z1[3] = rknos(p[i][x] + P2[i]/2,i);
			if(i==n) Z2[3] = rknos(p[i][x] + P2[i]/2,i);
		}
		
		for(i=1;i<=n;i++){
			P4[i] = rkmom(i,q[i][x] + Q3[i],q[i-1][x] + Q3[i-1],q[i+1][x] + Q3[i+1],zl[x] + Z1[3],zr[x] + Z2[3],p[1][x] + P3[1],p[n][x] + P3[n]);
			Q4[i] = rkpos(p[i][x] + P3[i]);
			if(i==1) Z1[4] = rknos(p[i][x] + P3[i],i);
			if(i==n) Z2[4] = rknos(p[i][x] + P3[i],i);
		}
		
		for(i=1;i<=n;i++){
			p[i][x+1] = p[i][x] + (P1[i]+(2*P2[i])+(2*P3[i])+P4[i])/6;
			q[i][x+1] = q[i][x] + (Q1[i]+(2*Q2[i])+(2*Q3[i])+Q4[i])/6;
			zl[x+1] = zl[x] + (Z1[1]+(2*Z1[2])+(2*Z1[3])+Z1[4])/6;
			zr[x+1] = zr[x] + (Z2[1]+(2*Z2[2])+(2*Z2[3])+Z2[4])/6;
		}
		
		while(x>=r/2){
			for(i=1;i<=n;i++){
				j[i] = j[i] + (-q[i][x+1] + q[i-1][x+1] + q[i+1][x+1])*p[i][x+1]*(2/r);
				T[i] = T[i] + pow(p[i][x+1],2)*(2/r);
			}
		}
	}
	
	J = 0;
	// Energy
	
	for(i=1;i<=n;i++){
		J = J + j[i]/n;
	}
	
	kappa = J*n/(Th-Tc);
	// Heat Coeffecient
	
	printf("%f",kappa);
	
	return 0;
}


