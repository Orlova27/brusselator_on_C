#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

#define PREF "brus2da45b75"
#define SUFF ".dat"

const float A = 4.5;
const float B = 7.5;
const float Dux = 2;
const float Dvx = 16;
const float Duy = 2;
const float Dvy = 16;
                       
float F1(float u, float v)
{
	return (A + u*u*v-(B+1)*u);
}

float F2(float u, float v)
{
	return (-u*u*v + B*u);

}

float runge_kutt(float h, float u, float v, int f)
{	
	float k11 = 0, k12 = 0, k13 = 0, k14 = 0;
	float k21 = 0, k22 = 0, k23 = 0, k24 = 0;
	k11 = h*F1(u, v);
	k21 = h*F2(u, v);
	k12 = h*F1(u + k11/2, v + k21/2);
	k22 = h*F2(u + k11/2, v + k21/2);
	k13 = h*F1(u + k12/2, v + k22/2);
	k23 = h*F2(u + k12/2, v + k22/2);
	k14 = h*F1(u + k13/2, v + k23/2);
	k24 = h*F2(u + k13/2, v + k23/2);
	u = u + (k11 + 2*k12 + 2*k13 + k14)/6;
	v = v + (k21 + 2*k22 + 2*k23 + k24)/6;

	if(f==1) return u;
	if(f==2) return v;
}


void solve_tridiagonal(float * x, int X, float * a, float * b, float * c, float* d) {
    /*
     solves Ax = d where A is a tridiagonal matrix consisting of vectors a, b, c
     x - contains the input vector x indexed from 0 to X - 1 inclusive
     X - number of equations (length of vector x)
     a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive
     b - the main diagonal, indexed from 0 to X - 1 inclusive
     c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive
     d - right part of equation, indexed from 0 to X-1
     */
     
    float q[X];
    float p[X];

    p[0] = -c[0] / b[0];
    q[0] = -d[0] / b[0];
    
    for (int ix = 1; ix <= X; ix++) {
        float m = 1.0 / (-b[ix-1] - a[ix-1] * p[ix - 1]);
        p[ix] = c[ix-1] * m;
        q[ix] = (-d[ix-1] + a[ix-1] * q[ix - 1]) * m;
    }
    x[X-1] = q[X];
    
    for (int ix = X - 2; ix>=0 ; ix--)
        x[ix] = q[ix+1]+ p[ix+1] * x[ix + 1];
}
                                                                                                                                                                                                                                                                                           
int main(){

float x0 = 0;
float x1 = 30;
float y0 = 0;
float y1 = 30;
float T = 10;
float t0 = 0;
float h = 0.1;
float tau = 0.001;
int Ny = (y1 - y0)/h;
int Nx = (x1 - x0) / h;
int Nt = (T - t0) / tau;

int count = 0;

float rux = Dux*tau/2/h/h;
float rvx = Dvx*tau/2/h/h;
float ruy = Duy*tau/2/h/h;
float rvy = Dvy*tau/2/h/h;

float aux[Nx+1]; float bux[Nx+1]; float cux[Nx+1];
float avx[Nx+1]; float bvx[Nx+1]; float cvx[Nx+1];
float auy[Ny+1]; float buy[Ny+1]; float cuy[Ny+1];
float avy[Ny+1]; float bvy[Ny+1]; float cvy[Ny+1];

aux[0] = 0.; avx[0] = 0.;
aux[Nx] = -2*rux;
avx[Nx] = -2*rvx;
auy[0] = 0.; avy[0] = 0.;
auy[Ny] = -2*ruy;
avy[Ny] = -2*rvy;

for(int i=1; i<Nx; i++){
	avx[i] = -rvx;
	aux[i] = -rux;
}
for(int i=1; i<Ny; i++){
	avy[i] = -rvy;
	auy[i] = -ruy;
}
for(int i=0; i<=Nx; i++){
	bux[i] = 1+2*rux;
	bvx[i] = 1+2*rvx;
}
for(int i=0; i<=Ny; i++){
	buy[i] = 1+2*ruy;
	bvy[i] = 1+2*rvy;
}
cux[Nx] = 0.; cux[0] = -2*rux;
cvx[Nx] = 0.; cvx[0] = -2*rvx;
cuy[Ny] = 0.; cuy[0] = -2*ruy;
cvy[Ny] = 0.; cvy[0] = -2*rvy;
for(int i=1; i<Nx; i++){
	cvx[i] = -rvx;
	cux[i] = -rux;
}
for(int i=1; i<Ny; i++){
	cvy[i] = -rvy;
	cuy[i] = -ruy;
}

float fu[4][Nx+1][Ny+1]; 
float fv[4][Nx+1][Ny+1];

for(int j = 0; j < 4; j++){
	for(int k = 0; k <= Ny; k++){
		for(int i = 0; i <= Nx; i++){
			fu[j][i][k] = 0.;
			fv[j][i][k] = 0.;
		}	
	}
}
 /* float noiseu = 0.005*(rand()%100);    
  float noisev = 0.005*(rand()%100);        */                                                                     
for(int i = 0; i <= Nx; i++){
for(int k = 0; k <= Ny; k++){
	fv[0][i][k] = B/A + 0.0001*(rand()%100); //+ noiseu; 
	fu[0][i][k] = A + 0.0001*(rand()%100); //+ noisev;
}}

// НАДО СОЗДАВАТЬ МАССИВЫ РАЗМЕРА MAX(NX,NY)+1
float urightpart[Nx+1];
float vrightpart[Nx+1];
float u[Nx+1];
float v[Nx+1];

for(int i = 0; i <= Nx; i++){
	v[i] = 0;
	u[i] = 0;
}

int number_of_files = Nt+1;

char buf[BUFSIZ];
FILE* fp; 

while(count*tau<T){
	
	for(int i = 0; i<=Nx; i++){
		for(int k = 0; k<=Ny; k++){
		fu[1][i][k] = runge_kutt(tau, fu[0][i][k], fv[0][i][k], 1);
		fv[1][i][k] = runge_kutt(tau, fu[0][i][k], fv[0][i][k], 2);
		}
	}
	for(int k = 0; k<=Ny; k++){
		urightpart[0] = fu[1][0][k] + 2*rux*(fu[0][1][k] - fu[0][0][k]);
		vrightpart[0] = fv[1][0][k] + 2*rvx*(fv[0][1][k] - fv[0][0][k]);
 		urightpart[Nx] = fu[1][Nx][k] + 2*rux*(fu[0][Nx-1][k] - fu[0][Nx][k]);
		vrightpart[Nx] = fv[1][Nx][k] + 2*rvx*(fv[0][Nx-1][k] - fv[0][Nx][k]);
		for(int i = 1; i < Nx; i++){
			urightpart[i] = fu[1][i][k]+rux*(fu[0][i+1][k] - 2*fu[0][i][k] + fu[0][i-1][k]); 
			vrightpart[i] = fv[1][i][k]+rvx*(fv[0][i+1][k] - 2*fv[0][i][k] + fv[0][i-1][k]);
		}
		solve_tridiagonal(u, Nx+1, aux, bux, cux, urightpart);
		solve_tridiagonal(v, Nx+1, avx, bvx, cvx, vrightpart);
		for(int i = 0; i <= Nx; i++){
			fu[2][i][k] = u[i];
			fv[2][i][k] = v[i];
		}
	}
	
	for(int i=0; i<=Nx; i++){
		urightpart[0] = fu[2][i][0] + 2*ruy*(fu[1][i][1] - fu[1][i][0]);
		vrightpart[0] = fv[2][i][0] + 2*rvy*(fv[1][i][1] - fv[1][i][0]);
 		urightpart[Ny] = fu[2][i][Ny] + 2*ruy*(fu[1][i][Ny-1] - fu[1][i][Ny]);
		vrightpart[Ny] = fv[2][i][Ny] + 2*rvy*(fv[1][i][Ny-1] - fv[1][i][Ny]);
		for(int k = 1; k < Ny; k++){
			urightpart[k] = fu[2][i][k]+ruy*(fu[1][i][k+1] - 2*fu[1][i][k] + fu[1][i][k-1]); 
			vrightpart[k] = fv[2][i][k]+rvy*(fv[1][i][k+1] - 2*fv[1][i][k] + fv[1][i][k-1]);
		}
		solve_tridiagonal(u, Ny+1, auy, buy, cuy, urightpart);
		solve_tridiagonal(v, Ny+1, avy, bvy, cvy, vrightpart);
		for(int k = 0; k<=Ny; k++){
		fu[3][i][k] = u[k];
		fv[3][i][k] = v[k];
		}
	}
		
		if(count == 9000){       
		sprintf(buf, "%s%d%s", PREF, count, SUFF); 
		fp = fopen(buf, "w");
			for(int i = 0; i <= Nx; i++){
			for(int k = 0; k <= Ny; k++){
				fprintf(fp, "%f %f %f %f %f \n", count*tau, i*h, k*h, fu[3][i][k], fv[3][i][k]);
			}
			fprintf(fp, "\n");	
			}
		fclose(fp);
		}			

		
	for(int i = 0; i <= Nx; i++){
	for(int k = 0; k<=Ny; k++){
		fu[0][i][k] = fu[3][i][k]; 
		fv[0][i][k] = fv[3][i][k];
	}}

	for(int i= 0; i <BUFSIZ; i++){
		buf[i] = 0;
	}
	count++;
}

return 0;
}








































