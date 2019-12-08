#include <stdio.h>
#include <stdlib.h>

const float A = 1;
const float B = 1.97;
const float Du = 1;
const float Dv = 7;
                       
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

FILE* fp;
fp = fopen("newsolver.dat", "w");

float a = 0;
float b = 182;
float T = 20;
float t0 = 0;
float h = 0.6;
float tau = 0.001;
int Nx = (b - a) / h;
int Nt = (T - t0) / tau;

int count = 0;

float ru = Du*tau/2/h/h;
float rv = Dv*tau/2/h/h;

float au[Nx+1]; float bu[Nx+1]; float cu[Nx+1];
float av[Nx+1]; float bv[Nx+1]; float cv[Nx+1];
au[0] = 0.; av[0] = 0.;
au[Nx] = -2*ru;
av[Nx] = -2*rv;

for(int i=1; i<Nx; i++){
	av[i] = -rv;
	au[i] = -ru;
}
for(int i=0; i<=Nx; i++){
	bu[i] = 1+2*ru;
	bv[i] = 1+2*rv;
}

cu[Nx] = 0.;
cv[Nx] = 0.;
cu[0] = -2*ru;
cv[0] = -2*rv;
for(int i=1; i<Nx; i++){
	cv[i] = -rv;
	cu[i] = -ru;
}

float fu[2][Nx+1]; 
float fv[2][Nx+1];

for(int j = 0; j < 2; j++){
	for(int i = 0; i <= Nx; i++){
		fu[j][i] = 0.;
		fv[j][i] = 0.;
	}
}
    float noiseu = 0.005*(rand()%100);    
    float noisev = 0.005*(rand()%100);                                                                             
for(int i = 0; i <= Nx; i++){
	fv[0][i] = B/A + noisev; 
	fu[0][i] = A + noiseu;
}

float urightpart[Nx+1];
float vrightpart[Nx+1];
float u[Nx+1];
float v[Nx+1];

for(int i = 0; i <= Nx; i++){
	v[i] = 0;
	u[i] = 0;
}


while(count*tau<T){
	
	for(int i = 0; i<=Nx; i++){
	fu[1][i] = runge_kutt(tau, fu[0][i], fv[0][i], 1);
	fv[1][i] = runge_kutt(tau, fu[0][i], fv[0][i], 2);
	}

	urightpart[0] = fu[1][0] + 2*ru*(fu[0][1] - fu[0][0]);
	vrightpart[0] = fv[1][0] + 2*rv*(fv[0][1] - fv[0][0]);
 	urightpart[Nx] = fu[1][Nx] + 2*ru*(fu[0][Nx-1] - fu[0][Nx]);
	vrightpart[Nx] = fv[1][Nx] + 2*rv*(fv[0][Nx-1] - fv[0][Nx]);

	for(int i = 1; i < Nx; i++){
		urightpart[i] = fu[1][i]+ru*(fu[0][i+1] - 2*fu[0][i] + fu[0][i-1]); 
		vrightpart[i] = fv[1][i]+rv*(fv[0][i+1] - 2*fv[0][i] + fv[0][i-1]);
	}

	solve_tridiagonal(u, Nx+1, au, bu, cu, urightpart);
	solve_tridiagonal(v, Nx+1, av, bv, cv, vrightpart);
	for(int i = 0; i < Nx+1; i++){
		fu[1][i] = u[i];
		fv[1][i] = v[i];
	}
	count++;
	for(int i = 0; i <= Nx; i++){
		fprintf(fp, "%f %f %f %f \n", (count-1)*tau, i*h, fu[1][i], fv[1][i]);
	}
	fprintf(fp, "\n");

	for(int i = 0; i <= Nx; i++){
		fu[0][i] = fu[1][i]; 
		fv[0][i] = fv[1][i];
	}
}	

fclose(fp);
return 0;
}

























































