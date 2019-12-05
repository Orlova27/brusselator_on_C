#include <stdio.h>
#include <stdlib.h>

#define mpi 3.1415926

void solve_tridiagonal(float * restrict const x, const int X, const float * restrict const a, const float * restrict const b, float * restrict const c) {
    /*
     solves Ax = v where A is a tridiagonal matrix consisting of vectors a, b, c
     x - initially contains the input vector v, and returns the solution x. indexed from 0 to X - 1 inclusive
     X - number of equations (length of vector x)
     a - subdiagonal (means it is the diagonal below the main diagonal), indexed from 1 to X - 1 inclusive
     b - the main diagonal, indexed from 0 to X - 1 inclusive
     c - superdiagonal (means it is the diagonal above the main diagonal), indexed from 0 to X - 2 inclusive
     
     Note: contents of input vector c will be modified, making this a one-time-use function (scratch space can be allocated instead for this purpose to make it reusable)
     Note 2: We don't check for diagonal dominance, etc.; this is not guaranteed stable
     */

    

    c[0] = c[0] / b[0];
    x[0] = x[0] / b[0];
    
   

    /* loop from 1 to X - 1 inclusive, performing the forward sweep */
    for (int ix = 1; ix < X; ix++) {
        const float m = 1.0f / (b[ix] - a[ix] * c[ix - 1]);
        c[ix] = c[ix] * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }
    
    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (int ix = X - 2; ix>=0 ; ix--)
        x[ix] -= c[ix] * x[ix + 1];
}

const float A = 1;
const float B = 1.5;
const float Du = 1;
const float Dv = 100;
                                                                                                                                                                                                                                                                                                                                                       
int main(){

float a = 0;
float b = 10;
float T = 1;
float t0 = 0;
float h = 1.0;
float tau = 0.004;
int Nx = (b - a) / h;
int Nt = (T - t0) / tau;

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

float fu[Nt+1][Nx+1]; 
float fv[Nt+1][Nx+1];

for(int j = 0; j <= Nt; j++){
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
	
	

for(int j = 1; j <= Nt; j++){

	urightpart[0] = fu[j-1][0] + 2*ru*(fu[j-1][1] - fu[j-1][0]) + tau*(A + fu[j-1][0]*fu[j-1][0]*fv[j-1][0]  - (B+1)*fu[j-1][0]);
	vrightpart[0] = fv[j-1][0] + 2*rv*(fv[j-1][1] - fv[j-1][0]) + tau*(-fu[j-1][0]*fu[j-1][0]*fv[j-1][0] + B*fu[j-1][0]);
 	urightpart[Nx] = fu[j-1][Nx] + 2*ru*(fu[j-1][Nx-1] - fu[j-1][Nx]) +  tau*(A + fu[j-1][Nx]*fu[j-1][Nx]*fv[j-1][Nx] - (B+1)*fu[j-1][Nx]);
	vrightpart[Nx] = fv[j-1][Nx] + 2*rv*(fv[j-1][Nx-1] - fv[j-1][Nx]) + tau*(-fu[j-1][Nx]*fu[j-1][Nx]*fv[j-1][Nx] + B*fu[j-1][Nx]);

	for(int i = 1; i < Nx; i++){
		urightpart[i] = fu[j-1][i]+ru*(fu[j-1][i+1] - 2*fu[j-1][i] + fu[j-1][i-1]) + tau*(A + fu[j-1][i]*fu[j-1][i]*fv[j-1][i] - (B+1)*fu[j-1][i]); 
		vrightpart[i] = fv[j-1][i]+rv*(fv[j-1][i+1] - 2*fv[j-1][i] + fv[j-1][i-1]) + tau*(-fu[j-1][i]*fu[j-1][i]*fv[j-1][i] + B*fu[j-1][i]);
	}

	solve_tridiagonal(urightpart, Nx+1, au, bu, cu);
	solve_tridiagonal(vrightpart, Nx+1, av, bv, cv);
	for(int i = 0; i < Nx+1; i++){
		fu[j][i] = urightpart[i];
		fv[j][i] = vrightpart[i];
	}

/*
	fu[j][0] = fu[j-1][0];// + tau*(A + fu[j][0]*fu[j][0]*fv[j][0]- (B+1)*fu[j][0] + 2*Du/h/h*(fu[j][1] - fu[j][0]));
	fu[j][Nx-1] = fu[j-1][Nx-1];// + tau*(A + fu[j][Nx-1]*fu[j][Nx-1]*fv[j][Nx-1] - (B+1)*fu[j][Nx-1] + 2*Du/h/h*(fu[j][Nx-2] - fu[j][Nx-1]));
	fv[j][0] = fv[j-1][0];// + tau*(-fu[j][0]*fu[j][0]*fv[j][0] + B*fu[j][0] + 2*Dv/h/h*(fv[j][1] - fv[j][0]));
	fv[j][Nx-1] = fv[j-1][Nx-1];// + tau*(-fu[j][Nx-1]*fu[j][Nx-1]*fv[j][Nx-1] + B*fu[j][Nx-1] + 2*Dv/h/h*(fv[j][Nx-2] - fv[j][Nx-1]));*/
}	

FILE* fp;
fp = fopen("brus15.dat", "w");
for(int j = 0; j <= Nt; j++){
	for(int i = 0; i <= Nx; i++){
		fprintf(fp, "%f %f %f %f \n", j*tau, i*h, fu[j][i], fv[j][i]);
	}
	fprintf(fp, "\n");
}
fclose(fp);
return 0;
}

























































