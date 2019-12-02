#include <stdio.h>

void solve_tridiagonal(long double * restrict const x, const int X, const long double * restrict const a, const long double * restrict const b, long double * restrict const c) {
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
        const long double m = 1.0f / (b[ix] - a[ix] * c[ix - 1]);
        c[ix] = c[ix] * m;
        x[ix] = (x[ix] - a[ix] * x[ix - 1]) * m;
    }
    
    /* loop from X - 2 to 0 inclusive (safely testing loop condition for an unsigned integer), to perform the back substitution */
    for (int ix = X - 2; ix>=0 ; ix--)
        x[ix] -= c[ix] * x[ix + 1];
}

const long double A = 1;
const long double B = 3;
const long double Du = 1;
const long double Dv = 100;
                                                                                                                                                                                                                                                                                                                                                       
int main(){

long double a = 0;
long double b = 1;
long double T = 1;
long double t0 = 0;
long double h = 0.01;
long double tau = 0.01;
int Nx = (b - a) / h;
int Nt = (T - t0) / tau;

long double ru = Du*tau/2/h/h;
long double rv = Dv*tau/2/h/h;

long double au[Nx-2]; long double bu[Nx-2]; long double cu[Nx-2];
long double av[Nx-2]; long double bv[Nx-2]; long double cv[Nx-2];
au[0] = 0.; av[0] = 0.;
for(int i=1; i<Nx-2; i++){
	av[i] = -rv;
	au[i] = -ru;
}
for(int i=0; i<Nx-2; i++){
	bu[i] = 1+2*ru;
	bv[i] = 1+2*rv;
}
cu[Nx-2] = 0.;
cv[Nx-2] = 0.;
for(int i=0; i<Nx-3; i++){
	cv[i] = -rv;
	cu[i] = -ru;
}

long double fu[Nt][Nx]; 
long double fv[Nt][Nx];

for(int j = 0; j < Nt; j++){
	for(int i = 0; i < Nx; i++){
		fu[j][i] = 0.;
		fv[j][i] = 0.;
	}
}
                                                                                              
for(int i = 0; i < Nx; i++){
	fv[0][i] = B/A;
	fu[0][i] = A;
}

long double urightpart[Nx-2];
long double vrightpart[Nx-2];

for(int j = 1; j < Nt; j++){
	for(int i = 1; i < Nx-1; i++){
		urightpart[i] = (1-tau*(B+1)-2*ru)*fu[j-1][i] + ru*(fu[j-1][i-1] + fu[j-1][i+1]) + tau*fu[j-1][i]*fu[j-1][i]*fv[j-1][i] + tau*A; 
		vrightpart[i] = (1 -2*rv)*fv[j-1][i] + tau*B*fu[j-1][i] - tau*fu[j-1][i]*fu[j-1][i]*fv[j-1][i];
	}
//	if(j == 3) for(int i=1; i < Nx-1; i++){ printf("%f ",urightpart[i]); }
	solve_tridiagonal(urightpart, Nx-2, au, bu, cu);
	solve_tridiagonal(vrightpart, Nx-2, av, bv, cv);
//	if(j == 3) for(int i=1; i < Nx-1; i++){ printf("%f ",urightpart[i]); }
	for(int i = 1; i < Nx-1; i++){
		fu[j][i] = urightpart[i];
		fv[j][i] = vrightpart[i];
	}

	fu[j][0] = fu[j][1];
	fu[j][Nx-1] = fu[j][Nx-2];
	fv[j][0] = fv[j][1];
	fv[j][Nx-1] = fv[j][Nx-2];


/*
	fu[j][0] = fu[j-1][0] + tau*(A + fu[j-1][0]*fu[j-1][0]*fv[j-1][0]- (B+1)*fu[j-1][0] + 2*Du/h/h*(fu[j-1][1] - fu[j-1][0]));
	fu[j][Nx-1] = fu[j-1][Nx-1] + tau*(A + fu[j-1][Nx-1]*fu[j-1][Nx-1]*fv[j-1][Nx-1] - (B+1)*fu[j-1][Nx-1] + 2*Du/h/h*(fu[j-1][Nx-2] - fu[j-1][Nx-1]));
	fv[j][0] = fv[j-1][0] + tau*(-fu[j-1][0]*fu[j-1][0]*fv[j-1][0] + B*fu[j-1][0] + 2*Dv/h/h*(fv[j-1][1] - fv[j-1][0]));
	fv[j][Nx-1] = fv[j-1][Nx-1] + tau*(-fu[j-1][Nx-1]*fu[j-1][Nx-1]*fv[j-1][Nx-1] + B*fu[j-1][Nx-1] + 2*Dv/h/h*(fv[j-1][Nx-2] - fv[j-1][Nx-1]));*/
}	

FILE* fp;
fp = fopen("brus.dat", "w");
for(int j = 0; j < Nt; j++){
	for(int i = 0; i < Nx; i++){
		fprintf(fp, "%Lf %Lf %Lf %Lf \n", j*tau, i*h, fu[j][i], fv[j][i]);
	}
}
fclose(fp);
return 0;
}

























































