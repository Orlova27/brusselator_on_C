#include <stdio.h>

/*
constant parameters
*/
const float A = 1;
const float B = 3;
const float Du = 8;
const float Dv = 10;
                                                                                                                                                                                                                                                                                                                                                       
int main() {

FILE* fp;
fp = fopen("raspru", "w");

float a = 0; // space grid start
float b = 1; // space grid end
float T = 20; // time interval
float t0 = 0; // time start
float h = 0.01; // space step
float tau = 0.01; // time step
int Nx = (b - a) / h; // number of space steps
int Nt = (T - t0) / tau; // number of time steps

float ru = Du*tau/2/h/h;  
float rv = Dv*tau/2/h/h;

float u[Nx*Nt];
float v[Nx*Nt];
/*
write to file solutions in format 
t x u v \n
starting from initial conditions u = A, v = B/A
*/

float t = 0; // time

for(int x = 0; x < Nx; x++)
{
	fprintf(fp, "%f %f", t, x*h);
	u[x] = A;
	fprintf(fp, " %f", u[x]);
	v[x] = B/A;
	fprintf(fp, " %f", v[x]);	
	fprintf(fp, "\n");
}

for(int k=1; k<Nt-1; k++) //Nt-1
{
	t += tau;
	u[Nx*k] = u[Nx*(k-1)] + tau*(A + u[Nx*(k-1)]*u[Nx*(k-1)]*v[Nx*(k-1)] - (B+1)*u[Nx*(k-1)]+2*Du/h/h*(u[Nx*(k-1)+1]-u[Nx*(k-1)]));
	u[Nx*k-1+Nx] = u[Nx*(k-1)-1+Nx] + tau*(A + u[Nx*(k-1)-1+Nx]*u[Nx*(k-1)-1+Nx]*v[Nx*(k-1)-1+Nx] - (B+1)*u[Nx*(k-1)-1+Nx] + Du*2/h/h*(u[Nx*(k-1)-2+Nx]-u[Nx*(k-1)-1+Nx])); 
	v[Nx*k] = v[Nx*(k-1)] + tau*(-u[Nx*(k-1)]*u[Nx*(k-1)]*v[Nx*(k-1)]+B*u[Nx*(k-1)]+Dv*2/h/h*(v[Nx*(k-1)+1]-v[Nx*(k-1)]));
	v[Nx*k-1+Nx] = v[Nx*(k-1)-1+Nx] + tau*(-u[Nx*(k-1)-1+Nx]*u[Nx*(k-1)-1+Nx]*v[Nx*(k-1)-1+Nx]+B*u[Nx*(k-1)-1+Nx]+Dv*2/h/h*(v[Nx*(k-1)-2+Nx]-v[Nx*(k-1)-1+Nx]));
	float rightpartu[Nx-2]; // right part for crank-nicholson for u
	float rightpartv[Nx-2]; // right part for crank-nicholson for v
	for(int i = 0; i<Nx-2; i++) {
		rightpartu[i] = tau*A + tau*v[(k-1)*Nx+i+1]*u[(k-1)*Nx+i+1]*u[(k-1)*Nx+i+1]+ru*(u[(k-1)*Nx+i+1-1]+u[(k-1)*Nx+i+1+1])+u[(k-1)*Nx+i+1]*(1-tau*(B+1)-2*ru);
		rightpartv[i] = (1-2*rv)*v[(k-1)*Nx+i+1]+tau*B*u[(k-1)*Nx+i+1]-tau*u[(k-1)*Nx+i+1]*u[(k-1)*Nx+i+1]*v[(k-1)*Nx+i+1];
	}
	float alphau[Nx-2];
	float betau[Nx-2];
	alphau[0] = ru/(1+2*ru);
	betau[0] = rightpartu[0]/(1+2*ru);
	for(int n=1; n<=Nx-4; n++){
			alphau[n] = ru/(-ru*alphau[n-1]+1+2*ru);
			betau[n] = (rightpartu[n-1]+ru*betau[n-1])/(-ru*alphau[n-1]+1+2*ru);
	} // coefficients in thomas method
	u[Nx*k-2+Nx] = (rightpartu[Nx-3]+ru*betau[Nx-4])/(1+2*ru-ru*alphau[Nx-4]); // last non-border solution
	for(int m=Nx-4; m>=0; m--){
		u[Nx*k+1+m] = u[Nx*k+2+m]*alphau[m]+betau[m];
	} // evaluate one layer for u
	float alphav[Nx-2];
	float betav[Nx-2];
	alphav[0] = rv/(1+2*rv);
	betav[0] = rightpartv[0]/(1+2*rv);
	for(int n=1; n<=Nx-4; n++){
			alphav[n] = rv/(-rv*alphav[n-1]+1+2*rv);
			betav[n] = (rightpartv[n-1]+rv*betav[n-1])/(-rv*alphav[n-1]+1+2*rv);
	} // coefficients in thomas method
	v[Nx*k-2+Nx] = (rightpartv[Nx-3]+rv*betav[Nx-4])/(1+2*rv-rv*alphav[Nx-4]); // last non-border solution
	for(int m=Nx-4; m>=0; m--){
		v[Nx*k+1+m] = v[Nx*k+2+m]*alphav[m]+betav[m];
	} // evaluate one layer for v

	for(int j = 0; j<Nx; j++) {
	fprintf(fp, "%f %f %f %f \n", t, j*h, u[Nx*k+j], v[Nx*k+j]);
	}

}
fclose(fp);
return 0;
}
