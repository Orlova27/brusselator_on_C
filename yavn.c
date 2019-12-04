#include <stdio.h>
#include <stdlib.h>

const float A = 1;
const float B = 1;
const float Du = 1;
const float Dv = 100;
                                                                                                                                                                                                                                                                                                                                                       
int main(){

float a = 0;
float b = 1;
float T = 1;
float t0 = 0;
float h = 0.01;
float tau = 0.01;
int Nx = (b - a) / h;
int Nt = (T - t0) / tau;

float ru = Du*tau/2/h/h;
float rv = Dv*tau/2/h/h;

float fu[Nt+1][Nx+1]; 
float fv[Nt+1][Nx+1];

for(int j = 0; j < Nt; j++){
	for(int i = 0; i < Nx; i++){
		fu[j][i] = 0.;
		fv[j][i] = 0.;
	}
}
                                                                                              
    double noiseu = 0.005*(rand()%100);    
    double noisev = 0.005*(rand()%100);                                                                             
for(int i = 0; i <= Nx; i++){
	fv[0][i] = B/A + noisev; 
	fu[0][i] = A + noiseu;
}


for(int j = 1; j <= Nt; j++){
for(int i = 1; i <Nx; i++){
	fu[j][i] = fu[j-1][i] + tau*(A+fu[j-1][i]*fu[j-1][i]*fv[j-1][i] - (B+1)*fu[j-1][i] + Du/h/h*(fu[j-1][i+1]+fu[j-1][i-1]-2*fu[j-1][i]));
	fv[j][i] = fv[j-1][i] + tau*(-fu[j-1][i]*fu[j-1][i]*fv[j-1][i] + B*fu[j-1][i] + Dv/h/h*(fv[j-1][i+1]+fv[j-1][i-1]-2*fv[j-1][i]));
}
	fu[j][0] = fu[j][1];
	fv[j][0] = fv[j][1];
	fv[j][Nx] = fv[j][Nx-1];
	fu[j][Nx] = fu[j][Nx-1];
}

FILE* fp;
fp = fopen("brus.dat", "w");
for(int j = 0; j <= Nt; j++){
	for(int i = 0; i <= Nx; i++){
		fprintf(fp, "%f %f %f %f \n", j*tau, i*h, fu[j][i], fv[j][i]);
	}
}
fclose(fp);
return 0;
}

























































