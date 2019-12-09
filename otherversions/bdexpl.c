#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const float A = 1;
const float B = 1.5;
const float Du = 1;
const float Dv = 100;
                                                                                                                                                                                                                                                                                                                                                       
int main(){

FILE* fp;
fp = fopen("explicit.dat", "w");

float a = 0;
float b = 182;
float T = 100;
float t0 = 0;
float h = 0.6;
float tau = 0.001;
int Nx = (b - a) / h;
int Nt = (T - t0) / tau;

float ru = Du*tau/2/h/h;
float rv = Dv*tau/2/h/h;

float fu[2][Nx+1]; 
float fv[2][Nx+1];
int stime;
long ltime;

ltime = time (NULL);
stime = (unsigned int) ltime/2;
srand(stime);
int count = 0;

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


while(count*tau<T){
for(int i = 1; i <Nx; i++){
	fu[1][i] = fu[0][i] + tau*(A+fu[0][i]*fu[0][i]*fv[0][i] - (B+1)*fu[0][i]  + Du/h/h*(fu[0][i+1]+fu[0][i-1]-2*fu[0][i]));
	fv[1][i] = fv[0][i] + tau*(-fu[0][i]*fu[0][i]*fv[0][i] + B*fu[0][i] + Dv/h/h*(fv[0][i+1]+fv[0][i-1]-2*fv[0][i]));
}
	fu[1][0] = fu[1][1];
	fv[1][0] = fv[1][1];
	fv[1][Nx] = fv[1][Nx-1];
	fu[1][Nx] = fu[1][Nx-1];
	count ++;

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

























































