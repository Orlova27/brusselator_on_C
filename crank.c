#include <stdio.h>

/*
constant parameters
*/
const float A = 1;
const float B = 3;
const float Du = 1;
const float Dv = 100;
                                                                                                                                                                                                                                                                                                                                                       
int main() {

FILE* fp;
/*
open file for writing solutions for u
*/
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
/*
float au[Nt][Nx-1];
for(int i = 0; i <= Nx-1; i++) { for(int j = 0; j< Nt; j++) {au[j][i] = -ru;}}
float bu[Nt][Nx];
for(int i = 0; i <= Nx; i++) { for(int j = 0; j < Nt; j++) {bu[j][i] = 1 + 2*ru;}} 
float cu[Nt][Nx-1];
for(int i = 0; i < Nx; i++){ for(int j = 0; j< Nt; j++) {cu[j][i] = -ru;}}
float du[Nt][Nx];
for(int j = 0; j < Nt; j++){for(int i = 0; i <= Nx; i++) {du[j][i] = 0;}}
float fu[Nt][Nx];
for(int j = 0; j < Nt; j++){for(int i = 0; i <= Nx; i++) { fu[j][i] = 0;}}
*/
/*
float av[Nt][Nx-1];
for(int i = 0; i <= Nx-1; i++) { for(int j = 0; j< Nt; j++) {av[j][i] = -rv;}}
float bv[Nt][Nx];
for(int i = 0; i <= Nx; i++) { for(int j = 0; j < Nt; j++) {bv[j][i] = 1 + 2*rv;}} 
float cv[Nt][Nx-1];
for(int i = 0; i <= Nx; i++){ for(int j = 0; j< Nt; j++) {cv[j][i] = -rv;}}*/
float dv[Nt][Nx];
for(int j = 0; j < Nt; j++){
for(int i = 0; i <= Nx; i++) {dv[j][i] = 0;}}

float fv[Nt][Nx];
for(int j = 0; j < Nt; j++){
for(int i = 0; i <= Nx; i++) { fv[j][i] = 0;}}


for(int i = 0; i <= Nt; i++) {fv[i][0] = A/B;}
for(int i = 0; i <= Nt; i++) {fu[i][0] = A;}                                             
                                                                                                                               

for(int j = 1; j < Nt; j++){                                                                                                                          
for (int i=1; i<Nx-1; i++) 
{
    du[j][i] = (1 - tau*(B+1) -2*ru)*fu[j][i] + ru*(fu[j][i-1] + fu[j][i+1]) + tau*fu[j][i]*fu[j][i]*fv[j][i] + tau*A;

    /*dv[j][i] = (1 -2*rv)*fv[j][i] + tau*B*fu[j][i] + tau*fu[j][i]*fu[j][i]*fv[j][i];*/

}}

for(int j = 0; j < Nt; j++){fu[j][0] = 1;}

for(int i = 1; i < Nx-1; i++){
for(int j = 0; j < Nt; j++){                                                                                                                                                                           
  float c_star[Nt][Nx];
  float d_star[Nt][Nx];
  /*for(int k = 0; k <= Nx -1; k++) {for(int l = 0; l < Nt; l++){ c_star[l][k] = 0;}}
  for(int k = 0; k <= Nx -1; k++) {for(int l = 0; l < Nt; l++){ d_star[l][k] = 0;}}*/
                                                                                                                                                      
  c_star[j][0] = cu[j][0] / bu[j][0];
  d_star[j][0] = du[j][0] / bu[j][0];
                                                                                                                                             
  for (int p=1; p<Nx; p++) 
	{
	    float m = 1.0 / (bu[j][p] - au[j][p] * c_star[j][p-1]);
	    c_star[j][p] = cu[j][p] * m;
	    d_star[j][p] = (du[j][p] - au[j][p] * d_star[j][p-1]) * m;
	}
                                                                                                                                               
  for (int p=Nx; p> 0; p--) 
	{
	if(p==1 ||p == Nx){fu[j+1][p] = fu[j+1][p-1];}	    
	fu[j+1][p] = d_star[j][p] - c_star[j][p] * du[j][p+1];
	}
	
}
}

for(int j = 0; j< Nt; j++){
for (int i=0; i<Nx; i++) 
{
	 fprintf(fp, "%f %f %f \n", j*0.01, i*0.01, fu[j][i]);
	
}
fprintf(fp, "\n");}
fclose(fp);
/* Now for v */
FILE* fpv;
fpv = fopen("rasprv", "w");
for(int i = 0; i <= Nx-1; i++) { for(int j = 0; j< Nt; j++) {au[j][i] = -rv;}}

for(int i = 0; i <= Nx; i++) { for(int j = 0; j < Nt; j++) {bu[j][i] = 1 + 2*rv;}} 

for(int i = 0; i <= Nx; i++){ for(int j = 0; j< Nt; j++) {cu[j][i] = -rv;}}

for(int j = 0; j < Nt; j++){for(int i = 0; i <= Nx; i++) {dv[j][i] = 0;}}

for(int j = 1; j < Nt; j++){for (int i=1; i<Nx-1; i++) 
{
    dv[j][i] = (1 -2*rv)*fv[j][i] + tau*B*fu[j][i] + tau*fu[j][i]*fu[j][i]*fv[j][i];

}}

for(int i = 1; i < Nx-1; i++){
for(int j = 0; j < Nt; j++){                                                                                                                                                                           
  float c_star[Nt][Nx];
  float d_star[Nt][Nx];
  /*for(int k = 0; k <= Nx -1; k++) {for(int l = 0; l < Nt; l++){ c_star[l][k] = 0;}}
  for(int k = 0; k <= Nx -1; k++) {for(int l = 0; l < Nt; l++){ d_star[l][k] = 0;}}*/
                                                                                                                                                      
  c_star[j][0] = cu[j][0] / bu[j][0];
  d_star[j][0] = dv[j][0] / bu[j][0];
                                                                                                                                             
  for (int p=1; p<Nx; p++) 
	{
	    float m = 1.0 / (bu[j][p] - au[j][p] * c_star[j][p-1]);
	    c_star[j][p] = cu[j][p] * m;
	    d_star[j][p] = (dv[j][p] - au[j][p] * d_star[j][p-1]) * m;
	}
            fu[j+1][1] = fu[j+1][0];                                                                                                                                   
  for (int p=Nx-1; p> 1; p--) 
	{
	    fu[j+1][p] = d_star[j][p] - c_star[j][p] * dv[j][p+1];
	}
	
	fu[j+1][Nx] = fu[j+1][Nx-1];  
}
}

for(int j = 0; j< Nt; j++){
for (int i=0; i<Nx; i++) 
{
	 fprintf(fpv, "%f \n", fu[j][i]);
	
}
fprintf(fpv, "\n");}
/*fprintf(fp, "\n");

for( int j = 0; j<Nt; j++){
for (int i=0; i<Nx; i++) 
{
	 fprintf(fp, "%f \n", fv[j][i]);
}
fprintf(fp, "\n");}
*/
fclose(fpv);
return 0;
}
