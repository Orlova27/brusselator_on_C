#include <stdio.h>
void solve_tridiagonal_in_place_destructive(float * restrict const x, const int X, const float * restrict const a, const float * restrict const b, float * restrict const c) {
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

int main(){
//float x[5][5] = {{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1}};
//float a[5][5] = {{0,-1,-1,-1,-1},{0,-1,-1,-1,-1},{0,-1,-1,-1,-1},{0,-1,-1,-1,-1},{0,-1,-1,-1,-1}};
//float b[5][5] = {{3,3,3,3,3},{3,3,3,3,3},{3,3,3,3,3},{3,3,3,3,3},{3,3,3,3,3}};
float x[5] = {1,1,1,1,1};
float a[5] = {0,-1,-1,-1,-1};
float b[5] = {3,3,3,3,3};
float c[5][5] = {{-1,-1,-1,-1,0},{2,-1,-1,-1,0},{3,-1,-1,-1,0},{4,-1,-1,-1,0},{5,-1,-1,-1,0}};

float (*ptr)[5] = c;
//printf("%f ", *(*(ptr)+10));
int k=1;

/*solve_tridiagonal_in_place_destructive(x, 5, a, b, *(ptr)+k*5);
for(int i=0; i<=4; i++){
	printf("%f ", x[i]);
};*/
for(int k=0; k<=4; k++){
	for(int i=0; i<=4; i++){
		c[k][i] = b[i];
	};
}

	for(int i=0; i<=4; i++){ 
		for(int j=0; j<=4; j++){
			printf("%f ", c[i][j]); 
		};
		printf("\n");
	};

//printf("%f ", **(c+3*sizeof(float)));
/*
for(int k=0; k<=20;k+=5){
	printf("%d \n", k);
	solve_tridiagonal_in_place_destructive(*(x+sizeof(float)*k), 5, *(a+sizeof(float)*k), *(b+sizeof(float)*k), *(c+sizeof(float)*k));
	for(int i=0; i<=4; i++){ 
		for(int j=0; j<=4; j++){
			printf("%f ", x[i][j]); 
		};
		printf("\n");
	};
}
*/
return 0;
}
