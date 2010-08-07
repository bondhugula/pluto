#include <stdio.h> 
#include <stdlib.h> 

#define N 5000
#define alpha 1
#define beta 1

double A[N][N];
double B[N][N];
double x[N];
double u1[N];
double u2[N];
double v2[N];
double v1[N];
double w[N];
double y[N];
double z[N];

void init_array()  
{  
  int i, j;
  for (i=0; i<N; i++) {
    u1[i] = i;
    u2[i] = (i+1)/N/2.0;
    v1[i] = (i+1)/N/4.0;
    v2[i] = (i+1)/N/6.0;
    y[i] = (i+1)/N/8.0;
    z[i] = (i+1)/N/9.0;
    x[i] = 0.0;
    w[i] = 0.0;
    for (j=0; j<N; j++) {
      A[i][j] = ((double) i*j)/N;
    }
  }
}  

int main() 
{ 
  init_array(); 
 
  /*@ profiled code @*/ 
 
  return w[0];   //needed to avoid dead code elimination 
} 
