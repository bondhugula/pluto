#include <stdio.h> 
#include <stdlib.h> 

#define T 400
#define N 400
double X[N][N+20];
double A[N][N+20];
double B[N][N+20];

void init_array()  
{  
  int i, j;
  for (i=0; i<N; i++) 
    for (j=0; j<N; j++) 
      {
	A[i][j] = (1+(i*j)%1024)/2.0;
	B[i][j] = (1+(i*j)%1024)/3.0;
	X[i][j] = (1+(i*j)%1024)/3.0;
      }
}  
 
int main() 
{ 
  init_array(); 
 
  /*@ profiled code @*/ 
 
  return B[0][0];   //needed to avoid dead code elimination 
} 
