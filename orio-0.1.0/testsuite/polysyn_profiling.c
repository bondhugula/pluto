#include <stdio.h>
#include <stdlib.h>

#define N 500
double A[N][N];
double L[N][N]; 
double U[N][N]; 

void init_array() 
{ 
  int i, j, k; 
  for (i=0; i<N; i++)
    for (j=0; j<N; j++) { 
      L[i][j] = 0.0; 
      U[i][j] = 0.0; 
    } 
  for (i=0; i<N; i++)
    for (j=0; j<=i; j++) { 
      L[i][j] = i+j+1; 
      U[j][i] = i+j+1; 
    } 
  for (i=0; i<N; i++) 
    for (j=0; j<N; j++) 
      for (k=0; k<N; k++) 
	A[i][j] += L[i][k]*U[k][j];
} 

int main()
{
  init_array();

  /*@ profiled code @*/

  return A[0][0];   //needed to avoid dead code elimination
}
