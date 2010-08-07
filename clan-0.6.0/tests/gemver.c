#include <stdio.h>
#define n 100

int main()
{ int i=0, j=0, N=n, alpha, beta;
  float A[n][n], u1[n], u2[n], v1[n], v2[n], w[n], x[n], y[n], z[n];

  /* GEMVER kernel */
#pragma scop
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      A[i][j] = A[i][j] + u1[i]*v1[j] + u2[i]*v2[j];

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      x[i] = x[i] + beta* A[i][j]*y[j];


  for (i=0; i<N; i++)
    x[i] = x[i] + z[i];

  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      w[i] = w[i] + alpha* A[i][j]*x[j];
#pragma endscop
  
  return 0;
}
