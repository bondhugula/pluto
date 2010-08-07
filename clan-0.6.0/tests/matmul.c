/* matmul.c  128*128 matrix multiply */
#include <stdio.h>
#define N 128

int main()
{
  int   i,j,k;
  float a[N][N], b[N][N], c[N][N]; 

  /* We read matrix a */
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      scanf(" %f",&a[i][j]);

  /* We read matrix b */
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
      scanf(" %f",&b[i][j]);

  /* c = a * b */
#pragma scop
  for(i=0; i<N; i++)
    for(j=0; j<N; j++)
    {
      c[i][j] = 0.0;
      for(k=0; k<N; k++)
        c[i][j] = c[i][j] + a[i][k]*b[k][j];
    }
#pragma endscop

  /* We print matrix c */
  for(i=0; i<N; i++)
  {
    for(j=0; j<N; j++)
      printf("%6.2f ",c[i][j]);
    printf("\n");
  }
  
  return 0;
}
