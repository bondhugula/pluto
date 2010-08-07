#include <stdio.h>
#define N 100

int main()
{ int i=0, j=0, n=100 ;
  float a[N+1][N+1], b[N+1], c[N+1], result ;

  /* ax-do kernel */
#pragma scop
  for (i=1;i<=n;i++)                 
    c[i] = 0 ;                       
  for (i=1;i<=n;i++)                
    for (j=1;j<=n;j++)              
      c[i] = c[i] + a[i][j] * b[j] ; 
#pragma endscop
  
  result = c[N-1];
  printf("fib[%d] = %d\n", N-1, result);

  return 0;
}
