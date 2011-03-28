#include "util.h"

#define length 10

double t_start, t_end;

//dynamic programming
int main (int argc, char* argv[])
{
  int i,j,k;
  int sum_c[length][length][length];
  int c[length][length];
  int W[length][length]; //input
  int out[1]; //output
  int iter;
#define NITER 10000

	IF_TIME(t_start = rtclock());
  for (iter=0; iter<NITER; iter++) {

  /* pluto start (length) */
#pragma scop
  for( i=0; i<=length-1; i++)
    for( j=0; j<=length-1; j++)
      c[i][j] = 0;

  for( i=0; i<=length-2; i++)
  {
    for( j=i+1; j<=length-1; j++)
    {
      sum_c[i][j][i] = 0;
      for( k=i+1; k<=j-1; k++)
      {
        sum_c[i][j][k] = sum_c[i][j][k-1] + c[i][k] + c[k][j];
      }
      c[i][j] = sum_c[i][j][j-1] + W[i][j];
    }
  }
#pragma endscop

  out[0] += c[0][length-1];
  }
    IF_TIME(t_end = rtclock());
    IF_TIME(printf("%0.6lfs\n", t_end - t_start));

}
