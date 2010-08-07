

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

double X[N][N+20];
double A[N][N+20];
double B[N][N+20];

void init_arrays()
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

double rtclock()
{
  struct timezone tzp;
  struct timeval tp;
  int stat;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main()
{
  init_arrays();

  double annot_t_start=0, annot_t_end=0, annot_t_total=0;
  int annot_i;

  for (annot_i=0; annot_i<REPS; annot_i++)
  {
    annot_t_start = rtclock();
    
    int i1,i2,t;
     
    for (t=0; t<=T-1; t++)
      {
	
	for (i1=0; i1<=N-1; i1++)
	  for (i2=1; i2<=N-1; i2++)
	    {
	      X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
	      B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
	    }

	for (i1=1; i1<=N-1; i1++)
	  for (i2=0; i2<=N-1; i2++)
	    {
	      X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
	      B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
	    }
      }

    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return ((int) B[0][0]); 

}
                                    
