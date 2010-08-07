

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define alpha 1
#define beta 1
double A[N][N +20];
double B[N][N +20];
double x[N];
double u1[N];
double u2[N];
double v1[N];
double v2[N];
double w[N];
double y[N];
double z[N];

void init_arrays()
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

    int i,j;
     
for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    B[i][j] = A[i][j] + u1[i]*v1[j] + u2[i]*v2[j];

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    x[i] = x[i] + beta* B[j][i]*y[j];


for (i=0; i<=N-1; i++)
  x[i] = x[i] + z[i];

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    w[i] = w[i] + alpha* B[i][j]*x[j];

    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);

  return ((int) w[0]); 

}
                                    

