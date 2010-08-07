

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

double a[T][N];

void init_input_vars()
{
  int i, j;
  for (i=0; i<T; i++) 
    for (j=0; j<N; j++) 
      a[i][j] = i+((double)j)/N;
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
  init_input_vars();


  double orio_t_start=0, orio_t_end=0, orio_t_total=0;
  int orio_i;
  int t,i,j;

  for (orio_i=0; orio_i<REPS; orio_i++)
  {
    orio_t_start = rtclock();
    
      
    for (t=1; t<=T-1; t++) 
      for (i=1; i<=N-2; i++) 
	a[t][i] = 0.333 * (a[t-1][i-1] + a[t-1][i] + a[t-1][i+1]);


    orio_t_end = rtclock();
    orio_t_total += orio_t_end - orio_t_start;
  }
  
  orio_t_total = orio_t_total / REPS;
  printf("%f\n", orio_t_total);
  
  return 0;
}
                                    
