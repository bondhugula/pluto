#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

#ifdef TIME
double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
#endif

#define NMAX 3000
#define MEASURE_TIME 1

static double a[NMAX][NMAX], c[NMAX][NMAX];

void dsyrk(long N) 
{
  int i,j,k;

#pragma scop
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=j; k<N; k++) {
        c[j][k] += a[i][j] * a[i][k];
      }
    }
  }
#pragma endscop
}


int main()
{
  double t_start, t_end;

  long N=NMAX;
  int i,j;

  for (i = 0; i < NMAX; i++) {
    for (j = 0; j < NMAX; j++) {
      c[i][j] = 0.0;
      a[i][j] = i*j*3.2345 / NMAX;
    }
  }

  IF_TIME(t_start = rtclock());

  dsyrk(N);

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

  if (fopen(".test", "r"))  {
    for (i = 0; i < NMAX; i++) {
      for (j = 0; j < NMAX; j++) {
        printf("%lf ", c[i][j]);
      }
    }
  }

  return 0;

}
