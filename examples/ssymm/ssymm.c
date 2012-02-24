#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define NMAX 500

#pragma declarations
double a[NMAX][NMAX], b[NMAX][NMAX], c[NMAX][NMAX];
#pragma enddeclarations

#define TIME 1

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}


int main()
{
  long N=NMAX;
  int i,j, k;
  double t_start, t_end;

  if (fopen(".test", "r"))  {
    N = N/2;
  }

  for (i = 0; i < NMAX; i++) {
    for (j = 0; j < NMAX; j++) {
      c[i][j] = 0.0;
      a[i][j] = b[i][j] = i*j*0.5 / NMAX;
    }
  }

  IF_TIME(t_start = rtclock());

#pragma scop
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<j-1; k++) {
        c[i][k] += a[j][k] * b[i][j];
        c[i][j] += a[j][j] * b[i][j];
      }
      c[i][j] += a[j][j] * b[i][j];
    }
  }
#pragma endscop

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

  if (fopen(".test", "r"))  {
      for (i = 0; i < NMAX; i++) {
          for (j = 0; j < NMAX; j++) {
              printf("%lf ", c[i][j]);
          }
          printf("\n");
      }
  }
  return 0;
}
