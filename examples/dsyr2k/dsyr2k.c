#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "decls.h"

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


void dsyr2k(long N) {
	int i,j,k;
	
#pragma scop
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			for (k=j; k<N; k++) {
				c[j][k] += a[i][j] * b[i][k] + b[i][j] * a[i][k];
			}
		}
	}
#pragma endscop
}

double t_start, t_end;

int main()
{
  long N=NMAX;
  int i,j;

  for (i = 0; i < NMAX; i++) {
    for (j = 0; j < NMAX; j++) {
      c[i][j] = 0.0;
      a[i][j] = 1.0+i*j*32.5/(N*N);
      b[i][j] = 1.0+i*j*3.4565/(N*N);
    }
  }

  IF_TIME(t_start = rtclock());

  dsyr2k(N);

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

  if (fopen(".test", "r")) {
    for (i = 0; i < NMAX; i++) {
      for (j = 0; j < NMAX; j++) {
        fprintf(stdout, "%lf ", c[i][j]);
      }
    }
        fprintf(stdout, "\n");
  }

  return 0;

}
