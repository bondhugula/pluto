#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define TIME 1

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


#define NMAX 2000

#define A_SIZE NMAX
#define B_SIZE NMAX
static double A[A_SIZE][A_SIZE];
static double B[B_SIZE][B_SIZE];
static double C[B_SIZE][B_SIZE];

void tmm(long Ni, long Nj, long Nk) {
	int i, j, k;
	
#pragma scop
	for(i = 0; i < Ni; i++) {
		for(j=i; j < Nj; j++) { 
			for(k=i;  k< Nk; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
#pragma endscop

}

int main(int argc, char *argv)
{
  double t_start, t_end;
	long N=NMAX;

	int i, j, k;

	for (i = 0; i < NMAX; i++) {
		for (j = 0; j < NMAX; j++) {
			C[i][j] = 0.0;
			A[i][j] = B[i][j] = (i+j)/2.0;
		}
	}


  IF_TIME(t_start = rtclock());

	tmm(N,N,N);

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

  if (fopen(".test", "r")) {
    for (i = 0; i < NMAX; i++) {
      for (j = 0; j < NMAX; j++) {
        fprintf(stdout, "%lf ", C[i][j]);
      }
    }

  }

  return 0;
}
