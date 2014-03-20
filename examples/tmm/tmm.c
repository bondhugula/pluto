#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>

#define NMAX 2000

#define A_SIZE NMAX
#define B_SIZE NMAX

#pragma declarations

double A[A_SIZE][A_SIZE];
double B[B_SIZE][B_SIZE];
double C[B_SIZE][B_SIZE];

#pragma enddeclarations

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


void print_array()
{
    int i, j;

    for (i = 0; i < NMAX; i++) {
      for (j = i; j < NMAX; j++) {
        fprintf(stderr, "%lf ", C[i][j]);
      }
      fprintf(stderr, "\n");
    }
}


int main(int argc, char *argv)
{
    double t_start, t_end;

	int i, j, k;

	for (i = 0; i < NMAX; i++) {
		for (j = 0; j < NMAX; j++) {
			C[i][j] = 0.0;
			A[i][j] = B[i][j] = (i+j)/2.0;
		}
	}


  IF_TIME(t_start = rtclock());

#pragma scop
	for(i = 0; i < NMAX; i++) {
		for(j=i; j < NMAX; j++) { 
			for(k=i;  k< NMAX; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
#pragma endscop



  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));

  if (fopen(".test", "r")) {
#ifdef MPI
      if (my_rank == 0) {
          print_array();
      }
#else
          print_array();
#endif
  }

  return 0;
}
