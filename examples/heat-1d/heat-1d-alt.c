/*
 * Discretized 1D heat equation stencil with non periodic boundary conditions
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#ifdef PERFCTR
#include "papiStdEventDefs.h"
#include <papi.h>
#include "papi_defs.h"
#endif

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

/*
 * N is the number of points
 * T is the number of timesteps
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 5000
#define T 500
#define BASE 1024
#endif

#define NUM_FP_OPS 4

/* Define our arrays */
#pragma declarations
double point[N+2];
double new_point[N+2];
#pragma enddeclarations

#define __PLACE_TO_INSERT_FORWARD_DECLARATIONS

/*
 * Print the array
 */
void print_points() {
  int i;

  for(i = 0; i < N+2; i++) {
    fprintf(stderr, "%e ", point[i]);
  }
  fprintf(stderr, "\n");
}

int main(int argc, char * argv[]) {
  long int t, i, j, k;
  //  double point[2][N+2];// declaration here overflows the activation record stack. Declare globally

  double t_start = 0.0, t_end = 0.0;

  double Total = 0.0;

#ifdef DEBUG
  printf("Number of points = %ld\nNumber of timesteps = %ld\n", N, T);
#endif

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 1; i < N+1; i++) {
    point[i] = 1.0 * (rand() % BASE);
  }
  point[0] = 0.0;
  point[N+1] = 0.0;

  IF_TIME(t_start = rtclock());

#pragma scop
  for (t = 0; t < T-1; t++) {
    for (i = 1; i < N+1; i++) {
      new_point[i] = 0.125 * (point[i+1] - 2.0 * point[i] + point[i-1]);
    }
    for (i = 1; i < N+1; i++) {
        point[i] = new_point[i];
    }
  }
#pragma endscop

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "time = %0.6lfs\n", t_end - t_start));

#ifdef __MPI
  if (my_rank == 0) {
#endif
    Total = 0.0;
    for(i = 0; i < N+2; i++) {
      Total += point[i];
    }
    printf("Sum(A): %e\n", Total);
#ifdef __MPI
  }
#endif

  if (fopen(".test", "r")) {
#ifdef __MPI
    if (my_rank == 0) {
      print_points();
    }
#else
    print_points();
#endif
  }

  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
