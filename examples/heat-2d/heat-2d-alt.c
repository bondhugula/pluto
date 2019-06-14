/*
 * Discretized 2D heat equation stencil with non periodic boundary conditions
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 *
 * This is the alternate version used for distributed memory code generation,
 * since accesses with modulo's aren't supported with distmem.
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

#define GLOBAL_MY_RANK
int my_rank;

/*
 * N is the number of points
 * T is the number of timesteps
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 2048L
#define T 512L
#endif

#define NUM_FP_OPS 10

#define __DECLARATION_OF_point double point[N+2][N+2]
#define __DECLARATION_OF_POINTER_TO_point(point) double (*point)[N+2]
#define __DECLARATION_OF_new_point double new_point[N+2][N+2]
#define __DECLARATION_OF_POINTER_TO_new_point(new_point) double (*new_point)[N+2]

#define __PLACE_TO_INSERT_FORWARD_DECLARATIONS

double rtclock()
{
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday (&Tp, &Tzp);
  if (stat != 0) printf("Error return from gettimeofday: %d",stat);
  return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

void init_array(double point[N+2][N+2], double new_point[N+2][N+2]) {
  const int BASE = 1024;

  long int i, j;

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 0; i < N+2; i++) {
    for (j = 0; j < N+2; j++) {
      point[i][j] = 0.0;
      new_point[i][j] = 0.0;
    }
  }

  for (i = 1; i < N+1; i++) {
    for (j = 1; j < N+1; j++) {
      point[i][j] = 1.0 * (rand() % BASE);
    }
  }
}

/*
 * Print the array
 */
void print_array(double point[N+2][N+2]) {
  int i, j;

  for(i = 0; i < N+2; i++) {
    for(j = 0; j < N+2; j++) {
      fprintf(stderr, "%e ", point[i][j]);
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}

void kernel(double point[N+2][N+2], double new_point[N+2][N+2]) {
  long int t, i, j;
#pragma scop
  for (t = 0; t < T; t++) {
    for (i = 1; i < N+1; i++) {
      for (j = 1; j < N+1; j++) {
        new_point[i][j] =   0.125 * (point[i+1][j] - 2.0 * point[i][j] + point[i-1][j])
                            + 0.125 * (point[i][j+1] - 2.0 * point[i][j] + point[i][j-1])
                            + point[i][j];
      }
    }

    for(i = 1; i < N+1; i++) {
      for(j = 1; j < N+1; j++) {
        point[i][j] = new_point[i][j];
      }
    }
  }
#pragma endscop
}

int main(int argc, char * argv[]) {
  double t_start, t_end;

/* Define our arrays */
  double (*point)[N+2][N+2];
  point = (double(*)[N+2][N+2]) malloc((N+2)*(N+2)*sizeof(double));
  double (*new_point)[N+2][N+2];
  new_point = (double(*)[N+2][N+2]) malloc((N+2)*(N+2)*sizeof(double));

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  init_array(*point, *new_point);

#ifdef DEBUG
  printf("Number of points = %ld\nNumber of timesteps = %ld\n", N, T);
#endif

#ifdef PERFCTR
  PERF_INIT;
#endif

  IF_TIME(t_start = rtclock());

  kernel(*point, *new_point);

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "time = %0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
  PERF_EXIT;
#endif

  if (fopen(".test", "r")) {
#ifdef __MPI
    if (my_rank == 0) {
      print_array(*point);
    }
#else
    print_array(*point);
#endif
  }

  free(point);
  free(new_point);

  return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
