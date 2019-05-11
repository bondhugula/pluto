#include <math.h>
#include <omp.h>
#define ceild(n, d) ceil(((double)(n)) / ((double)(d)))
#define floord(n, d) floor(((double)(n)) / ((double)(d)))
#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define M 2048
#define N 2048
#define K 2048
#define alpha 1
#define beta 1
double A[M][K + 13];
double B[K][N + 13];
double C[M][N + 13];

#ifdef PERFCTR
#include "papi_defs.h"
#include <papi.h>
#endif

#include <sys/time.h>
#include <unistd.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array() {
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i][j] = (i + j);
      B[i][j] = (double)(i * j);
      C[i][j] = 0.0;
    }
  }
}

void print_array() {
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      fprintf(stderr, "%lf ", C[i][j]);
      if (j % 80 == 79)
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
}

double rtclock() {
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday(&Tp, &Tzp);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}
double t_start, t_end;

int main() {
  int i, j, k;
  register double s;

  init_array();

#ifdef PERFCTR
  PERF_INIT;
#endif

  IF_TIME(t_start = rtclock());

  int t1, t2, t3, t4, t5, t6;

  register int lb, ub, lb1, ub1, lb2, ub2;
  register int lbv, ubv;

  omp_set_nested(1);
  omp_set_num_threads(2);
  /* Generated from PLUTO-produced CLooG file by CLooG v0.14.1 64 bits in 0.01s.
   */
  if ((M >= 1) && (N >= 1) && (K >= 1)) {
    lb1 = 0;
    ub1 = floord(M - 1, 32);
#pragma omp parallel for shared(lb1, ub1) private(lb2, ub2, t1, t2, t3, t4,    \
                                                  t5, t6)
    for (t1 = lb1; t1 <= ub1; t1++) {
      lb2 = 0;
      ub2 = floord(N - 1, 32);
#pragma omp parallel for shared(t1, lb1, ub1, lb2, ub2) private(t2, t3, t4,    \
                                                                t5, t6)
      for (t2 = lb2; t2 <= ub2; t2++) {
        for (t3 = 0; t3 <= floord(K - 1, 32); t3++) {
          for (t4 = max(0, 32 * t1); t4 <= min(M - 1, 32 * t1 + 31); t4++) {
            for (t5 = max(0, 32 * t3); t5 <= min(K - 1, 32 * t3 + 31); t5++) {
              {
                lbv = max(0, 32 * t2);
                ubv = min(N - 1, 32 * t2 + 31);
#pragma ivdep
#pragma vector always
                for (t6 = lbv; t6 <= ubv; t6++) {
                  C[t4][t6] = C[t4][t6] + A[t4][t5] * B[t5][t6];
                  ;
                }
              }
            }
          }
        }
      }
    }
  }
  /* End of CLooG code */

  IF_TIME(t_end = rtclock());
  IF_TIME(printf("%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
  PERF_EXIT;
#endif

#ifdef TEST
  print_array();
#endif
  return 0;
}
