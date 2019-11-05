#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#ifdef OPENBLAS
#include "cblas.h"
#elif BLIS
#include "blis/blis.h"
#elif MKL
#include "mkl.h"
#endif

#ifndef NUM_REPS
#define NUM_REPS 1
#endif

#define M 2048
#define N 2048
#define K 2048

#define alpha 1
#define beta 1

double A[M][K];
double B[K][N];
double C[M][N];

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_matrices() {
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      A[i][j] = (i + j);
      B[i][j] = (double)(i * j);
      C[i][j] = 0.0;
    }
  }
}

void print_matrix() {
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
  int i, j;
  int LDA, LDB, LDC;
  double *a, *b, *c;

  LDA = M;
  LDB = N;
  LDC = M;

  init_matrices();

  IF_TIME(t_start = rtclock());

  double _alpha = alpha;
  double _beta = beta;

  for (int t = 0; t < NUM_REPS; ++t) {
#ifdef BLIS
    bli_dgemm(BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE, M, N, K, &_alpha, &A[0][0],
              N, 1, &B[0][0], N, 1, &_beta, &C[0][0], N, 1);
#else
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, alpha,
                &A[0][0], LDA, &B[0][0], LDB, beta, &C[0][0], LDC);
#endif
  }

  IF_TIME(t_end = rtclock());
  IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));
  IF_TIME(fprintf(stdout, "%0.2lf GFLOPS\n",
                  2.0 * NUM_REPS * M * N * K / (t_end - t_start) / 1E9));

  if (fopen(".test", "r")) {
    print_matrix();
  }

  return 0;
}
