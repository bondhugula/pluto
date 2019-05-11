// CHECK: T(S1): (t, t+i, 2t+i+j)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#define N 1000
#define T 500
double a[N][N];

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
      a[i][j] = i * i + j * j;
    }
  }
}

void print_array() {
  int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      fprintf(stderr, "%0.2lf ", a[i][j]);
      if (j % 80 == 20)
        fprintf(stderr, "\n");
    }
  }
  fprintf(stderr, "\n");
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

int main() {
  int i, j, k, t;

  double t_start, t_end;

  init_array();

#ifdef PERFCTR
  PERF_INIT;
#endif

  IF_TIME(t_start = rtclock());

#pragma scop
  for (t = 0; t <= T - 1; t++) {
    for (i = 1; i <= N - 2; i++) {
      for (j = 1; j <= N - 2; j++) {
        a[i][j] = (a[i - 1][j - 1] + a[i - 1][j] + a[i - 1][j + 1] +
                   a[i][j - 1] + a[i][j] + a[i][j + 1] + a[i + 1][j - 1] +
                   a[i + 1][j] + a[i + 1][j + 1]) /
                  9.0;
      }
    }
  }
#pragma endscop

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
