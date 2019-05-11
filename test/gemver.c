// CHECK: Output written
#include <stdio.h>

#define N 4000

#define alpha 1
#define beta 1

double A[N][N];
double B[N][N];

double x[N];
double u1[N];
double u2[N];
double v2[N];
double v1[N];
double w[N];
double y[N];
double z[N];

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
    u1[i] = i;
    u2[i] = (i + 1) / N / 2.0;
    v1[i] = (i + 1) / N / 4.0;
    v2[i] = (i + 1) / N / 6.0;
    y[i] = (i + 1) / N / 8.0;
    z[i] = (i + 1) / N / 9.0;
    x[i] = 0.0;
    w[i] = 0.0;
    for (j = 0; j < N; j++) {
      A[i][j] = ((double)i * j) / N;
    }
  }
}

void print_array() {
  int i, j;

  for (i = 0; i < N; i++) {
    fprintf(stderr, "%0.2lf ", w[i]);
    if (i % 80 == 20)
      fprintf(stderr, "\n");
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

#define beta 1
#define alpha 1

main() {
  double t_start, t_end;
  int i, j;

  init_array();

  IF_TIME(t_start = rtclock());

  /* pluto start (N,alpha,beta) */
#pragma scop
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      x[i] = x[i] + beta * A[j][i] * y[j];

  for (i = 0; i < N; i++)
    x[i] = x[i] + z[i];

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      w[i] = w[i] + alpha * A[i][j] * x[j];
#pragma endscop

  /* pluto end */
  IF_TIME(t_end = rtclock());
  IF_TIME(printf("%0.6lfs\n", t_end - t_start));

#ifdef TEST
  print_array();
#endif

  return 0;
}
