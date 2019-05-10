// CHECK: T(S1): (t, 2t+i, 0)
// CHECK: T(S2): (t, 2t+i+1, 1)
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#define N 2000000
#define T 1000

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double rtclock() {
  struct timezone Tzp;
  struct timeval Tp;
  int stat;
  stat = gettimeofday(&Tp, &Tzp);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

#pragma declarations
double a[N];
double b[N];
#pragma enddeclarations

double t_start, t_end;

int main() {
  int t, i, j;

  init_array();

  IF_TIME(t_start = rtclock());

#pragma scop
  for (t = 0; t < T; t++) {
    for (i = 2; i < N - 1; i++) {
      b[i] = 0.33333 * (a[i - 1] + a[i] + a[i + 1]);
    }
    for (i = 2; i < N - 1; i++) {
      a[i] = b[i];
    }
  }
#pragma endscop

  IF_TIME(t_end = rtclock());
  IF_TIME(printf("%0.6lfs\n", t_end - t_start));

#ifdef TEST
  print_array();
#endif

  return 0;
}
