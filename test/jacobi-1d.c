#define N 1000
#define T 1000

#include <stdio.h>

int main() {
  unsigned t, i;
  double A[T][N];
#pragma scop
  for (t = 1; t < T; t++) {
    for (i = 1; i < N - 1; i++) {
      A[t][i] = A[t - 1][i - 1] + A[t - 1][i] + A[t - 1][i + 1];
    }
  }
#pragma endscop
  printf("%lf\n", A[T-1][0]);
  return 0;
}
