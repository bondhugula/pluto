#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main() {
  int i, j, k, l, m, n, t;

  int i1, i2;

#pragma scop
  for (t = 0; t < T; t++) {

    for (i1 = 0; i1 < N; i1++) {
      for (i2 = 1; i2 < N; i2++) {
        X[i1][i2] = X[i1][i2] - X[i1][i2 - 1] * A[i1][i2] / B[i1][i2 - 1];
        B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2 - 1];
      }
    }

    for (i1 = 0; i1 < N; i1++) {
      X[i1][N - 1] = X[i1][N - 1] / B[i1][N - 1];
    }

    for (i1 = 0; i1 < N; i1++) {
      for (i2 = 0; i2 <= N - 2; i2++) {
        X[i1][N - i2 - 2] =
            (X[i1][N - 2 - i2] - X[i1][N - 2 - i2 - 1] * A[i1][N - i2 - 3]) /
            B[i1][N - 3 - i2];
      }
    }

    for (i1 = 1; i1 < N; i1++) {
      for (i2 = 0; i2 < N; i2++) {
        X[i1][i2] = X[i1][i2] - X[i1 - 1][i2] * A[i1][i2] / B[i1 - 1][i2];
        B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1 - 1][i2];
      }
    }

    for (i2 = 0; i2 < N; i2++) {
      X[N - 1][i2] = X[N - 1][i2] / B[N - 1][i2];
    }

    for (i1 = 0; i1 <= N - 2; i1++) {
      for (i2 = 0; i2 < N; i2++) {
        X[N - 2 - i1][i2] =
            (X[N - 2 - i1][i2] - X[N - i1 - 3][i2] * A[N - 3 - i1][i2]) /
            B[N - 2 - i1][i2];
      }
    }
  }
#pragma endscop

  return 0;
}
