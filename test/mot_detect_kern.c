#include <stdio.h>

#define m 4
#define n 4
#define M 32
#define N 32

// motion detection kernel
int main(int argc, char *argv[]) {
  int i, j, k, l;
  int optDelta[(M - m + 1) * (N - n + 1) + 1];
  int Delta[M - m + 1][N - n + 1][(2 * m + 1) * (2 * n + 1)];
  int A[(m + 1) * (m + 1)][(n + 1) * (n + 1)]; // input
  int opt[1];                                  // output

  optDelta[0] = 0;
  for (i = m; i <= M; i++) {
    for (j = n; j <= N; j++) {
      Delta[i][j][0] = 0;
      for (k = i - m; k <= i + m; k++) {
        for (l = j - n; l <= j + n; l++) {
          Delta[i][j][(2 * n + 1) * k - (2 * n + 1) * i + l - j +
                      (2 * m * n + m + n + 1)] =
              A[i][j] - A[k][l] +
              Delta[i][j][(2 * n + 1) * k - (2 * n + 1) * i + l - j +
                          (2 * m * n + m + n)];
        }
      }
      optDelta[(N - n + 1) * i + j - (m * N - m * n + m + n + 1)] =
          Delta[i][j][(2 * m + 1) * (2 * n + 1)] +
          optDelta[(N - n + 1) * i + j - (m * N - m * n + m + n)];
    }
  }
  opt[0] = optDelta[(M - m + 1) * (N - n + 1)];
}
