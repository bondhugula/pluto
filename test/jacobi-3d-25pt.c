// CHECK: T(S1): (t, 4t+i, 4t+j, 4t+k)
// TILE-PARALLEL: T(S1): (4t-i, 4t+i, 4t+j, 4t+k)

int u[2][1003][1003][1003];

int main() {
  int N = 1000, T = 1000;
  int i, j, t, k;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 4; i <= N - 1; i++) {
      for (j = 4; j <= N - 1; j++) {
        for (k = 4; k <= N - 1; k++) {
          /* function doesn't represent an accurate discretization, just
           * representative */
          u[t % 2][i][j][k] =
              u[(t - 1) % 2][i - 4][j][k] + u[(t - 1) % 2][i + 4][j][k] +
              u[(t - 1) % 2][i - 3][j][k] + u[(t - 1) % 2][i + 3][j][k] +
              u[(t - 1) % 2][i - 2][j][k] + u[(t - 1) % 2][i + 2][j][k] +
              u[(t - 1) % 2][i - 1][j][k] + u[(t - 1) % 2][i + 1][j][k] +
              u[(t - 1) % 2][i][j - 4][k] + u[(t - 1) % 2][i][j + 4][k] +
              u[(t - 1) % 2][i][j - 3][k] + u[(t - 1) % 2][i][j + 3][k] +
              u[(t - 1) % 2][i][j - 2][k] + u[(t - 1) % 2][i][j + 2][k] +
              u[(t - 1) % 2][i][j - 1][k] + u[(t - 1) % 2][i][j + 1][k] +
              u[(t - 1) % 2][i][j][k - 4] + u[(t - 1) % 2][i][j][k + 4] +
              u[(t - 1) % 2][i][j][k - 3] + u[(t - 1) % 2][i][j][k + 3] +
              u[(t - 1) % 2][i][j][k - 2] + u[(t - 1) % 2][i][j][k + 2] +
              u[(t - 1) % 2][i][j][k - 1] + u[(t - 1) % 2][i][j][k + 1] +
              u[(t - 1) % 2][i][j][k];
        }
      }
    }
  }
#pragma endscop

  return (int)u[T - 1][1][1][1];
}
