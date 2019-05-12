
int u[2][1003][1003];

int main() {
  int N = 1000, T = 1000;
  int i, j, t;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 4; i <= N - 1; i++) {
      for (j = 4; j <= N - 1; j++) {
        /* function doesn't represent an accurate discretization, just
         * representative */
        u[t % 2][i][j] = u[(t - 1) % 2][i - 1][j] + u[(t - 1) % 2][i][j + 1] +
                         u[(t - 1) % 2][i + 1][j] + u[(t - 1) % 2][i][j - 1] +
                         u[(t - 1) % 2][i - 4][j] + u[(t - 1) % 2][i + 4][j] +
                         u[(t - 1) % 2][i - 3][j] + u[(t - 1) % 2][i + 3][j] +
                         u[(t - 1) % 2][i - 2][j] + u[(t - 1) % 2][i + 2][j] +
                         u[(t - 1) % 2][i][j - 4] + u[(t - 1) % 2][i][j + 4] +
                         u[(t - 1) % 2][i][j - 3] + u[(t - 1) % 2][i][j + 3] +
                         u[(t - 1) % 2][i][j - 2] + u[(t - 1) % 2][i][j + 2] +
                         u[(t - 1) % 2][i][j];
      }
    }
  }
#pragma endscop

  return (int)u[T - 1][1];
}
