int u[100][1000][1000];

int main() {
  int N = 1000, T = 100;
  int i, j, t;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 1; i <= N - 1; i++) {
      for (j = 1; j <= N - 1; j++) {
        u[t][i][j] = u[(t - 1)][i - 1][j] + u[(t - 1)][i][j + 1] +
                     u[(t - 1)][i + 1][j] + u[(t - 1)][i][j - 1];
      }
    }
  }
#pragma endscop

  return (int)u[T - 1][1][1];
}
