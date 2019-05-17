#define N 10

int main(void) {
  int i, j;
  int **a, **b;

#pragma scop
  for (i = 0; i < 3 * N - 1; i++)
    for (j = 0; j < N; j++) {
      if ((i + j >= N - 1) && (i + j <= 3 * N - 2))
        a[i][j] = 0;
      if ((i + j >= 2 * N - 1) && (i + j <= 4 * N - 2))
        b[i][j] = a[i - N][j];
    }
#pragma endscop
  return 0;
}
