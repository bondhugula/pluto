
int main() {
  int i, j;
  long N = 1000;

  int a[1000][1000];
#pragma scop
  for (i = 1; i < N; i++) {
    for (j = 0; j < N; j++) {
      a[i][j] += a[i - 1][j] + a[i - 1][N - 1 - j];
    }
  }
#pragma endscop

  return 0;
}
