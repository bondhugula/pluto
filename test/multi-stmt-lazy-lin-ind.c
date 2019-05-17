/* A stripped-down version of a kernel obtained from Tomofomo Yuki */
void kernel(int M, int N, int D[M][N]) {
  int i, j;
  int Dp[N][M];

#pragma scop
  for (i = 1; i < M; i++) {
    for (j = 1; j < N; j++) {
      D[i][j] = D[i - 1][j - 1];
      Dp[i][j] = D[j][i];
    }
  }
#pragma endscop
}
