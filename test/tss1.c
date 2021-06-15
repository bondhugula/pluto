// CHECK: T(S1): (i, j, k)
// TILE-PARALLEL: T(S1): (0, i/32, j/512, i, j, 0, 0)
// TILE-PARALLEL: T(S2): (1, i/32, j/512, k/256, i, k, j)

#define M 2048
#define N 2048
#define K 2048

double A[M][K + 13];
double B[K][N + 13];
double C[M][N + 13];

int main() {
  int i, j, k;
  register double s;

#pragma scop
  for (i = 0; i < M; i++)
    for (j = 0; j < N; j++) {
      C[i][j] = 0.0f;
      for (k = 0; k < K; k++)
        C[i][j] = C[i][j] + A[i][k] * B[k][j];
    }
#pragma endscop

  return 0;
}
