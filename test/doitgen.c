// CHECK: Output written

#define N 25

double A[N][N][N];
double sum[N][N][N];
double C4[N][N];

int main() {
  int p, q, r, s;

#pragma scop
  for (r = 0; r < N; r++) {
    for (q = 0; q < N; q++) {
      for (p = 0; p < N; p++) {
        sum[r][q][p] = 0;
        for (s = 0; s < N; s++) {
          sum[r][q][p] = sum[r][q][p] + A[r][q][s] * C4[s][p];
        }
      }
      for (p = 0; p < N; p++) {
        A[r][q][p] = sum[r][q][p];
      }
    }
  }
#pragma endscop

  return 0;
}
