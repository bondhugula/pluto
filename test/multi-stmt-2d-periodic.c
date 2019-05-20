// CHECK: 2i-N = 0
// CHECK: 2j-N = 0
// CHECK: [iss] Splitting S1 into 2 statements
// CHECK: [iss] Splitting S2 into 2 statements
// CHECK: T(S1): (t, t+i, t+j)
// CHECK: T(S2): (t, t-i+N, t+j)
// CHECK: T(S3): (t, t+i, t-j+N)
// CHECK: T(S4): (t, t-i+N, t-j+N)
// CHECK: T(S5): (t, t+i, t+j)
// CHECK: T(S6): (t, t-i+N, t+j)
// CHECK: T(S7): (t, t+i, t-j+N)
// CHECK: T(S8): (t, t-i+N, t-j+N)
//
int u[2][1000][1000];
int v[2][1000][1000];

int main() {
  int N = 1000, T = 1000;
  int i, j, t;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 0; i <= N - 1; i++) {
      for (j = 0; j <= N - 1; j++) {
        u[t % 2][i][j] =
            u[(t - 1) % 2][i == 0 ? N - 1 : i - 1][j] + v[(t - 1) % 2][i][j] +
            u[(t - 1) % 2][i == N - 1 ? 0 : i + 1][j] +
            u[(t - 1) % 2][i == 0 ? N - 1 : i - 1][j == 0 ? N - 1 : j - 1] +
            u[(t - 1) % 2][i == N - 1 ? 0 : i + 1][j == N - 1 ? 0 : j + 1];
        v[t % 2][i][j] =
            v[(t - 1) % 2][i == 0 ? N - 1 : i - 1][j] + v[(t - 1) % 2][i][j] +
            v[(t - 1) % 2][i == N - 1 ? 0 : i + 1][j] +
            v[(t - 1) % 2][i == 0 ? N - 1 : i - 1][j == 0 ? N - 1 : j - 1] +
            u[(t - 1) % 2][i == N - 1 ? 0 : i + 1][j == N - 1 ? 0 : j + 1];
      }
    }
  }
#pragma endscop

  return (int)u[T - 1][1];
}
