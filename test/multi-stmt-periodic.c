
int u[2][1000];
int v[2][1000];
int w[2][1000];

// CHECK: 2i-N = 0
// CHECK: [iss] Splitting S1 into 2 statements
// CHECK: [iss] Splitting S2 into 2 statements
int main() {
  int N = 1000, T = 1000;
  int i, t;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 0; i <= N - 1; i++) {
      u[t % 2][i] = v[(t - 1) % 2][i == 0 ? N - 1 : i - 1] + u[(t - 1) % 2][i] +
                    u[(t - 1) % 2][i == N - 1 ? 0 : i + 1];
      v[t % 2][i] = v[(t - 1) % 2][i == 0 ? N - 1 : i - 1] + v[(t - 1) % 2][i] +
                    u[(t - 1) % 2][i == N - 1 ? 0 : i + 1];
    }
  }
#pragma endscop

  return (int)u[T - 1][1];
}
