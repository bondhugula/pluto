
int u[2][1000];

#define N 500
#define T 1000

int main() {
  int i, t;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 0; i <= 2 * N - 1; i++) {
      u[t % 2][i] = u[(t - 1) % 2][i == 0 ? 2 * N - 1 : i - 1] +
                    u[(t - 1) % 2][i] +
                    u[(t - 1) % 2][i == 2 * N - 1 ? 0 : i + 1];
    }
  }
#pragma endscop

  return (int)u[T - 1][1];
}
