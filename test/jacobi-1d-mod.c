
int u[2][1000];

int main() {
  int N = 1000, T = 1000;
  int i, t;

#pragma scop
  for (t = 1; t <= T - 1; t++) {
    for (i = 1; i <= N - 2; i++) {
      u[t % 2][i] =
          u[(t - 1) % 2][i - 1] + u[(t - 1) % 2][i] + u[(t - 1) % 2][i + 1];
    }
  }
#pragma endscop

  return (int)u[T - 1][1];
}
