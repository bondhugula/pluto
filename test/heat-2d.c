// CHECK: T(S1): (t, t+i, t+j)
// TILE-PARALLEL: T(S1): ((t-i)/32+(t+i)/32, (t+i)/32, (t+j)/32, t, t+i, t+j)
#define N 4000L
#define T 1000L

/* Define our arrays */
double A[2][N + 2][N + 2];

int main(int argc, char *argv[]) {
#pragma scop
  for (int t = 0; t < T; t++) {
    for (int i = 1; i < N + 1; i++) {
      for (int j = 1; j < N + 1; j++) {
        A[(t + 1) % 2][i][j] =
            0.125 * (A[t % 2][i + 1][j] - 2.0 * A[t % 2][i][j] +
                     A[t % 2][i - 1][j]) +
            0.125 * (A[t % 2][i][j + 1] - 2.0 * A[t % 2][i][j] +
                     A[t % 2][i][j - 1]) +
            A[t % 2][i][j];
      }
    }
  }
#pragma endscop

  return 0;
}
