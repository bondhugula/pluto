// CHECK: 2i-_N = 0
// CHECK: 2j-_N = 0
// CHECK: [iss] Splitting S1 into 4 statements
// CHECK-POST-TILE: [pluto] After post-tile distribution
// CHECK-POST-TILE: T(S1): ((t-i)/32, (t+i)/32, (t+j)/32, t, 1, t+i, t+j)
// CHECK-POST-TILE: T(S2): ((t+i-_N)/32, (t-i+_N)/32, (t+j)/32, t, 2, t-i+_N, t+j)
// CHECK-POST-TILE: T(S3): ((t-i)/32, (t+i)/32, (t-j+_N)/32, t, 3, t+i, t-j+_N)
// CHECK-POST-TILE: T(S4): ((t+i-_N)/32, (t-i+_N)/32, (t-j+_N)/32, t, 4, t-i+_N, t-j+_N)

#define N 1600L
#define T 500L

/* Define our arrays */
double A[2][N][N];

int main(int argc, char *argv[]) {
  short _N = N - 1;
  int _T = T;

#pragma scop
  for (int t = 0; t < _T; t++) {
    for (int i = 0; i <= _N; i++) {
      for (int j = 0; j <= _N; j++) {
        A[(t + 1) % 2][i][j] =
            0.125 * (A[t % 2][i == _N ? 0 : i + 1][j] - 2.0 * A[t % 2][i][j] +
                     A[t % 2][(i == 0) ? _N : (i - 1)][j]) +
            0.125 * (A[t % 2][i][j == _N ? 0 : j + 1] - 2.0 * A[t % 2][i][j] +
                     A[t % 2][i][(j == 0) ? _N : (j - 1)]) +
            A[t % 2][i][j];
      }
    }
  }
#pragma endscop
  return 0;
}
