//
// CHECK: T(S1): (t, t+j, t)
// CHECK: T(S2): (t, t+j, t+i)
// CHECK: T(S3): (t, t+j, t+i)
// CHECK: T(S4): (t, t+j+1, t+i+1)
//
// After intra-tile optimize
// CHECK: T(S1): (t, t, t+j)
// CHECK: T(S2): (t, t+i, t+j)
// CHECK: T(S3): (t, t+i, t+j)
// CHECK: T(S4): (t, t+i+1, t+j+1)
// CHECK: Output written
//
#define tmax 128
#define nx 2048
#define ny 2048

double ex[nx][ny + 1];
double ey[nx + 1][ny];
double hz[nx][ny];

int main() {
  int t, i, j;

#pragma scop
  for (t = 0; t < tmax; t++) {
    for (j = 0; j < ny; j++)
      ey[0][j] = t;
    for (i = 1; i < nx; i++)
      for (j = 0; j < ny; j++)
        ey[i][j] = ey[i][j] - 0.5 * (hz[i][j] - hz[i - 1][j]);
    for (i = 0; i < nx; i++)
      for (j = 1; j < ny; j++)
        ex[i][j] = ex[i][j] - 0.5 * (hz[i][j] - hz[i][j - 1]);
    for (i = 0; i < nx; i++)
      for (j = 0; j < ny; j++)
        hz[i][j] = hz[i][j] -
                   0.7 * (ex[i][j + 1] - ex[i][j] + ey[i + 1][j] - ey[i][j]);
  }
#pragma endscop

  return 0;
}
