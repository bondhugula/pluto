// CHECK: T(S1): (t-j, t+j, t)
// CHECK: loop types (loop, loop, loop)
// CHECK: T(S2): (t-j, t+j, t+i)
// CHECK: loop types (loop, loop, loop)
// CHECK: T(S3): (t-j, t+j, t+i)
// CHECK: loop types (loop, loop, loop)
// CHECK: T(S4): (t-j, t+j+1, t+i+1)
// CHECK: loop types (loop, loop, loop)
//
// CHECK: Output written

/* Tests the greedy cost model with clustering in the dfp approach. The greedy
 * cost model finds the right permutation that enables fusion of all statements,
 * and thus after skewing, diamond tiling is also enabled. The hyperplane types
 * are also checked because after skewing, the hyperplane type changes from
 * scalar to loop at the last level for the first statement. */
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
