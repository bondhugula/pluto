// TILE-PARALLEL: T(S1): (i, j)
// TILE-PARALLEL: [Pluto] After tiling:
// TILE-PARALLEL: T(S1): (i/32, j/32, i, j)
// TILE-PARALLEL: [Pluto] After tile scheduling:
// TILE-PARALLEL: T(S1): (i/32+j/32, j/32, i, j)
#pragma scop

for (i = 1; i < N; i++) {
  for (j = 1; j < N; j++) {
    a[i][j] = a[i - 1][j] + a[i][j - 1];
  }
}
#pragma endscop
