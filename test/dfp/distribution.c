/* Checks for complete loop distribution when there is no dependence between two
 * SCCs. */
// CHECK: T(S1): (1, i, j, k)
// CHECK: T(S2): (0, i, j, k)
// [Pluto] After tiling:
// CHECK: T(S1): (1, i/32, j/32, k/32, i, j, k)
// CHECK: T(S2): (0, i/32, j/32, k/32, i, j, k)
// [Pluto] After intra-tile optimize
// CHECK: T(S1): (1, i/32, j/32, k/32, i, k, j)
// CHECK: T(S2): (0, i/32, j/32, k/32, i, k, j)
// CHECK: Output written

#pragma scop
/* E := A*B */
for (i = 0; i < _PB_NI; i++)
  for (j = 0; j < _PB_NJ; j++) {
    for (k = 0; k < _PB_NK; ++k)
      E[i][j] += A[i][k] * B[k][j];
  }
/* F := C*D */
for (i = 0; i < _PB_NJ; i++)
  for (j = 0; j < _PB_NL; j++) {
    for (k = 0; k < _PB_NM; ++k)
      F[i][j] += C[i][k] * D[k][j];
  }
#pragma endscop
}
