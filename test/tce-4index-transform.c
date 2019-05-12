// CHECK: Output written
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#include "decls.h"

#include "util.h"

int main() {
  int a, q, r, s, p, b, c, d;

  init_array();

#pragma scop
  /* pluto start (N) */
  for (a = 0; a < N; a++)
    for (q = 0; q < N; q++)
      for (r = 0; r < N; r++)
        for (s = 0; s < N; s++)
          for (p = 0; p < N; p++)
            T1[a][q][r][s] = T1[a][q][r][s] + A[p][q][r][s] * C4[p][a];

  for (a = 0; a < N; a++)
    for (b = 0; b < N; b++)
      for (r = 0; r < N; r++)
        for (s = 0; s < N; s++)
          for (q = 0; q < N; q++)
            T2[a][b][r][s] = T2[a][b][r][s] + T1[a][q][r][s] * C3[q][b];

  for (a = 0; a < N; a++)
    for (b = 0; b < N; b++)
      for (c = 0; c < N; c++)
        for (s = 0; s < N; s++)
          for (r = 0; r < N; r++)
            T3[a][b][c][s] = T3[a][b][c][s] + T2[a][b][r][s] * C2[r][c];

  for (a = 0; a < N; a++)
    for (b = 0; b < N; b++)
      for (c = 0; c < N; c++)
        for (d = 0; d < N; d++)
          for (s = 0; s < N; s++)
            B[a][b][c][d] = B[a][b][c][d] + T3[a][b][c][s] * C1[s][d];
            /* pluto end */
#pragma endscop

#ifdef TEST
  print_array();
#endif
  return 0;
}
