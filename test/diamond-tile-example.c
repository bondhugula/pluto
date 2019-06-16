// Example - thanks to Tobias Grosser.
// Dependences are (1, 1) and (3, -1). The diamond tiling hyperplanes are (1,1)
// and (1,-1).
//
// TILE-PARALLEL: T(S1): (t+i, t-i)
// TILE-PARALLEL: After tiling:
// TILE-PARALLEL: T(S1): ((t+i)/32, (t-i)/32, t+i, t-i)
#pragma scop
for (t = 5; t < T - 1; t++) {
  for (i = 5; i < N - 5; i++) {
    A[t + 1][i] = 0.125 * (A[t][i - 1] + A[t - 2][i + 1]);
  }
}
#pragma endscop
