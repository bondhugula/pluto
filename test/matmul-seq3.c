

#pragma scop
// C = AB
for (i = 0; i < N; ++i)
for (j = 0; j < N; ++j)
for (k = 0; k < N; ++k)
C[i][j] += A[i][k] * B[k][j];
// F = DE
for (i = 0; i < N; ++i)
for (j = 0; j < N; ++j)
for (k = 0; k < N; ++k)
F[i][j] += D[i][k] * E[k][j];
// G = CF
for (i = 0; i < N; ++i)
for (j = 0; j < N; ++j)
for (k = 0; k < N; ++k)
G[i][j] += C[i][k] * F[k][j];
#pragma endscop
