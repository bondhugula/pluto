CONSTANT n;

for (i = 0; i < n; i++)
  for (k = 0; k < n; k++)
    for (j = 0; j < n; j++)
C[i, j] = C[i, j] +
          A[i, k] *
              B[k, j] for (i = 0; i < n; i++) for (k = 0; k < n;
                                                   k++) for (j = 0; j < n; j++)
                  D[i, j] = D[i, j] + E[i, k] * C[k, j]
