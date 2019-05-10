#pragma scop
for (i = 0; i < N; ++i)
  for (j = 0; j < N; ++j) {
    p1 = (matrix[i][j] == 1);
    irif(p1) for (k = 0; k < N; ++k) {
      p2 = (matrix[j][k] == 1);
      irif(p2) matrix[i][k] = 1;
    }
  }
#pragma endscop
