#pragma scop
  for (i = 1; i <= n; i++) {
    for (j = 1; j <= m; j++) {
      data[i][j] -= mean[j];
	}
  }

  for (j1 = 1; j1 <= m; j1++) {
    for (j2 = j1; j2 <= m; j2++) {
      for (i = 1; i <= n; i++)
        symmat[j1][j2] += data[i][j1] * data[i][j2];
    }
  }
#pragma endscop
