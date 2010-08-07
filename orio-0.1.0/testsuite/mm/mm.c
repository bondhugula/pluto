for (i=0; i<=M-1; i++ )
  for (j=0; j<=N-1; j++ )
    for (k=0; k<=K-1; k++ )
      C[i][j]=C[i][j]+A[i][k]*B[k][j];
