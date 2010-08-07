
/* pluto start (T,N) */
for (t=0; t<=T-1; t++)
  {

    for (i1=0; i1<=N-1; i1++)
      for (i2=1; i2<=N-1; i2++)
	{
          X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
          B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
	}

    for (i1=1; i1<=N-1; i1++)
      for (i2=0; i2<=N-1; i2++)
	{
          X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
          B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
	}
  }
/* pluto end */
