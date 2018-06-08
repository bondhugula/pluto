
for (t=0; t<=T-1; t++)
    for (i=1; i<=N-2; i++)
        for (j=1; j<=N-2; j++)
            A[i][j] = (A[i-1][j-1] + A[i-1][j] + A[i-1][j+1]
                       + A[i][j-1] + A[i][j] + A[i][j+1]
                       + A[i+1][j-1] + A[i+1][j] + A[i+1][j+1])/9.0;
