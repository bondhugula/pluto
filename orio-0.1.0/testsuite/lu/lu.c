
for (k=0; k<=N-1; k++) {
    for (j=k+1; j<=N-1; j++)
        A[k][j] = A[k][j]/A[k][k];
    for(i=k+1; i<=N-1; i++)
        for (j=k+1; j<=N-1; j++)
            A[i][j] = A[i][j] - A[i][k]*A[k][j];
}



