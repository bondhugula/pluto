#pragma scop
for (i=1; i<N; i++)
    for (j=1; j<N; j++)
        a[i][j] = a[i-1][j] + a[i][j-1];
#pragma endscop
