#pragma scop
for(i=0; i<M; i++)
    for(j=0; j<N; j++)
        for(k=0; k<K; k++)
            C[i][j] = _C_[i][j] + A[i][k] * B[k][j];
#pragma endscop
