#pragma scop
for (i=0; i<N; i++) {
    for (j=0; j<M; j++) {
        s += a[i][j];
    }
}
#pragma endscop
