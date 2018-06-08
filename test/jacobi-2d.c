
#pragma scop
for(t=1; t<=T-1; t++) {
    for(i=1; i<=N-1; i++) {
        for(j=1; j<=N-1; j++) {
            u[t][i][j] = u[t-1][i-1][j] + u[t-1][i][j+1] + u[t-1][i+1][j] + u[t-1][i][j-1];
        }
    }
}
#pragma endscop
