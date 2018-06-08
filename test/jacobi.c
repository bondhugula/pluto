
#pragma scop
for(t=1; t<=T-1; t++) {
    for(i=1; i<=N-1; i++) {
        u[t][i] = 0.333*(u[t-1][i-1] + u[t-1][i] + u[t-1][i+1]);
    }
}
#pragma endscop
