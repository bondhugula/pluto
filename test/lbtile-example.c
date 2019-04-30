/*  Thanks to Tobias Grosser */
#pragma scop
for (t = 5; t < T-1; t++) {
    for (i = 5; i < N-5; i++) {
        A[t+1][i] =   0.125 * (A[t][i-1] + A[t-2][i+1]);
    }
} 
#pragma endscop
