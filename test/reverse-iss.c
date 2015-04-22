#pragma scop
for (i=0; i<N; i++) {
    a[i] = 2*a[N-1-i];
}
#pragma endscop


