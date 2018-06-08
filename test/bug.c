#pragma scop
for (i = 1; i < N-1; i++) {
    A[i] = A[i] * 0.3333333333;
}
#pragma endscop
