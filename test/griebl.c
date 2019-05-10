#pragma scop
for (i = 0; i < n; i++)
A[2 * i] = 42

    for (j = 0; j < n; j++) B[j] = A[3 * j]
#pragma endscop
