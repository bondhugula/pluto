// CHECK: T(S1): (j, i)
#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < M; j++) {
    a[j][i] = a[j][i] + 1;
  }
}
#pragma endscop
