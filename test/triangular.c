#pragma scop
for (i = 0; i < N; i++) {
  for (j = 0; j < i; j++) {
    a[i] += a[j];
  }
}

#pragma endscop
