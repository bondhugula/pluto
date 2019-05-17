#pragma scop
for (i = 0; i < N; i += 2) {
  a[i] = b[i];
}
#pragma endscop
