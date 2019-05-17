
#pragma scop
for (t = 0; t < T; t++) {
  for (i = 2; i < N - 1; i++) {
    b[i] = log10(0.33333 * (a[i - 1] + cos(a[i]) + a[i + 1]));
  }
  for (j = 2; j < N - 1; j++) {
    a[j] = b[j];
  }
}
#pragma endscop
