
#pragma scop
for (i = 2; i <= N; i++) {
  a[i] = b[i] + c[i];
  b[i + 2] = a[i - 1] + c[i - 1];
  a[i + 1] = b[i + 3] + 1;
}
#pragma endscop
