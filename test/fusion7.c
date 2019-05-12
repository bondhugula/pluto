// CHECK: Output written
#pragma scop
for (i = 0; i < N; i++) {
  a[i] = b[i];
}
for (i = 0; i < M; i++) {
  a[i] = b[i];
}

#pragma endscop
