// CHECK: Output written

#pragma scop
for (t = M; t <= N; t++) {
  a[i] = b[i];
}
#pragma endscop
