// CHECK: Output written

#pragma scop
for (i = 0; i < N; i++) {
  x[i] = a[i] + b[i];
}
s = 0.0;
for (i = 0; i < N; i++) {
  s = s + x[i] * 3;
}
#pragma endscop
