// CHECK: Output written

#pragma scop
for (i = 0; i < 100; i++) {
  for (j = 0; j < 4; j++) {
    A[4 * i + j] = 2 * A[4 * i + j] + 2;
  }
}
#pragma endscop
