// CHECK: Output written

#pragma scop
for (x = 0; x < 100; x++) {
  for (y = 0; y < 100; y++) {
    A[x][y] = 1;
  }
}

for (x = 0; x < 100; x++) {
  for (y = 0; y < 100; y++) {
    B[x][y] = A[x + 1][y] + A[x][y + 1];
  }
}
#pragma endscop
