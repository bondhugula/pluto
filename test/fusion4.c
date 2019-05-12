// CHECK: Output written

#pragma scop
for (x = 0; x < 10000; x++) {
  for (y = 0; y < 10000; y++) {
    B[x][y] = x + y;
  }
}

for (x = 0; x < 10000; x++) {
  for (y = 0; y < 10000; y++) {
    C[x][y] = B[x - 1][y - 1] + B[x][y - 1];
  }
}
#pragma endscop
