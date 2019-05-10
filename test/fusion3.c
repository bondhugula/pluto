// CHECK: Output written
#pragma scop
for (y = 0; y < 10000; y++) {
  for (x = 0; x < 10000; x++) {
    B[y][x] = x + y;
  }
}

for (x = 0; x < 10000; x++) {
  for (y = 0; y < 10000; y++) {
    C[y][x] = B[y - 1][x - 1] + B[y][x];
  }
}
#pragma endscop
