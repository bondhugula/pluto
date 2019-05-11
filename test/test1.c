#pragma scop

for (x = 0; x < 100; x++) {
  for (z = 0; z < 4; z++) {
    A[4 * x + z] = 1;
  }
}

for (x = 0; x < 100; x++) {
  for (z = 0; z < 4; z++) {
    B[4 * x + z] = A[4 * x + z];
  }
}
#pragma endscop

write();
