
#pragma scop
/* pluto start (N) */

for (x = 0; x < 100; x++) {
  A[2 * (x)] = 1;
}

for (x = 0; x < 100; x++) {
  B[2 * (x)] = A[2 * (x) + 2];
}
/* pluto end */

#pragma endscop
