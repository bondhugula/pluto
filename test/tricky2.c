// CHECK: Output written

#pragma scop
for (n = 0; n < N; ++n) {
  autocorr += 1;
}
for (n = 0; n < M; ++n) {
  temp3 += 1;
}
#pragma endscop
