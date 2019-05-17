// CHECK: Output written
#pragma scop
for (i = 0; i <= -M; i++) {
  s += 1;
}
#pragma endscop
