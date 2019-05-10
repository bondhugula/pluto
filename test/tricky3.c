// CHECK: Output written

#pragma scop
for (p = 0; p < pointc; ++p) {
  dist_min = 0;
  for (k = 0; k < clusterc; ++k) {
    dist = 0;
    kmin = 0;
  }
  for (d = 0; d < dims; ++d) {
    clusterv = 0;
  }
}
#pragma endscop
