#include <float.h>
#include <stdint.h>

void k_classify_local(float *pointv, uint32_t pointc, uint32_t point_lda,
                      float *clusterv, float *clusterv2, uint32_t *counterv,
                      uint32_t clusterc, uint32_t dims) {
  uint32_t p, k, d, kmin;
  float dist_min, dist, tmp;

  // A conservative (accesses) version

#pragma scop
  for (p = 0; p < pointc; ++p) {
    kmin = -1;
    dist_min = FLT_MAX;
    for (k = 0; k < clusterc; ++k) {
      dist = 0.0;
      for (d = 0; d < dims; ++d) {
        tmp = clusterv[0] - pointv[0];
        dist += tmp * tmp;
      }
      dist_min = (dist < dist_min) ? dist : dist_min;
      kmin = (dist < dist_min) ? k : kmin;
    }
    counterv[kmin]++;
    for (d = 0; d < dims; ++d) {
      clusterv2[0] += pointv[0];
    }
  }
#pragma endscop
}
