#include <float.h>
#include <stdint.h>

void k_classify_local(float *pointv, uint32_t pointc, uint32_t point_lda,
                      float *clusterv, float *clusterv2, uint32_t *counterv,
                      uint32_t clusterc, uint32_t dims) {
  uint32_t p, k, d, kmin, dist_min, tmp, dist;

#define clusterv(d, k)                                                         \
  { clusterv[d * clusterc + k] }
#define pointv(d, p)                                                           \
  { pointv[d * point_lda + p] }
#define clusterv2(d, kmin)                                                     \
  { clusterv2[d * clusterc + kmin] }
#define counterv(kmin)                                                         \
  { couterv[kmin] }

#pragma scop
  for (p = 0; p < pointc; ++p) {
    kmin = -1;
    dist_min = FLT_MAX;
    for (k = 0; k < clusterc; ++k) {
      dist[k] = 0.0;
      for (d = 0; d < dims; ++d) {
        tmp = clusterv[d][k] - pointv(d, p);
        dist[k] += tmp * tmp;
      }
      dist_min[k] = (dist[k] < dist_min) ? dist[k] : dist_min;
      kmin[k] = (dist[k] < dist_min) ? k : kmin;
    }
    counterv[0] = counterv[0] + 1;
    for (d = 0; d < dims; ++d) {
      clusterv2[d][kmin] += pointv[d][p];
    }
  }
#pragma endscop
}
