// CHECK: Output written
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <assert.h>

#include "decls.h"

#include "util.h"

int main() {
  int i, j, k, l, t;

  double t_start, t_end;

  init_array();

  IF_TIME(t_start = rtclock());

#pragma scop
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      x1[i] = x1[i] + a[i][j] * y_1[j];
    }
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      x2[i] = x2[i] + a[j][i] * y_2[j];
    }
  }
#pragma endscop

  IF_TIME(t_end = rtclock());
  IF_TIME(printf("%0.6lfs\n", t_end - t_start));

#ifdef TEST
  print_array();
#endif
  return 0;
}
