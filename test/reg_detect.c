#include <math.h>
#include <stdio.h>

#include "decls.h"
#include "util.h"

double t_start, t_end;

// simplified regularity detection algorithm
int main(int argc, char *argv[]) {
  int i, j, cnt;
  int sum_tang[maxgrid][maxgrid];
  int mean[maxgrid][maxgrid];
  int diff[maxgrid][maxgrid][length];
  int sum_diff[maxgrid][maxgrid][length];
  int tangent[maxrgc];        // input
  int path[maxgrid][maxgrid]; // output

  for (j = 0; j <= maxgrid - 1; j++) {
    sum_tang[j][j] = tangent[(maxgrid + 1) * j];
    for (i = j + 1; i <= maxgrid - 1; i++) {
      sum_tang[j][i] = sum_tang[j][i - 1] + tangent[i + maxgrid * j];
    }
  }
  IF_TIME(t_start = rtclock());

  long s = 0;
  int t;
#define NITER 1000
  for (t = 0; t < NITER; t++) {

    /* pluto start (maxgrid,length,maxrgc) */
    for (j = 0; j <= maxgrid - 1; j++) {
      for (i = j; i <= maxgrid - 1; i++) {
        for (cnt = 0; cnt <= length - 1; cnt++) {
          diff[j][i][cnt] = sum_tang[j][i];
        }
      }
    }

    for (j = 0; j <= maxgrid - 1; j++) {
      for (i = j; i <= maxgrid - 1; i++) {
        sum_diff[j][i][0] = diff[j][i][0];
        for (cnt = 1; cnt <= length - 1; cnt++) {
          sum_diff[j][i][cnt] = sum_diff[j][i][cnt - 1] + diff[j][i][cnt];
        }
        mean[j][i] = sum_diff[j][i][length - 1];
      }
    }

    for (i = 0; i <= maxgrid - 1; i++) {
      path[0][i] = mean[0][i];
    }

    for (j = 1; j <= maxgrid - 1; j++) {
      for (i = j; i <= maxgrid - 1; i++) {
        path[j][i] = path[j - 1][i - 1] + mean[j][i];
      }
    }
    /* pluto end */
    s += path[maxgrid - 1][1];
  }

  IF_TIME(t_end = rtclock());
  IF_TIME(printf("%0.6lfs\n", t_end - t_start));
}
