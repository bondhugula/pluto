#include "demosaic_func.h"

int N = 128;
#pragma parameter N 128 1000
int M = 128;
#pragma parameter M 128 1000

int main(void) {
  int x, Z;
  static int A[400];
  static int B[400];

#pragma scop
  for (x = 0; x < 2; x++) {
    for (z = 0; z < 2; z++) {
      A[2 * x + z] = 1;
    }
  }

  for (x = 0; x < 2; x++) {
    for (z = 0; z < 2; z++) {
      B[2 * x + z] = A[2 * x + z];
    }
  }
#pragma endscop

  writeImage(c);

  return (0);
}
