#include <stdint.h>

#ifndef _PLUTO_MATRIX_H
#define _PLUTO_MATRIX_H
/* A matrix */
struct plutoMatrix {
  /* The values */
  int64_t **val;

  unsigned nrows;
  unsigned ncols;

  /* Pre-allocated number of rows */
  int alloc_nrows;
  int alloc_ncols;
};
typedef struct plutoMatrix PlutoMatrix;

#endif
