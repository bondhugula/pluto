#include <stdint.h>

typedef struct plutoContext PlutoContext;

#ifndef _PLUTO_MATRIX_H
#define _PLUTO_MATRIX_H
/* A matrix */
struct pluto_matrix {
  /* The values */
  int64_t **val;

  unsigned nrows;
  unsigned ncols;

  /* Pre-allocated number of rows */
  int alloc_nrows;
  int alloc_ncols;

  PlutoContext *context;
};
typedef struct pluto_matrix PlutoMatrix;

#endif
