/******************************************************************************
 *               libpluto -  A library version of Pluto
 ******************************************************************************
 *
 * Copyright (C) 2012 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE in the
 * top-level directory for details.
 *
 * This file is part of libpluto.
 *
 */
#ifndef _PLUTO_MATRIX_H
#define _PLUTO_MATRIX_H

#include <stdint.h>

typedef struct plutoContext PlutoContext;

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
