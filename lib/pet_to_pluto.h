/*
 * Pluto: An automatic parallelier and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This software is available under the MIT license. Please see LICENSE
 * in the top-level directory for details.
 *
 * This file is part of libpluto.
 *
 */
#ifndef _PET_TO_PLUTO_H_
#define _PET_TO_PLUTO_H_

#include "pet.h"

typedef struct plutoProg PlutoProg;
typedef struct plutoContext PlutoContext;

#if defined(__cplusplus)
extern "C" {
#endif

PlutoProg *pet_to_pluto_prog(struct pet_scop *pscop, isl_ctx *,
                             PlutoContext *context);

#if defined(__cplusplus)
}
#endif

#endif
