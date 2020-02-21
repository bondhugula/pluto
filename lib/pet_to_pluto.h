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
