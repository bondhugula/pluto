#ifndef _PET_TO_PLUTO_H_
#define _PET_TO_PLUTO_H_

#include "pet.h"
#include "pluto.h"

#if defined(__cplusplus)
extern "C" {
#endif

PlutoProg *pet_to_pluto_prog(struct pet_scop *pscop, isl_ctx *, PlutoOptions *);

#if defined(__cplusplus)
}
#endif

#endif
