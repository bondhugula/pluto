#ifndef _OSL_PLUTO_H_
#define _OSL_PLUTO_H_

#include "osl/scop.h"

typedef struct plutoContext PlutoContext;
typedef struct plutoProg PlutoProg;

#if defined(__cplusplus)
extern "C" {
#endif

PlutoProg *osl_scop_to_pluto_prog(osl_scop_p scop, PlutoContext *context);

void pluto_populate_scop(osl_scop_p scop, PlutoProg *prog,
                         PlutoContext *context);

#if defined(__cplusplus)
}
#endif

#endif // _OSL_PLUTO_H_
