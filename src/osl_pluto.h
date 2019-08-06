#ifndef _OSL_PLUTO_H_
#define _OSL_PLUTO_H_

#include "osl/scop.h"

#include "program.h"

PlutoProg *osl_scop_to_pluto_prog(osl_scop_p scop, PlutoOptions *options);

void pluto_populate_scop(osl_scop_p scop, PlutoProg *prog,
                         PlutoOptions *options);

#endif
