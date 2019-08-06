#ifndef _PLUTO_CODEGEN_IF_H_

#define _PLUTO_CODEGEN_IF_H_

#include "program.h"
#include "osl/extensions/loop.h"

osl_loop_p pluto_get_vector_loop_list(const PlutoProg *prog);
osl_loop_p pluto_get_parallel_loop_list(const PlutoProg *prog,
                                        int vloopsfound);
#endif
