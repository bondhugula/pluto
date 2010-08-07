#ifndef CLOOG_ISL_H
#define CLOOG_ISL_H

#ifndef CLOOG_ISL
#define CLOOG_ISL
#endif
#include <cloog/cloog.h>
#include <cloog/isl/domain.h>

CloogState *cloog_isl_state_malloc(struct isl_ctx *ctx);

#endif /* define _H */
