#ifndef CLOOG_ISL_BACKEND_H
#define CLOOG_ISL_BACKEND_H

#include <isl_constraint.h>

struct cloogbackend {
	struct isl_ctx	*ctx;
	unsigned	ctx_allocated : 1;
};

#endif /* define _H */
