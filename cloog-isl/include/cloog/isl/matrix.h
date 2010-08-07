#ifndef CLOOG_ISL_MATRIX_H
#define CLOOG_ISL_MATRIX_H

#include <cloog/isl/backend.h>

#if defined(__cplusplus)
extern "C" 
  {
#endif 

typedef struct isl_basic_set CloogConstraintSet;

struct cloogequalities {
	int			  n;
	unsigned		  total_dim;
	CloogConstraintSet	**constraints;
	int			 *types;
};
typedef struct cloogequalities CloogEqualities;

typedef struct isl_constraint *CloogConstraint;

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
