#ifndef CLOOG_STATE_H
#define CLOOG_STATE_H

#if defined(CLOOG_POLYLIB)
#include <cloog/polylib/backend.h>
#elif defined(CLOOG_ISL)
#include <cloog/isl/backend.h>
#else
struct cloogbackend;
#endif
typedef struct cloogbackend CloogBackend;

#if defined(__cplusplus)
extern "C" {
#endif 

struct cloogstate {
  CloogBackend *backend;

  cloog_int_t one;
  cloog_int_t negone;

  int block_allocated;
  int block_freed;
  int block_max;

  int domain_allocated;
  int domain_freed;
  int domain_max;

  int loop_allocated;
  int loop_freed;
  int loop_max;

  int statement_allocated;
  int statement_freed;
  int statement_max;
};
typedef struct cloogstate CloogState;

CloogState *cloog_core_state_malloc(void);
CloogState *cloog_state_malloc(void);

void cloog_core_state_free(CloogState *state);
void cloog_state_free(CloogState *state);

#if defined(__cplusplus)
}
#endif 

#endif
