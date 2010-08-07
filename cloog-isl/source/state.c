#include <stdlib.h>
#include "../include/cloog/cloog.h"

/**
 * Allocate state and initialize backend independent part.
 */
CloogState *cloog_core_state_malloc(void)
{
  CloogState *state;

  state = (CloogState *)malloc(sizeof(CloogState));
  if (!state) 
    cloog_die("memory overflow.\n");

  state->backend = NULL;

  cloog_int_init(state->one);
  cloog_int_set_si(state->one, 1);
  cloog_int_init(state->negone);
  cloog_int_set_si(state->negone, -1);

  state->block_allocated = 0;
  state->block_freed = 0;
  state->block_max = 0;

  state->domain_allocated = 0;
  state->domain_freed = 0;
  state->domain_max = 0;

  state->loop_allocated = 0;
  state->loop_freed = 0;
  state->loop_max = 0;

  state->statement_allocated = 0;
  state->statement_freed = 0;
  state->statement_max = 0;

  return state;
}

/**
 * Free state.
 */
void cloog_core_state_free(CloogState *state)
{
  cloog_int_clear(state->one);
  cloog_int_clear(state->negone);
  free(state);
}
