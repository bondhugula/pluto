#include <stdlib.h>
#include <cloog/polylib/cloog.h>

/**
 * Allocate complete state.
 */
CloogState *cloog_state_malloc(void)
{
	CloogState *state = cloog_core_state_malloc();
	state->backend = (CloogBackend *)malloc(sizeof(CloogBackend));
	if (!state->backend) 
		cloog_die("memory overflow.\n");
	state->backend->MAX_RAYS = 50;
	return state;
}

/**
 * Free state and backend independent parts.
 */
void cloog_state_free(CloogState *state)
{
	free(state->backend);
	cloog_core_state_free(state);
}
