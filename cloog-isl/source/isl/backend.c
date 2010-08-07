#include <cloog/isl/cloog.h>

/**
 * Allocate and initialize full state.
 */
CloogState *cloog_state_malloc(void)
{
	return cloog_isl_state_malloc(NULL);
}

/**
 * Allocate and initialize full state for isl backend.
 */
CloogState *cloog_isl_state_malloc(struct isl_ctx *ctx)
{
	CloogState *state;
	int allocated = !ctx;

	state = cloog_core_state_malloc();
	if (!ctx)
		ctx = isl_ctx_alloc();
	state->backend = isl_alloc_type(ctx, CloogBackend);
	state->backend->ctx = ctx;
	state->backend->ctx_allocated = allocated;
	return state;
}

/**
 * Free state and backend independent parts.
 */
void cloog_state_free(CloogState *state)
{
	if (state->backend->ctx_allocated)
		isl_ctx_free(state->backend->ctx);
	free(state->backend);
	cloog_core_state_free(state);
}
