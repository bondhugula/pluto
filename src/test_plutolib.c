#include "pluto/libpluto.h"
#include "isl/union_set.h"
#include "isl/union_map.h"

int main() {
	isl_ctx *ctx = isl_ctx_alloc();
	isl_union_set *domains = isl_union_set_read_from_str(ctx,
		"{ Stmt_for_body41[i0, i1]: 0 <= i0, i1 <= 10}");
	isl_union_map *deps = isl_union_map_read_from_str(ctx, "{ }");

	PlutoOptions *options = pluto_options_alloc();
    options->tile = 1;

	isl_union_map *schedule = pluto_schedule(domains, deps, options);

	isl_union_map_dump(schedule);

    isl_union_set_free(domains);
    isl_union_map_free(deps);
    isl_union_map_free(schedule);

    pluto_options_free(options);

    isl_ctx_free(ctx);
}
