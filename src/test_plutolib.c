#include "pluto/libpluto.h"
#include "isl/union_set.h"
#include "isl/union_map.h"

int main() {
	isl_ctx *ctx = isl_ctx_alloc();
	isl_union_set *domains = isl_union_set_read_from_str(ctx,
		"{ Stmt_for_body8[i0, i1, i2] : i0 >= 0 and i0 <= 1023 and i1 >= 0 and i1 <= 1023 and i2 >= 0 and i2 <= 1023; Stmt_for_body3[i0, i1] : i0 >= 0 and i0 <= 1023 and i1 >= 0 and i1 <= 1023 }");
	isl_union_map *deps = isl_union_map_read_from_str(ctx, "{ S_1_w_0[i0, i1] -> S_0_r_0[i0, i1, 0] : i0 >= 0 and i0 <= 1023 and i1 >= 0 and i1 <= 1023; S_0_w_0[i0, i1, i2] -> S_0_r_0[i0, i1, 1 + i2] : i0 >= 0 and i0 <= 1023 and i1 >= 0 and i1 <= 1023 and i2 >= 0 and i2 <= 1022 }");

	PlutoOptions *options = pluto_options_alloc();
    options->tile = 1;
    options->parallel = 1;

	isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");

    isl_union_set_free(domains);
    isl_union_map_free(deps);
    isl_union_map_free(schedule);

    pluto_options_free(options);

    isl_ctx_free(ctx);
}
