#include "pluto/libpluto.h"
#include "isl/union_set.h"
#include "isl/union_map.h"

int main() {
	isl_ctx *ctx = isl_ctx_alloc();
    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            "{S_1[i0, i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= 99; S_0[i0] : i0 >= 0 and i0 <= 99; S_2[i0] : i0 >= 0 and i0 <= 99 }");
    isl_union_map *deps = isl_union_map_read_from_str(ctx, "{S_1[i0, 99] -> S_0[1 + i0] : i0 >= 0 and i0 <= 98; S_1[i0, i1] -> S_1[i0, 1 + i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= 98; S_1[i0, 99] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_1[i0, 0] : i0 >= 0 and i0 <= 99; S_2[i0] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_2[i0] : i0 >= 0 and i0 <= 99; S_1[i0, 99] -> S_2[i0] : i0 >= 0 and i0 <= 99 }");

	PlutoOptions *options = pluto_options_alloc();
    options->tile = 1;
    options->parallel = 1;

	isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");
    isl_printer_free(printer);


    isl_union_set_free(domains);
    isl_union_map_free(deps);
    isl_union_map_free(schedule);

    pluto_options_free(options);

    isl_ctx_free(ctx);
}
