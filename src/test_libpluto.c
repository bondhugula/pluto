#include "pluto/libpluto.h"
#include "isl/union_set.h"
#include "isl/union_map.h"

int another() 
{
	isl_ctx *ctx = isl_ctx_alloc();
    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_1[i0, i1] : i0 >= 0 and i0 <= p_0 and i1 >= 0 and i1 <= p_3 and p_2 >= 0; S_0[i0] : i0 >= 0 and i0 <= p_0}");
    isl_union_map *deps = isl_union_map_read_from_str(ctx, "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_0[i0] -> S_1[o0, o1] : (exists (e0 = [(p_7)/8]: 8o1 = -p_5 + p_7 + 8192i0 - 8192o0 and 8e0 = p_7 and i0 >= 0 and o0 <= p_0 and 8192o0 >= -8p_3 - p_5 + p_7 + 8192i0 and 8192o0 <= -p_5 + p_7 + 8192i0 and p_2 >= 0 and o0 >= 1 + i0)); S_1[i0, i1] -> S_0[o0] : (exists (e0 = [(p_1)/8], e1 = [(p_4)/8], e2 = [(-p_1 + p_7)/8184]: 8192o0 = p_5 - p_7 + 8192i0 + 8i1 and 8e0 = p_1 and 8e1 = p_4 and 8184e2 = -p_1 + p_7 and i1 >= 0 and 8i1 <= 8192p_0 - p_5 + p_7 - 8192i0 and 8184i1 >= 1024 + 1024p_1 - 1023p_5 - p_7 - 8380416i0 and p_2 >= 0 and p_7 <= -1 + p_5 and 8i1 >= 1 + 8p_3 + p_4 - p_5 - 8192i0 and i1 <= p_3 and i0 >= 0 and 8i1 >= 8192 - p_5 + p_7))}");

	PlutoOptions *options = pluto_options_alloc();
    options->tile = 1;
    options->parallel = 1;
    options->debug = 1;

	isl_union_map *schedule = pluto_schedule(domains, deps, options);

    if (schedule) {
        isl_printer *printer = isl_printer_to_file(ctx, stdout);
        isl_printer_print_union_map(printer, schedule);
        printf("\n");
        isl_printer_free(printer);

        // Check if the schedule can be applied to the domain.
        domains = isl_union_set_apply(domains, schedule);
    }else{
        printf("No schedule\n");
    }


    isl_union_set_free(domains);
    isl_union_map_free(deps);

    pluto_options_free(options);

    isl_ctx_free(ctx);
}

int main() {
	isl_ctx *ctx = isl_ctx_alloc();
    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            "[n] -> {S_1[i0, i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= 99; S_0[i0] : i0 >= 0 and i0 <= 99; S_2[i0] : i0 >= 0 and i0 <= 99 }");
    isl_union_map *deps = isl_union_map_read_from_str(ctx, "[n] -> {S_1[i0, 99] -> S_0[1 + i0] : i0 >= 0 and i0 <= 98; S_1[i0, i1] -> S_1[i0, 1 + i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= 98; S_1[i0, 99] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_1[i0, 0] : i0 >= 0 and i0 <= 99; S_2[i0] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_2[i0] : i0 >= 0 and i0 <= 99; S_1[i0, 99] -> S_2[i0] : i0 >= 0 and i0 <= 99 }");

	PlutoOptions *options = pluto_options_alloc();
    options->tile = 1;
    options->parallel = 1;

	isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");
    isl_printer_free(printer);

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedule);

    isl_union_set_free(domains);
    isl_union_map_free(deps);

    pluto_options_free(options);

    isl_ctx_free(ctx);

    another();
}
