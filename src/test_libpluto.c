#include "pluto/libpluto.h"
#include "isl/union_set.h"
#include "isl/union_map.h"

PlutoOptions *options;

/* 
 * Each test case should name statements S_0, S_1, ...
 */
int test2() 
{
    printf("\n\n*** TEST CASE 2 ***\n\n");
    isl_ctx *ctx = isl_ctx_alloc();
    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_1[i0, i1] : i0 >= 0 and i0 <= p_0 and i1 >= 0 and i1 <= p_3 and p_2 >= 0; S_0[i0] : i0 >= 0 and i0 <= p_0}");
    isl_union_map *deps = isl_union_map_read_from_str(ctx, "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_0[i0] -> S_1[o0, o1] : (exists (e0 = [(p_7)/8]: 8o1 = -p_5 + p_7 + 8192i0 - 8192o0 and 8e0 = p_7 and i0 >= 0 and o0 <= p_0 and 8192o0 >= -8p_3 - p_5 + p_7 + 8192i0 and 8192o0 <= -p_5 + p_7 + 8192i0 and p_2 >= 0 and o0 >= 1 + i0)); S_1[i0, i1] -> S_0[o0] : (exists (e0 = [(p_1)/8], e1 = [(p_4)/8], e2 = [(-p_1 + p_7)/8184]: 8192o0 = p_5 - p_7 + 8192i0 + 8i1 and 8e0 = p_1 and 8e1 = p_4 and 8184e2 = -p_1 + p_7 and i1 >= 0 and 8i1 <= 8192p_0 - p_5 + p_7 - 8192i0 and 8184i1 >= 1024 + 1024p_1 - 1023p_5 - p_7 - 8380416i0 and p_2 >= 0 and p_7 <= -1 + p_5 and 8i1 >= 1 + 8p_3 + p_4 - p_5 - 8192i0 and i1 <= p_3 and i0 >= 0 and 8i1 >= 8192 - p_5 + p_7))}");

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

    isl_ctx_free(ctx);
}

void test1()
{
    isl_ctx *ctx = isl_ctx_alloc();
    isl_union_set *domains;
    isl_union_map *deps; 

    printf("\n\n*** TEST CASE 1 ***\n\n");
  
    domains = isl_union_set_read_from_str(ctx,
            "[n] -> {S_1[i0, i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= 99; S_0[i0] : i0 >= 0 and i0 <= 99; S_2[i0] : i0 >= 0 and i0 <= 99 }");
    deps = isl_union_map_read_from_str(ctx, "[n] -> {S_1[i0, 99] -> S_0[1 + i0] : i0 >= 0 and i0 <= 98; S_1[i0, i1] -> S_1[i0, 1 + i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= 98; S_1[i0, 99] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_1[i0, 0] : i0 >= 0 and i0 <= 99; S_2[i0] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_2[i0] : i0 >= 0 and i0 <= 99; S_1[i0, 99] -> S_2[i0] : i0 >= 0 and i0 <= 99 }");

    isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");
    isl_printer_free(printer);

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedule);

    isl_union_set_free(domains);
    isl_union_map_free(deps);


    isl_ctx_free(ctx);
}

/*
crasges when dependece vector is in the 3rd quadrant (-1, -1)
*/
void crash_negative_dep_vector() 
{
    printf("\n\n*** TEST CASE 4 ***\n\n");

    isl_ctx *ctx = isl_ctx_alloc();

    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 < i1 <= R - 1; }");
    isl_union_map *deps = isl_union_map_read_from_str(ctx,
        "[R, T] -> {"
        /* crashes with this dependence */
        "S_0[i0, i1] -> S_0[i0 - 1, i1 - 1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; }");
        /* does not crash with this dependence */
        //"S_0[i0, i1] -> S_0[i0 - 1, i1 + 1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; }");

    isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedule);

    isl_union_set_free(domains);
    isl_union_map_free(deps);


    isl_ctx_free(ctx);

}

/*
Works with options->lbtile = 1. Does not work with options->libtile = 0 &&
    options->partlbtile = 1.
i0: time
i1: row
*/
void test_diamond_tiling() 
{
    printf("\n\n*** TEST CASE 3 ***\n\n");

    isl_ctx *ctx = isl_ctx_alloc();

    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; }");
    isl_union_map *deps = isl_union_map_read_from_str(ctx,
        "[R, T] -> {"
        "S_0[i0, i1] -> S_0[i0 + 1, i1 - 1] : 0 <= i0 <= T - 1 and 1 <= i1 <= R - 2; "
        "S_0[i0, i1] -> S_0[i0 + 1, i1 + 1] : 0 <= i0 <= T - 1 and 1 <= i1 <= R - 2; }");
    isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");
    isl_printer_free(printer);

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedule);

    isl_union_set_free(domains);
    isl_union_map_free(deps);


    isl_ctx_free(ctx);
}

void test4() {
    isl_ctx *ctx = isl_ctx_alloc();
    isl_union_set *domains = isl_union_set_read_from_str(ctx,
            " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; S_1[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; }");
    isl_union_map *deps = isl_union_map_read_from_str(ctx,
        "[R, T] -> {"
        "S_0[i0, i1] -> S_1[i0 - 1, i1 - 1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; }"
        "S_0[i0, i1] -> S_1[i0 - 1, i1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; }"
        "S_0[i0, i1] -> S_1[i0 - 1, i1 + 1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; }");
    isl_union_map *schedule = pluto_schedule(domains, deps, options);

    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    isl_printer_print_union_map(printer, schedule);
    printf("\n");
    isl_printer_free(printer);

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedule);

    isl_union_set_free(domains);
    isl_union_map_free(deps);


    isl_ctx_free(ctx);
}

int main() 
{
    options = pluto_options_alloc();
    options->tile = 1;
    options->parallel = 1;
    options->debug = 0;
    options->moredebug = 0;
    options->islsolve = 1;
    options->partlbtile = 1;
    options->lbtile = 1;

    test1();
    test2();
    test_diamond_tiling();
    crash_negative_dep_vector();

    pluto_options_free(options);
}
