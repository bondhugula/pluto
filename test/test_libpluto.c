#include "pluto/libpluto.h"
#include "isl/union_map.h"
#include "isl/union_set.h"

PlutoOptions *options;

void run_test_and_cleanup(const char *domains_str, const char *deps_str,
                          PlutoOptions *options, __isl_take isl_ctx *ctx) {

  isl_union_set *domains = isl_union_set_read_from_str(ctx, domains_str);
  isl_union_map *deps = isl_union_map_read_from_str(ctx, deps_str);

  isl_union_map *schedule = pluto_schedule(domains, deps, options);
  if (schedule) {
    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    printer = isl_printer_print_union_map(printer, schedule);
    printf("\n");
    isl_printer_free(printer);

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedule);
  } else {
    printf("No schedule found!\n");
  }

  // Test remapping.
  // TODO: test it properly instead of just checking it doesn't crash.
  Remapping remapping;
  pluto_get_remapping_str(domains_str, deps_str, options, &remapping);

  isl_union_set_free(domains);
  isl_union_map_free(deps);

  isl_ctx_free(ctx);
}

// Each test case should name statements S_0, S_1, ...

// CHECK-LABEL: *** TEST CASE 1 ***
// CHECK: T(S1): (i0, 0, 0)
// CHECK: T(S2): (i0, 1, i1)
// CHECK: T(S3): (i0, 2, 0)
void test1() {
  printf("\n\n*** TEST CASE 1 ***\n\n");

  isl_ctx *ctx = isl_ctx_alloc();
  const char *domains_str =
      "[n] -> {S_1[i0, i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 <= "
      "99; S_0[i0] : i0 >= 0 and i0 <= 99; S_2[i0] : i0 >= 0 and i0 <= 99 "
      "}";
  const char *deps_str =
      "[n] -> {S_1[i0, 99] -> S_0[1 + i0] : i0 >= 0 and i0 <= 98; S_1[i0, "
      "i1] -> S_1[i0, 1 + i1] : i0 >= 0 and i0 <= 99 and i1 >= 0 and i1 "
      "<= 98; S_1[i0, 99] -> S_1[1 + i0, 0] : i0 >= 0 and i0 <= 98; "
      "S_0[i0] -> S_1[i0, 0] : i0 >= 0 and i0 <= 99; S_2[i0] -> S_1[1 + "
      "i0, 0] : i0 >= 0 and i0 <= 98; S_0[i0] -> S_2[i0] : i0 >= 0 and i0 "
      "<= 99; S_1[i0, 99] -> S_2[i0] : i0 >= 0 and i0 <= 99 }";

  run_test_and_cleanup(domains_str, deps_str, options, ctx);
}

// CHECK-LABEL: *** TEST CASE 2 ***
// CHECK: T(S1): (i0, 1024i0)
// CHECK: T(S2): (i0, 1024i0+i1)
int test2() {
  printf("\n\n*** TEST CASE 2 ***\n\n");
  isl_ctx *ctx = isl_ctx_alloc();
  const char *domains_str =
      "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_1[i0, i1] : i0 >= 0 and "
      "i0 <= p_0 and i1 >= 0 and i1 <= p_3 and p_2 >= 0; S_0[i0] : i0 >= "
      "0 and i0 <= p_0}";
  const char *deps_str =
      "[p_0, p_1, p_2, p_3, p_4, p_5, p_7] -> { S_0[i0] -> S_1[o0, o1] : "
      "(exists (e0 = [(p_7)/8]: 8o1 = -p_5 + p_7 + 8192i0 - 8192o0 and "
      "8e0 = p_7 and i0 >= 0 and o0 <= p_0 and 8192o0 >= -8p_3 - p_5 + "
      "p_7 + 8192i0 and 8192o0 <= -p_5 + p_7 + 8192i0 and p_2 >= 0 and o0 "
      ">= 1 + i0)); S_1[i0, i1] -> S_0[o0] : (exists (e0 = [(p_1)/8], e1 "
      "= [(p_4)/8], e2 = [(-p_1 + p_7)/8184]: 8192o0 = p_5 - p_7 + 8192i0 "
      "+ 8i1 and 8e0 = p_1 and 8e1 = p_4 and 8184e2 = -p_1 + p_7 and i1 "
      ">= 0 and 8i1 <= 8192p_0 - p_5 + p_7 - 8192i0 and 8184i1 >= 1024 + "
      "1024p_1 - 1023p_5 - p_7 - 8380416i0 and p_2 >= 0 and p_7 <= -1 + "
      "p_5 and 8i1 >= 1 + 8p_3 + p_4 - p_5 - 8192i0 and i1 <= p_3 and i0 "
      ">= 0 and 8i1 >= 8192 - p_5 + p_7))}";

  run_test_and_cleanup(domains_str, deps_str, options, ctx);
}

// CHECK-LABEL: *** TEST CASE 3 ***
// CHECK: T(S1): (i0-i1, i0+i1)
void test_diamond_tiling() {
  printf("\n\n*** TEST CASE 3 ***\n\n");

  isl_ctx *ctx = isl_ctx_alloc();

  const char *domains_str = " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; }";
  // Dependences (1,-1) and (1,1).
  const char *deps_str =
      "[R, T] -> {"
      "S_0[i0, i1] -> S_0[i0 + 1, i1 - 1] : 0 <= i0 <= T - 1 and 1 <= i1 "
      "<= R - 2; "
      "S_0[i0, i1] -> S_0[i0 + 1, i1 + 1] : 0 <= i0 <= T - 1 and 1 <= i1 "
      "<= R - 2; }";

  run_test_and_cleanup(domains_str, deps_str, options, ctx);
}

// CHECK-LABEL: *** TEST CASE 4 ***
// CHECK: T(S1): (i0+i1, i0)
void test4() {
  printf("\n\n*** TEST CASE 4 ***\n\n");

  isl_ctx *ctx = isl_ctx_alloc();

  const char *domains_str =
      " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 < i1 <= R - 1; }";
  // Dependence (1, -1)
  const char *deps_str = "[R, T] -> {"
                         "S_0[i0, i1] -> S_0[i0 + 1, i1 - 1] : 1 "
                         "<= i0 <= T and 1 <= i1 <= R - 2; }";
  run_test_and_cleanup(domains_str, deps_str, options, ctx);
}

// CHECK-LABEL: *** TEST CASE 5 ***
// CHECK: T(S1): (i0, i1, 0)
// CHECK: T(S2): (i0+1, i1+1, 1)
void test5() {
  printf("\n\n*** TEST CASE 5 ***\n\n");
  isl_ctx *ctx = isl_ctx_alloc();
  const char *domains_str =
      " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; "
      "S_1[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; }";
  // Inter-statement deps; benefits from shifting.
  const char *deps_str =
      "[R, T] -> {"
      "S_0[i0, i1] -> S_1[i0 - 1, i1 - 1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; "
      "}"
      "S_0[i0, i1] -> S_1[i0 - 1, i1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; }"
      "S_0[i0, i1] -> S_1[i0 - 1, i1 + 1] : 1 <= i0 <= T and 1 <= i1 <= R - 2; "
      "}";

  run_test_and_cleanup(domains_str, deps_str, options, ctx);
}

// CHECK-LABEL: *** TEST CASE 6
// CHECK: T(S1): (i0, i1+i2, i1)
void test6_diamond_tiling_with_scalar_dimension() {
  printf("\n\n*** TEST CASE 6\n\n");
  isl_ctx *ctx = isl_ctx_alloc();
  const char *domains_str =
      " [R, T] -> { S_0[2, i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1;}";
  // Dependences (1,-1), (1,0), and (1,1)
  const char *deps_str =
      "[R, T] -> {"
      "S_0[2, i0, i1] -> S_0[2, i0 + 1, i1 - 1] : 1 <= i0 <= T and 1 <= "
      "i1 <= R - 2; }"
      "S_0[2, i0, i1] -> S_0[2, i0 + 1, i1] : 1 <= i0 <= T and 1 <= i1 <= "
      "R - 2; }"
      "S_0[2, i0, i1] -> S_0[2, i0 + 1, i1 + 1] : 1 <= i0 <= T and 1 <= "
      "i1 <= R - 2; }";

  run_test_and_cleanup(domains_str, deps_str, options, ctx);
}

int main() {
  options = pluto_options_alloc();
  options->tile = 1;
  options->parallel = 1;
  options->debug = 0;
  options->moredebug = 0;
  options->islsolve = 1;
  options->diamondtile = 1;
  options->fulldiamondtile = 0;

  test1();
  test2();
  test_diamond_tiling();
  test4();
  test5();
  test6_diamond_tiling_with_scalar_dimension();

  pluto_options_free(options);
}
