/******************************************************************************
 *               libpluto -  A library version of Pluto                       *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2012 Uday Bondhugula                                         *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License version 2.1 as published by the Free Software Foundation.          *
 *                                                                            *
 * This library is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
 * Lesser General Public License for more details.                            *
 *                                                                            *
 * A copy of the GNU Lesser General Public Licence can be found in the file
 * `LICENSE.LGPL2' in the top-level directory of this distribution.
 *
 * This file is part of libpluto.
 *
 */
#include "pluto/pluto.h"
#include "isl/union_map.h"
#include "isl/union_set.h"

void run_test_and_cleanup(const char *domains_str, const char *deps_str,
                          PlutoOptions *options, __isl_take isl_ctx *ctx) {

  isl_union_set *domains = isl_union_set_read_from_str(ctx, domains_str);
  isl_union_map *deps = isl_union_map_read_from_str(ctx, deps_str);

  isl_union_map *schedules =
      pluto_transform(isl_union_set_copy(domains), deps, NULL, NULL, options);
  if (schedules) {
    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    printer = isl_printer_print_union_map(printer, schedules);
    printf("\n");
    isl_printer_free(printer);

    // Check if the schedule can be applied to the domain.
    domains = isl_union_set_apply(domains, schedules);
  } else {
    printf("No schedule found!\n");
    isl_union_map_free(schedules);
  }
  isl_union_set_free(domains);

  // Test remapping.
  // TODO: test it properly instead of just checking it doesn't crash.
  Remapping remapping;
  pluto_get_remapping_str(domains_str, deps_str, options, &remapping);
  pluto_remapping_free(remapping);

  isl_ctx_free(ctx);
}

// Each test case should name statements S_0, S_1, ...

// CHECK-LABEL: *** TEST CASE 1 ***
// CHECK: T(S1): (i0, 0, 0)
// CHECK: T(S2): (i0, 1, i1)
// CHECK: T(S3): (i0, 2, 0)
void test1(PlutoOptions *options) {
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
void test2(PlutoOptions *options) {
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
void test_diamond_tiling(PlutoOptions *options) {
  printf("\n\n*** TEST CASE 3 ***\n\n");

  isl_ctx *ctx = isl_ctx_alloc();

  const char *domains_str =
      " [R, T] -> { S_0[i0, i1] : 0 <= i0 <= T and 0 <= i1 <= R - 1; }";
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
void test4(PlutoOptions *options) {
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
void test5(PlutoOptions *options) {
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
void test6_diamond_tiling_with_scalar_dimension(PlutoOptions *options) {
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

// CHECK-LABEL: *** TEST CASE test_lib_pluto_schedule ***
// CHECK: T(S1): (i0, i1, 0)
// CHECK: T(S2): (i0, i2, i1)
void test_lib_pluto_schedule(PlutoOptions *options) {
  printf("\n\n*** TEST CASE test_lib_pluto_schedule ***\n\n");

  isl_ctx *ctx = isl_ctx_alloc();

  isl_union_map *schedules = isl_union_map_read_from_str(
      ctx,
      "[p1, p2, p0] -> { S_0[i0, i1] -> [o0, o1, o2, o3] : o0 = 0 and o1 = "
      "0 and o2 = i0 and o3 = 0 and i0 >= 0 and i0 < p0 and i1 >= 0 and "
      "i1 < p1; S_1[i0, i1, i2] -> [o0, o1, o2, o3] : o0 = 1 and o1 = i0 "
      "and o2 = i1 and o3 = i2 and i0 >= 0 and i0 < p1 and i1 >= 0 and i1 "
      "< p2 and i2 >= 0 and i2 <= 100 }");

  isl_union_map *reads = isl_union_map_read_from_str(
      ctx,
      "[p1, p2, p0] -> { S_0[i0, i1]->M0[o0, o1] : o0 = i0 and o1 = i1;    "
      " S_1[i0, i1, i2]->M0[o0, o1] : o0 = i0 and o1 = i2 }");

  isl_union_map *writes = isl_union_map_read_from_str(
      ctx,
      "[p1, p2, p0] -> { S_1[i0, i1, i2]->M1[o0, o1] : o0 = i2 and o1 = i0}");

  isl_union_map *result = pluto_schedule(schedules, reads, writes, options);

  if (result) {
    isl_printer *printer = isl_printer_to_file(ctx, stdout);
    printer = isl_printer_print_union_map(printer, result);
    printf("\n");
    isl_printer_free(printer);
  }

  isl_union_map_free(result);
  isl_ctx_free(ctx);
}

int main() {
  PlutoOptions *options = pluto_options_alloc();
  options->tile = 1;
  options->parallel = 1;
  options->debug = 0;
  options->moredebug = 0;
  options->islsolve = 1;
  options->diamondtile = 1;
  options->fulldiamondtile = 0;
  options->isldepaccesswise = 0;

  test1(options);
  test2(options);
  test_diamond_tiling(options);
  test4(options);
  test5(options);
  test6_diamond_tiling_with_scalar_dimension(options);

  options->tile = 0;
  options->diamondtile = 0;
  options->fulldiamondtile = 0;
  test_lib_pluto_schedule(options);

  pluto_options_free(options);
}
