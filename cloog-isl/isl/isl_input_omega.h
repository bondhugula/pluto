#ifndef ISL_INPUT_OMEGA_H
#define ISL_INPUT_OMEGA_H

#include <stdio.h>
#include <isl_set.h>

struct isl_basic_set *isl_basic_set_read_from_file_omega(
		struct isl_ctx *ctx, FILE *input);
struct isl_basic_set *isl_basic_set_read_from_str_omega(
		struct isl_ctx *ctx, const char *str);
struct isl_basic_map *isl_basic_map_read_from_file_omega(
		struct isl_ctx *ctx, FILE *input);
struct isl_basic_map *isl_basic_map_read_from_str_omega(
		struct isl_ctx *ctx, const char *str);

#endif
