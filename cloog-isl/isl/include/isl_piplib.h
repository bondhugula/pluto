#ifndef ISL_PIPLIB_H
#define ISL_PIPLIB_H

#include <isl_ctx.h>
#include <isl_int.h>
#include <isl_map.h>
#ifndef ISL_PIPLIB
#error "no piplib"
#endif

#include <piplib/piplibMP.h>

void isl_seq_cpy_to_pip(Entier *dst, isl_int *src, unsigned len);
void isl_seq_cpy_from_pip(isl_int *dst, Entier *src, unsigned len);

PipMatrix *isl_basic_map_to_pip(struct isl_basic_map *bmap, unsigned pip_param,
			 unsigned extra_front, unsigned extra_back);

#endif
