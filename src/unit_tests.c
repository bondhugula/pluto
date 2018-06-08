/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the file
 * `LICENSE' in the top-level directory of this distribution.
 *
 * unit_tests.c
 *
 * This file is meant to test certain functions in Pluto in standalone
 * mode
 *
 */
#include <stdio.h>
#include <assert.h>

#include "constraints.h"
#include "program.h"


void test_rank() {
    while (!feof(stdin)) {
        PlutoMatrix *mat = pluto_matrix_input(stdin);
        printf("Rank of matrix is %d\n", pluto_matrix_get_rank(mat));
        pluto_matrix_free(mat);
    }
}

