/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software: you can redistribute it and/or modify
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
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "math_support.h"
#include "constraints.h"
#include "pluto.h"
#include "program.h"

#include <isl/constraint.h>
#include <isl/mat.h>
#include <isl/set.h>
#include "candl/candl.h"


#define CONSTRAINTS_SIMPLIFY_THRESHOLD 5000
#define MAX_FARKAS_CST  2000

/**
 *
 * Each constraint row has the following format
 *
 *      [dep distance bound | mapping coeff.s for S1, S2,... |constant]
 * Size:[       npar+1      | (nvar+1)*nstmts                | 1      ]
 *
 * npar - number of parameters in whole program
 * nvar - number of parameters in whole program
 */

/* Builds validity constraints for a dependence */
static void compute_permutability_constraints_dep(Dep *dep, PlutoProg *prog)
{
    PlutoConstraints *cst;
    PlutoConstraints *tiling_valid_cst, *bounding_func_cst;
    int nstmts, nvar, npar, src_stmt, dest_stmt, j, k, r;
    int src_offset, dest_offset;
    PlutoMatrix *phi;
    Stmt **stmts;

    nvar = prog->nvar;
    npar = prog->npar;
    stmts = prog->stmts;
    nstmts = prog->nstmts;

    /* IMPORTANT: It's assumed that all statements are of dimensionality nvar */

    IF_DEBUG(printf("[pluto] compute permutability constraints: Dep %d\n", dep->id+1););

    dest_stmt = dep->dest;
    src_stmt = dep->src;

    IF_MORE_DEBUG(printf("[pluto] compute permutability constraints: Dep %d\n", dep->id+1););

    /* Convert everything to >= 0 form */
    PlutoConstraints *dpoly = pluto_constraints_dup(dep->dpolytope);

    if (src_stmt != dest_stmt) {
        phi = pluto_matrix_alloc(2*nvar+npar+1, 2*(nvar+1)+1);
        pluto_matrix_set(phi, 0);

        for (r=0; r<nvar; r++) {
            /* Source stmt */
            phi->val[r][r] = -1;
        }
        for (r=nvar; r<2*nvar; r++) {
            /* Dest stmt */
            phi->val[r][(nvar+1)+(r-nvar)] = 1;
        }
        /* No parametric shifts: all zero for 2*nvar to 2*nvar+npar */

        /* Translation coefficients */
        phi->val[2*nvar+npar][(nvar+1)+nvar] = 1;
        phi->val[2*nvar+npar][nvar] = -1;
    }else{
        phi = pluto_matrix_alloc(2*nvar+npar+1, (nvar+1)+1);
        pluto_matrix_set(phi, 0);

        for (r=0; r<nvar; r++) {
            /* Source stmt */
            phi->val[r][r] = -1;
        }
        for (r=nvar; r<2*nvar; r++) {
            /* Dest stmt */
            phi->val[r][r-nvar] = 1;
        }
        /* No parametric shifts: so all zero for 2*nvar to 2*nvar+npar-1 */

        /* Translation coefficients cancel out;
         * so nothing for 2*nvar+npar */
    }

    /* Apply Farkas lemma for tiling validity constraints */
    tiling_valid_cst = farkas_lemma_affine(dpoly, phi);
    
    pluto_matrix_free(phi);

    if (src_stmt != dest_stmt) {
        phi = pluto_matrix_alloc(2*nvar+npar+1, npar+1+2*(nvar+1)+1);
        pluto_matrix_set(phi, 0);

        for (r=0; r<nvar; r++) {
            /* Source stmt */
            phi->val[r][npar+1+r] = 1;
        }
        for (r=nvar; r<2*nvar; r++) {
            /* Dest stmt */
            phi->val[r][npar+1+(nvar+1)+(r-nvar)] = -1;
        }
        for (r=2*nvar; r<2*nvar+npar; r++) {
            /* for \vec{u} - parametric bounding function */
            phi->val[r][r-2*nvar] = 1;
        }

        /* Translation coefficients of statements */
        phi->val[2*nvar+npar][npar+1+nvar] = 1;
        phi->val[2*nvar+npar][npar+1+(nvar+1)+nvar] = -1;
        /* for w */
        phi->val[2*nvar+npar][npar] = 1;
    }else{
        phi = pluto_matrix_alloc(2*nvar+npar+1, npar+1+(nvar+1)+1);
        pluto_matrix_set(phi, 0);

        for (r=0; r<nvar; r++) {
            /* Source stmt */
            phi->val[r][npar+1+r] = 1;
        }
        for (r=nvar; r<2*nvar; r++) {
            /* Dest stmt */
            phi->val[r][npar+1+(r-nvar)] = -1;
        }
        for (r=2*nvar; r<2*nvar+npar; r++) {
            /* for u */
            phi->val[r][r-2*nvar] = 1;
        }
        /* Statement's translation coefficients cancel out */

        /* for w */
        phi->val[2*nvar+npar][npar] = 1;
    }

    /* Apply Farkas lemma for bounding function constraints */
    bounding_func_cst = farkas_lemma_affine(dpoly, phi);

    pluto_matrix_free(phi);
    pluto_constraints_free(dpoly);

    /* Aggregate permutability and bounding function constraints together in
     * global format */

    /* Initialize everything to zero */
    cst = pluto_constraints_alloc(tiling_valid_cst->nrows + bounding_func_cst->nrows, CST_WIDTH);
    cst->ncols = CST_WIDTH;

    for (k=0; k<tiling_valid_cst->nrows+bounding_func_cst->nrows; k++)   {
        for (j=0; j<cst->ncols; j++)  {
            cst->val[k][j] = 0;
        }
    }

    src_offset = npar+1+src_stmt*(nvar+1);
    dest_offset = npar+1+dest_stmt*(nvar+1);

    /* Permutability constraints */
    if (!IS_RAR(dep->type)) {
        for (k=0; k<tiling_valid_cst->nrows; k++)   {
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows+k][src_offset+j] = tiling_valid_cst->val[k][j];
                if (src_stmt != dest_stmt) {
                    cst->val[cst->nrows+k][dest_offset+j] = tiling_valid_cst->val[k][nvar+1+j];
                }
            }
            /* constant part */
            if (src_stmt == dest_stmt)  {
                cst->val[cst->nrows+k][cst->ncols-1] = tiling_valid_cst->val[k][nvar+1];
            }else{
                cst->val[cst->nrows+k][cst->ncols-1] = tiling_valid_cst->val[k][2*nvar+2];
            }
        }
        cst->nrows = tiling_valid_cst->nrows;
    }

    if (!options->nodepbound)   {
        /* Add bounding constraints */
        src_offset = npar+1+src_stmt*(nvar+1);
        dest_offset = npar+1+dest_stmt*(nvar+1);

        for (k=0; k<bounding_func_cst->nrows; k++)   {
            for (j=0; j<npar+1; j++)  {
                cst->val[cst->nrows+k][j] = bounding_func_cst->val[k][j];
            }
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows+k][src_offset+j] = bounding_func_cst->val[k][npar+1+j];
                if (src_stmt != dest_stmt) {
                    cst->val[cst->nrows+k][dest_offset+j] = bounding_func_cst->val[k][npar+1+nvar+1+j];
                }
            }
            /* constant part */
            if (src_stmt == dest_stmt)  {
                cst->val[cst->nrows+k][cst->ncols-1] = bounding_func_cst->val[k][npar+1+nvar+1];
            }else{
                cst->val[cst->nrows+k][cst->ncols-1] = bounding_func_cst->val[k][npar+1+2*nvar+2];
            }
        }
        cst->nrows += bounding_func_cst->nrows;
    }


    /* Coefficients of those variables that don't appear in the outer loop
     * are useless */
    for (k=0; k<nvar; k++)    {
        if (!stmts[src_stmt]->is_orig_loop[k])  {
            for (j=0; j < cst->nrows; j++)   {
                cst->val[j][src_offset+k] = 0;
            }
        }
        if (src_stmt != dest_offset && !stmts[dest_stmt]->is_orig_loop[k])  {
            for (j=0; j < tiling_valid_cst->nrows+bounding_func_cst->nrows; j++)   {
                cst->val[j][dest_offset+k] = 0;
            }
        }
    }

    PlutoConstraints *bounding_cst = NULL;

    /* Copy only the bounding constraints */
    if (!options->nodepbound) {

		bounding_cst = pluto_constraints_alloc(bounding_func_cst->nrows, CST_WIDTH);
		bounding_cst->ncols = CST_WIDTH;
		bounding_cst->nrows = bounding_func_cst->nrows;

		for (k=0; k<bounding_func_cst->nrows; k++)   {
			for (j=0; j<(bounding_cst)->ncols; j++)  {
				(bounding_cst)->val[k][j] = 0;
			}
		}

		assert(cst->ncols == bounding_cst->ncols);

		for (k=0; k<bounding_func_cst->nrows; k++)   {
			for (j=0; j<bounding_cst->ncols; j++)  {
				bounding_cst->val[k][j] = cst->val[tiling_valid_cst->nrows+k][j];
			}
		}
    }

    pluto_constraints_free(tiling_valid_cst);
    pluto_constraints_free(bounding_func_cst);

    free(dep->valid_cst);
    dep->valid_cst = cst;

    free(dep->bounding_cst);
    dep->bounding_cst = bounding_cst;
}

/* This function itself is NOT thread-safe for the same PlutoProg */
PlutoConstraints *get_permutability_constraints(PlutoProg *prog)
{
    int i, inc, nstmts, nvar, npar, ndeps;
    PlutoConstraints *globcst;
    Dep **deps;

    nstmts = prog->nstmts;
    ndeps = prog->ndeps;
    deps = prog->deps;
    nvar = prog->nvar;
    npar = prog->npar;

    int total_cst_rows = 0;

    /* Compute the constraints and store them */
    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];

        if (options->rar == 0 && IS_RAR(dep->type)) {
            continue;
        }

        if (dep->valid_cst == NULL) {
            /* First time, get the constraints */
            compute_permutability_constraints_dep(dep, prog);

            IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", 
                        i+1, dep->valid_cst->nrows));
            total_cst_rows += dep->valid_cst->nrows;
            // IF_DEBUG(fprintf(stdout, "Constraints for dep: %d\n", i+1));
            // IF_DEBUG(pluto_constraints_pretty_print(stdout, dep->valid_cst));
        }
    }

    if (!prog->globcst) {
        prog->globcst = pluto_constraints_alloc(total_cst_rows, CST_WIDTH);
    }

    globcst = prog->globcst;

    globcst->ncols = CST_WIDTH;
    globcst->nrows = 0;

    /* Add constraints to globcst */
    for (i = 0, inc = 0; i < ndeps; i++) {
        Dep *dep = deps[i];

		/* print_polylib_visual_sets("BB_cst", dep->bounding_cst); */

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

        /* For debugging (skip deps listed here) */
        FILE *fp = fopen("skipdeps.txt", "r");
        if (fp) {
            int num;
            int found = 0;
            while (!feof(fp)) {
                fscanf(fp, "%d", &num);
                if (i == num-1) {
                    found = 1;
                    break;
                }
            }
            fclose(fp);
            if (found) {
                printf("Skipping dep %d\n", num);
                continue;
            }
        }

        /* Note that dependences would be marked satisfied (in
         * pluto_auto_transform) only after all possible independent solutions
         * are found to the formulation
         */
        if (dep_is_satisfied(dep) && dep->bounding_cst){
			/* Add only the bounding constraints when a dep is satisfied */
            pluto_constraints_add(globcst, dep->bounding_cst);
            //if(options->data_dist){
            //pluto_constraints_add(globcst, dep_bounding_cst[i]);
            //pluto_constraints_add(globcst, depcst[i]);
            //}
			continue;
        }

        /* Subsequent calls can just use the old ones */
        pluto_constraints_add(globcst, dep->valid_cst);
        /* print_polylib_visual_sets("global", dep->valid_cst); */

        IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", i+ 1,
                    globcst->nrows));
        /* This is for optimization as opposed to for correctness. We will
         * simplify constraints only if it crosses the threshold: at the time
         * it crosses the threshold or at 1000 increments thereon. Simplifying
         * each time slows down Pluto whenever there are few hundreds of
         * dependences. Not simplifying at all also leads to a slow down
         * because it leads to a large globcst and a number of constraits in
         * it are redundant */
        if (globcst->nrows >= CONSTRAINTS_SIMPLIFY_THRESHOLD + (1000*inc) &&
                globcst->nrows - dep->valid_cst->nrows < 
                CONSTRAINTS_SIMPLIFY_THRESHOLD + (1000*inc)) {
            pluto_constraints_simplify(globcst);
            IF_DEBUG(fprintf(stdout,
                        "After dep: %d; num_constraints_simplified: %d\n", i+1,
                        globcst->nrows));
            if (globcst->nrows >= CONSTRAINTS_SIMPLIFY_THRESHOLD + (1000*inc)) {
                inc++;
            }
        }
    }

    pluto_constraints_simplify(globcst);

    IF_DEBUG(fprintf(stdout, "After all dependences: num constraints: %d, num variables: %d\n",
                globcst->nrows, globcst->ncols - 1));
    IF_DEBUG2(pluto_constraints_pretty_print(stdout, globcst));

    return globcst;
}


/*
 * Construct a PlutoMatrix with the same content as the given isl_mat.
 */
/*
 * Returns linear independence constraints for a single statement.
 *
 * In particular, if H contains the first rows of an affine transformation,
 * then return a constraint on the coefficients of the next row that
 * ensures that this next row is linearly independent of the first rows.
 * Furthermore, the constraint is constructed in such a way that it allows
 * for a solution when combined with the other constraints on the coefficients
 * (currcst), provided any such constraint can be constructed.
 *
 * We do this by computing a basis for the null space of H and returning
 * a constraint that enforces the sum of these linear expressions
 * over the coefficients to be strictly greater than zero.
 * In this sum, some of the linear expressions may be negated to ensure
 * that a solution exists.
 *
 * The return value is a list of constraints, the first *orthonum corresponding
 * to the linear expressions that form a basis of the null space
 * and the final constraint the actual linear independence constraint.
 *
 * If the null space is 0-dimensional, *orthonum will be zero and the return
 * value is NULL
 */
PlutoConstraints **get_stmt_ortho_constraints(Stmt *stmt, const PlutoProg *prog,
        const PlutoConstraints *currcst, int *orthonum)
{
    int i, j, k, p, q;
    PlutoConstraints **orthcst;
    isl_ctx *ctx;
    isl_mat *h;
    isl_basic_set *isl_currcst;

    IF_DEBUG(printf("[pluto] get_stmt_ortho constraints S%d\n", stmt->id+1););

    int nvar = prog->nvar;
    int npar = prog->npar;
    int nstmts = prog->nstmts;
    HyperplaneProperties *hProps = prog->hProps;

    if (pluto_stmt_get_num_ind_hyps(stmt) >= stmt->dim_orig) {
        *orthonum = 0;
        return NULL;
    }

    /* Get rid of the variables that don't appear in the domain of this
     * statement and also beta rows */
    for (i = 0, p = 0; i < nvar; i++) {
        if (stmt->is_orig_loop[i]) {
            p++;
        }
    }

    assert(stmt->trans != NULL);

    for (j = 0, q = 0; j < stmt->trans->nrows; j++) {
        if (hProps[j].type != H_SCALAR) {
            q++;
        }
    }

    ctx = isl_ctx_alloc();
    assert(ctx);

    h = isl_mat_alloc(ctx, q, p);

    p=0; 
    q=0;
    for (i=0; i<nvar; i++) {
        if (stmt->is_orig_loop[i])    {
            q=0;
            for (j=0; j<stmt->trans->nrows; j++) {
                /* Skip rows of h that are zero */
                if (hProps[j].type != H_SCALAR)   {
                    h = isl_mat_set_element_si(h, q, p, stmt->trans->val[j][i]);
                    q++;
                }
            }
            p++;
        }
    }

    h = isl_mat_right_kernel(h);

    PlutoMatrix *ortho = pluto_matrix_from_isl_mat(h);

    isl_mat_free(h);

    orthcst = (PlutoConstraints **) malloc((nvar+1)*sizeof(PlutoConstraints *)); 

    for (i=0; i<nvar+1; i++)  {
        orthcst[i] = pluto_constraints_alloc(1, CST_WIDTH);
        orthcst[i]->ncols = CST_WIDTH;
    }

    /* All non-negative orthant only */
    /* An optimized version where the constraints are added as
     * c_1 >= 0, c_2 >= 0, ..., c_n >= 0, c_1+c_2+..+c_n >= 1
     *
     * basically only look in the orthogonal space where everything is
     * non-negative
     *
     * All of these constraints are added later to 
     * the global constraint matrix
     */

    /* Normalize ortho first */
    for (j=0; j<ortho->ncols; j++)    {
        if (ortho->val[0][j] == 0) continue;
        int colgcd = abs(ortho->val[0][j]);
        for (i=1; i<ortho->nrows; i++)    {
            if (ortho->val[i][j] == 0)  break;
            colgcd = gcd(colgcd,abs(ortho->val[i][j]));
        }
        if (i == ortho->nrows)   {
            if (colgcd > 1)    {
                for (k=0; k<ortho->nrows; k++)    {
                    ortho->val[k][j] /= colgcd;
                }
            }
        }
    }
    // printf("Ortho matrix\n");
    // pluto_matrix_print(stdout, ortho); 

    isl_currcst = isl_basic_set_from_pluto_constraints(ctx, currcst);

    assert(p == ortho->nrows);
    p=0;
    for (i=0; i<ortho->ncols; i++) {
        isl_basic_set *orthcst_i;

        j=0;
        for (q=0; q<nvar; q++) {
            if (stmt->is_orig_loop[q])    {
                orthcst[p]->val[0][npar+1+(stmt->id)*(nvar+1)+q] = ortho->val[j][i];
                j++;
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = -1;
        orthcst_i = isl_basic_set_from_pluto_constraints(ctx, orthcst[p]);
        orthcst[p]->val[0][CST_WIDTH-1] = 0;

        orthcst_i = isl_basic_set_intersect(orthcst_i,
                isl_basic_set_copy(isl_currcst));
        if (isl_basic_set_fast_is_empty(orthcst_i) 
                || isl_basic_set_is_empty(orthcst_i)) {
            pluto_constraints_negate_row(orthcst[p], 0);
        }
        isl_basic_set_free(orthcst_i);
        p++;
        /* assert(p<=nvar-1); */
    }

    // pluto_matrix_print(stdout, stmt->trans);

    if (p > 0)  {
        /* Sum of all of the above is the last constraint */
        for(j=0; j<CST_WIDTH; j++)  {
            for (i=0; i<p; i++) {
                orthcst[p]->val[0][j] += orthcst[i]->val[0][j];
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = -1;
        p++;
    }

    *orthonum = p;

    IF_DEBUG2(printf("Ortho constraints for S%d; %d sets\n", stmt->id+1, *orthonum));
    for (i=0; i<*orthonum; i++) {
        // print_polylib_visual_sets("li", orthcst[i]);
        // IF_DEBUG2(pluto_constraints_print(stdout, orthcst[i]));
    }

    /* Free the unnecessary ones */
    for (i=p; i<nvar+1; i++)    {
        pluto_constraints_free(orthcst[i]);
    }

    pluto_matrix_free(ortho);
    isl_basic_set_free(isl_currcst);
    isl_ctx_free(ctx);

    return orthcst;
}


/*
 * Check whether the dependence is satisfied at level 'level'
 * (works whether the dep is const or non-const, inter-stmt or
 * self edge
 */
bool dep_satisfaction_test(Dep *dep, PlutoProg *prog, int level)
{
    PlutoConstraints *cst;
    int j, src_dim, dest_dim, npar;
    bool is_empty;

    npar = prog->npar;

    Stmt *src_stmt = prog->stmts[dep->src];
    Stmt *dest_stmt = prog->stmts[dep->dest];

    src_dim = src_stmt->dim;
    dest_dim = dest_stmt->dim;

    assert(src_stmt->trans != NULL);
    assert(dest_stmt->trans != NULL);
    assert(level < src_stmt->trans->nrows);
    assert(level < dest_stmt->trans->nrows);

    cst = pluto_constraints_alloc(2*(1+dep->dpolytope->nrows), 
            src_dim+dest_dim+npar+1);

    /*
     * constraint format
     * \phi(src) - \phi (dest) >= 0
     * (reverse of satisfaction)
     */

    cst->is_eq[0] = 0;
    for (j=0; j<src_dim; j++)    {
        cst->val[0][j] = src_stmt->trans->val[level][j];
    }
    for (j=src_dim; j<src_dim+dest_dim; j++)    {
        cst->val[0][j] = -dest_stmt->trans->val[level][j-src_dim];
    }
    for (j=src_dim+dest_dim; j<src_dim+dest_dim+npar+1; j++)    {
        cst->val[0][j] = 
            src_stmt->trans->val[level][j-dest_dim] - dest_stmt->trans->val[level][j-src_dim];
    }

    cst->nrows = 1;

    pluto_constraints_add(cst, dep->dpolytope);

    /* if no solution exists, the dependence is satisfied, i.e., no points
     * satisfy \phi(src) - \phi(dest) <= 0 */
    is_empty = pluto_constraints_is_empty(cst);
    pluto_constraints_free(cst);

    return is_empty;
}

/* Direction vector component at level 'level'
 */
DepDir get_dep_direction(const Dep *dep, const PlutoProg *prog, int level)
{
    PlutoConstraints *cst;
    int j, src, dest;

    int npar = prog->npar;
    Stmt **stmts = prog->stmts;

    src = dep->src;
    dest = dep->dest;

    Stmt *src_stmt = stmts[dep->src];
    Stmt *dest_stmt = stmts[dep->dest];

    int src_dim = src_stmt->dim;
    int dest_dim = dest_stmt->dim;

    assert(level < stmts[src]->trans->nrows);
    assert(level < stmts[dest]->trans->nrows);

    cst = pluto_constraints_alloc(2 * (2 + dep->dpolytope->nrows),
            (src_dim + dest_dim) + npar + 1);

    /*
     * Check for zero
     *
     * To test \phi (dest) - \phi(src) = 0, we try
     *
     * \phi(dest) - \phi(src) >= 1
     */
    cst->is_eq[0] = 0;
    for (j = 0; j < src_dim; j++) {
        cst->val[0][j] = -stmts[src]->trans->val[level][j];
    }
    for (j = src_dim; j < src_dim + dest_dim; j++) {
        cst->val[0][j] = stmts[dest]->trans->val[level][j - src_dim];
    }
    for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
        cst->val[0][j] = -stmts[src]->trans->val[level][j - dest_dim] +
            stmts[dest]->trans->val[level][j - src_dim];
    }
    cst->val[0][src_dim + dest_dim + npar] =
        -stmts[src]->trans->val[level][src_dim + npar] +
        stmts[dest]->trans->val[level][dest_dim + npar] - 1;
    cst->nrows = 1;

    pluto_constraints_add(cst, dep->dpolytope);

    bool is_empty = pluto_constraints_is_empty(cst);

    if (is_empty) {
        for (j = 0; j < src_dim; j++) {
            cst->val[0][j] = stmts[src]->trans->val[level][j];
        }
        for (j = src_dim; j < src_dim + dest_dim; j++) {
            cst->val[0][j] = -stmts[dest]->trans->val[level][j - src_dim];
        }
        for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
            cst->val[0][j] = stmts[src]->trans->val[level][j - dest_dim] 
                - stmts[dest]->trans->val[level][j - src_dim];
        }
        cst->val[0][src_dim + dest_dim + npar] =
            stmts[src]->trans->val[level][src_dim + npar] 
            - stmts[dest]->trans->val[level][dest_dim + npar] - 1;
        cst->nrows = 1;

        pluto_constraints_add(cst, dep->dpolytope);

        is_empty = pluto_constraints_is_empty(cst);

        /* If no solution exists, all points satisfy \phi (dest) - \phi (src) = 0 */
        if (is_empty) {
            pluto_constraints_free(cst);
            return DEP_ZERO;
        }
    }

    /*
     * Check for PLUS
     * Constraint format
     * \phi(dest) - \phi (src) <= -1
     * (reverse of plus)
     */

    for (j = 0; j < src_dim; j++) {
        cst->val[0][j] = stmts[src]->trans->val[level][j];
    }
    for (j = src_dim; j < src_dim + dest_dim; j++) {
        cst->val[0][j] = -stmts[dest]->trans->val[level][j - src_dim];
    }
    for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
        cst->val[0][j] = stmts[src]->trans->val[level][j - dest_dim] 
            - stmts[dest]->trans->val[level][j - src_dim];
    }
    cst->val[0][src_dim + dest_dim + npar] =
        stmts[src]->trans->val[level][src_dim + npar] -
        stmts[dest]->trans->val[level][dest_dim + npar] - 1;

    cst->nrows = 1;

    pluto_constraints_add(cst, dep->dpolytope);

    is_empty = pluto_constraints_is_empty(cst);

    if (is_empty) {
        pluto_constraints_free(cst);
        return DEP_PLUS;
    }

    /*
     * Check for MINUS
     *
     * Constraint format
     * \phi(dest) - \phi (src) >= 1
     * reverse of minus, we alraedy know that it's not zero
     */

    for (j = 0; j < src_dim; j++) {
        cst->val[0][j] = -stmts[src]->trans->val[level][j];
    }
    for (j = src_dim; j < src_dim + dest_dim; j++) {
        cst->val[0][j] = stmts[dest]->trans->val[level][j - src_dim];
    }
    for (j = src_dim + dest_dim; j < src_dim + dest_dim + npar; j++) {
        cst->val[0][j] = -stmts[src]->trans->val[level][j - dest_dim] 
            + stmts[dest]->trans->val[level][j - src_dim];
    }
    cst->val[0][src_dim + dest_dim + npar] =
        -stmts[src]->trans->val[level][src_dim + npar] +
        stmts[dest]->trans->val[level][dest_dim + npar] - 1;
    cst->nrows = 1;

    pluto_constraints_add(cst, dep->dpolytope);

    is_empty = pluto_constraints_is_empty(cst);
    pluto_constraints_free(cst);

    if (is_empty) {
        return DEP_MINUS;
    }

    /* Neither ZERO, nor PLUS, nor MINUS, has to be STAR */
    return DEP_STAR;
}
