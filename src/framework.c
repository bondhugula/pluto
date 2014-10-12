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

static void eliminate_farkas_multipliers(PlutoConstraints *farkas_cst, int num_elim);



/**
 *
 * Each constraint row is represented as follows
 *
 *      [comm. vol bound | mapping coeff.s for S1, S2,... |constant]
 * Size:[    npar+1      | (nvar+1)*nstmts                | 1      ]
 *
 * npar - number of parameters in whole program
 * nvar - number of parameters in whole program
 *
 */
#if 0
static PlutoConstraints *get_permutability_constraints_uniform_dep (Dep *dep)
{
    int cst_offset;
    int j, dest_stmt;
    PlutoConstraints *cst;

    /* constant dependences */
    /* uniform self-edge, no need to apply farkas */
    dest_stmt = dep->dest;


    cst_offset = npar+1+dest_stmt*(nvar+1);

    cst = constraints_alloc(2, CST_WIDTH);
    cst->ncols = CST_WIDTH;

    if (!IS_RAR(dep->type)) {
        cst->nrows = 2;
        /* Tiling legality constraint */
        for (j=0; j<nvar; j++)  {
            cst->val[0][cst_offset+j] = -dep->h->val[j][nvar+npar];
        }
        /* Translation coefficient */
        cst->val[0][cst_offset+nvar]=0;

        /* Add bounding function */
        for (j=0; j<npar; j++)  {
            cst->val[1][j] = 0;
        }
        cst->val[1][npar] = 1;
        for (j=cst_offset; j<cst_offset+nvar; j++)  {
            cst->val[1][j] = -cst->val[0][j];
        }
        cst->val[1][cst_offset+nvar]=0;
    }else{
        /* Add bounding function */
        for (j=0; j<npar; j++)  {
            cst->val[0][j] = 0;
        }
        cst->val[0][npar] = 1;
        for (j=cst_offset; j<cst_offset+nvar; j++)  {
            cst->val[0][j] = dep->h->val[j-cst_offset][nvar+npar];
        }
        cst->val[0][cst_offset+nvar]=0;
        cst->nrows=1;
    }

    return cst;
}
#endif


/* Builds validity constraints for a non-uniform dependence */
static void compute_permutability_constraints_dep(Dep *dep, PlutoProg *prog)
{
    PlutoConstraints *farkas_cst, *comm_farkas_cst, *cst;
    int src_stmt, dest_stmt, j, k;
    int src_offset, dest_offset;

    int nvar = prog->nvar;
    int npar = prog->npar;
    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    dest_stmt = dep->dest;
    src_stmt = dep->src;

    /* Convert everything to >= 0 form */
    PlutoConstraints *dpoly = pluto_constraints_to_pure_inequalities_single(dep->dpolytope);

    /* Non-uniform dependence - farkas lemma comes in */
    /* Apply farkas lemma, eliminate farkas multipliers using
     * fourier-motzkin 
     * 
     * -- farkas_cst format for legality --
     * [ mapping coeff for src | ... for dest |farkas multipliers|constant]
     * SIZE: [nvar+1 | nvar+1 | dep.dpolytope->nrows+1 | 1]
     *
     * -- farkas_cst format for bounding function --
     * [bounding func | mapping coeff for src | ... for dest |farkas multipliers|constant]
     * SIZE: [npar+1| nvar+1 | nvar+1 | dep.dpolytope->nrows+1 | 1]
     *
     */
    if (src_stmt != dest_stmt)  {
        /* Inter-statement non-uniform dep */
        farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, 2*nvar+2+dpoly->nrows+2);
        comm_farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, npar+1+2*nvar+2+dpoly->nrows+2);

        farkas_cst->nrows = (2*nvar+npar+1)+1+dpoly->nrows+1;
        farkas_cst->ncols = 2*(nvar+1)+dpoly->nrows+2;

        comm_farkas_cst->nrows = (2*nvar+npar+1)+1+dpoly->nrows+1;
        comm_farkas_cst->ncols = npar+1+2*(nvar+1)+dpoly->nrows+2;
    }else{
        /* Intra-statement non-uniform dependence */
        farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, nvar+1+dpoly->nrows+2);
        comm_farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, npar+1+nvar+1+dpoly->nrows+2);

        farkas_cst->nrows = (2*nvar+npar+1)+1+dpoly->nrows+1;
        farkas_cst->ncols = (nvar+1)+dpoly->nrows+2;

        comm_farkas_cst->nrows = (2*nvar+npar+1)+1+dpoly->nrows+1;
        comm_farkas_cst->ncols = npar+1+(nvar+1)+dpoly->nrows+2;
    }


    /* Initialize all to zero */
    for (j=0; j<farkas_cst->nrows; j++)  {
        for (k=0; k<farkas_cst->ncols; k++)  {
            farkas_cst->val[j][k] = 0;
        }
    }

    for (j=0; j<comm_farkas_cst->nrows; j++)  {
        for (k=0; k<comm_farkas_cst->ncols; k++)  {
            comm_farkas_cst->val[j][k] = 0;
        }
    }

    if (src_stmt != dest_stmt)  {

        /* Add tiling legality constraints */
        for (j=0; j<2*nvar+npar+1; j++)  {
            if (j < nvar)   {
                /* src stmt coeff */
                farkas_cst->val[j][j] = -1;
            }else if (j < 2*nvar)   {
                /* dest stmt coeff */
                farkas_cst->val[j][j+1] = 1;
            }else if (j < 2*nvar+npar)  {
                /* Do nothing - all coeff multipliers stay zero */
                /* since structure parameters not in our affine mappings */
            }else{
                /* j = 2*nvar+npar */
                /* Translation coefficients in the affine mappings */
                farkas_cst->val[j][nvar] = -1;
                farkas_cst->val[j][2*nvar+1] = 1;
                /* \lambda_0 */
                farkas_cst->val[j][farkas_cst->ncols-2] = -1;
            } 

            /* Set coeff's for farkas multipliers (all except \lambda_0) */
            for (k=2*nvar+2; k<2*nvar+2+dpoly->nrows; k++)  {
                /* Note that dep polytope is dpolytope->nrows x (2*nvar+npar+1) */
                farkas_cst->val[j][k] = -dpoly->val[k-2*nvar-2][j];
            }
            farkas_cst->val[j][farkas_cst->ncols-1] = 0;
        }

        /* Since the above are equalities - add sigma negative */
        for (k=0; k<farkas_cst->ncols; k++)    {
            farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++)  {
                farkas_cst->val[2*nvar+npar+1][k] -= farkas_cst->val[j][k];
            }
        }

        /* All Farkas multipliers are non-negative */
        for (j=0; j<dpoly->nrows+1; j++)  {
            for (k=0; k<dpoly->nrows+1; k++)  {
                farkas_cst->val[2*nvar+npar+2+j][2*nvar+2+k] = ((j==k)?1:0);
            }
        }

        /* Bounding function constraints */
        for (k=0; k<npar; k++)  {
            comm_farkas_cst->val[2*nvar+k][k] = 1;
        }

        comm_farkas_cst->val[2*nvar+npar][npar] = 1;

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=0; k<farkas_cst->ncols-dpoly->nrows-2; k++)
                comm_farkas_cst->val[j][npar+1+k] = -farkas_cst->val[j][k];

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=farkas_cst->ncols-dpoly->nrows-2; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];

        /* Add opp inequality since the above were equalities */
        for (k=0; k<comm_farkas_cst->ncols; k++)    {
            comm_farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++) {
                comm_farkas_cst->val[2*nvar+npar+1][k] -= comm_farkas_cst->val[j][k];
            }
        }

        for (j=2*nvar+npar+2; j<farkas_cst->nrows; j++)
            for (k=0; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];
        
        eliminate_farkas_multipliers(farkas_cst, farkas_cst->ncols-2*nvar-3);
        eliminate_farkas_multipliers(comm_farkas_cst, comm_farkas_cst->ncols-npar-1-2*nvar-3);

        /* constraints_print(stdout, farkas_cst); */

    }else{
        /* Source stmt == Dest stmt */

        for (j=0; j<2*nvar+npar+1; j++)  {
            if (j < nvar)   {
                /* src stmt coeff */
                farkas_cst->val[j][j] = -1;
            }else if (j < 2*nvar)   {
                /* dest stmt coeff */
                farkas_cst->val[j][j-nvar] = 1;
            }else if (j < 2*nvar+npar)  {
                /* Do nothing - all coeff multipliers stay zero */
                /* NOTE: structure parameters not in our affine mappings */
            }else{
                /* Translation coefficient gets subtracted out */
                farkas_cst->val[j][nvar] = 0;
                farkas_cst->val[j][farkas_cst->ncols-2] = -1;
            } 

            /* Set coeff's for farkas multipliers */
            for (k=nvar+1; k<nvar+1+dpoly->nrows; k++)  {
                farkas_cst->val[j][k] = -dpoly->val[k-nvar-1][j];
            }
            farkas_cst->val[j][farkas_cst->ncols-1] = 0;
        }

        /* Since the above are equalities - add sigma negative */
        for (k=0; k<farkas_cst->ncols; k++)    {
            farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++)  {
                farkas_cst->val[2*nvar+npar+1][k] -= farkas_cst->val[j][k];
            }
        }

        /* All farkas multipliers are positive */
        for (j=0; j<dpoly->nrows+1; j++)  {
            for (k=0; k<dpoly->nrows+1; k++)  {
                farkas_cst->val[2*nvar+npar+2+j][nvar+1+k] = ((j==k)?1:0);
            }
        }

        /* Bounding function constraints */
        for (k=0; k<npar; k++)  {
            comm_farkas_cst->val[2*nvar+k][k] = 1;
        }

        comm_farkas_cst->val[2*nvar+npar][npar] = 1;

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=0; k<farkas_cst->ncols-dpoly->nrows-2; k++)
                comm_farkas_cst->val[j][npar+1+k] = -farkas_cst->val[j][k];

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=farkas_cst->ncols-dpoly->nrows-2; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];

        /* Add opp inequality since the above were equalities */
        for (k=0; k<comm_farkas_cst->ncols; k++)    {
            comm_farkas_cst->val[2*nvar+npar+1][k] = 0;
            for (j=0; j<2*nvar+npar+1; j++) {
                comm_farkas_cst->val[2*nvar+npar+1][k] -= comm_farkas_cst->val[j][k];
            }
        }

        for (j=2*nvar+npar+2; j<farkas_cst->nrows; j++)
            for (k=0; k<farkas_cst->ncols; k++)
                comm_farkas_cst->val[j][npar+1+k] = farkas_cst->val[j][k];

        eliminate_farkas_multipliers(farkas_cst, farkas_cst->ncols-nvar-2);
        eliminate_farkas_multipliers(comm_farkas_cst, comm_farkas_cst->ncols-npar-1-nvar-2);

        /* constraints_print(stdout, farkas_cst); */
    }

    /* Aggregate permutability and bounding function constraints together in
     * global format format */

    /* Initialize everything to zero */
    cst = pluto_constraints_alloc(farkas_cst->nrows + comm_farkas_cst->nrows, CST_WIDTH);
    cst->ncols = CST_WIDTH;

    for (k=0; k<farkas_cst->nrows+comm_farkas_cst->nrows; k++)   {
        for (j=0; j<cst->ncols; j++)  {
            cst->val[k][j] = 0;
        }
    }

    src_offset = npar+1+src_stmt*(nvar+1);
    dest_offset = npar+1+dest_stmt*(nvar+1);

    /* Permutability constraints */
    if (!IS_RAR(dep->type)) {
        for (k=0; k<farkas_cst->nrows; k++)   {
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows+k][src_offset+j] = farkas_cst->val[k][j];
                if (src_stmt != dest_stmt) {
                    cst->val[cst->nrows+k][dest_offset+j] = farkas_cst->val[k][nvar+1+j];
                }
            }
            /* constant part */
            if (src_stmt == dest_stmt)  {
                cst->val[cst->nrows+k][cst->ncols-1] = farkas_cst->val[k][nvar+1];
            }else{
                cst->val[cst->nrows+k][cst->ncols-1] = farkas_cst->val[k][2*nvar+2];
            }
        }
        cst->nrows = farkas_cst->nrows;
    }

    if (!options->nobound)   {
        /* Add bounding constraints */
        src_offset = npar+1+src_stmt*(nvar+1);
        dest_offset = npar+1+dest_stmt*(nvar+1);

        for (k=0; k<comm_farkas_cst->nrows; k++)   {
            for (j=0; j<npar+1; j++)  {
                cst->val[cst->nrows+k][j] = comm_farkas_cst->val[k][j];
            }
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows+k][src_offset+j] = comm_farkas_cst->val[k][npar+1+j];
                if (src_stmt != dest_stmt) cst->val[cst->nrows+k][dest_offset+j] = comm_farkas_cst->val[k][npar+1+nvar+1+j];
            }
            /* constant part */
            if (src_stmt == dest_stmt)  {
                cst->val[cst->nrows+k][cst->ncols-1] = comm_farkas_cst->val[k][npar+1+nvar+1];
            }else{
                cst->val[cst->nrows+k][cst->ncols-1] = comm_farkas_cst->val[k][npar+1+2*nvar+2];
            }
        }
        cst->nrows += comm_farkas_cst->nrows;
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
            for (j=0; j < farkas_cst->nrows+comm_farkas_cst->nrows; j++)   {
                cst->val[j][dest_offset+k] = 0;
            }
        }
    }

    PlutoConstraints *bounding_cst = NULL;

    /* Copy only the bounding constraints */
    if (!options->nobound) {

		bounding_cst = pluto_constraints_alloc(comm_farkas_cst->nrows, CST_WIDTH);
		bounding_cst->ncols = CST_WIDTH;
		bounding_cst->nrows = comm_farkas_cst->nrows;

		for (k=0; k<comm_farkas_cst->nrows; k++)   {
			for (j=0; j<(bounding_cst)->ncols; j++)  {
				(bounding_cst)->val[k][j] = 0;
			}
		}

		assert(cst->ncols == bounding_cst->ncols);

		for (k=0; k<comm_farkas_cst->nrows; k++)   {
			for (j=0; j<bounding_cst->ncols; j++)  {
				bounding_cst->val[k][j] = cst->val[farkas_cst->nrows+k][j];
			}
		}

    }

    pluto_constraints_free(farkas_cst);
    pluto_constraints_free(comm_farkas_cst);
    pluto_constraints_free(dpoly);

    free(dep->valid_cst);
    dep->valid_cst = cst;

    free(dep->bounding_cst);
    dep->bounding_cst = bounding_cst;
}


/* This function itself is NOT thread-safe for the same PlutoProg */
PlutoConstraints *get_permutability_constraints(Dep **deps, int ndeps, 
        PlutoProg *prog)
{
    int i, nstmts, nvar, npar;
    PlutoConstraints *globcst ;

    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;
    globcst = prog->globcst;

    int total_cst_rows = 0;

    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

        if (dep->valid_cst == NULL) {
            /* First time, get the constraints */
            compute_permutability_constraints_dep(dep, prog);

            IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", 
                        i+1, dep->valid_cst->nrows));
            total_cst_rows += dep->valid_cst->nrows;
            /* IF_DEBUG(fprintf(stdout, "Constraints for dep: %d\n", i+1)); */
            /* IF_DEBUG(pluto_constraints_pretty_print(stdout, depcst[i])); */
        }
    }

    if (!globcst) {
        globcst = pluto_constraints_alloc(total_cst_rows, CST_WIDTH);
    }
    globcst->ncols = CST_WIDTH;
    globcst->nrows = 0;

    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];

		print_polylib_visual_sets("BB_cst", dep->bounding_cst);

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

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
        print_polylib_visual_sets("global", dep->valid_cst);

        IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", i+ 1,
                    globcst->nrows));
        if (globcst->nrows >= 0.7 * MAX_CONSTRAINTS) {
            IF_DEBUG(fprintf(stdout,
                        "After dep: %d; num_constraints_simplified: %d\n", i + 1,
                        globcst->nrows));
        }
        pluto_constraints_simplify(globcst);
        /* pluto_constraints_pretty_print(stdout, globcst); */
    }

    pluto_constraints_simplify(globcst);

    IF_DEBUG(fprintf(stdout, "After all dependences: num constraints: %d, \
                num variables: %d\n", globcst->nrows, globcst->ncols-1));

    return globcst;
}


/*
 * Eliminates the last num_elim variables from farkas_cst -- these are the
 * farkas multipliers
 */
static void eliminate_farkas_multipliers(PlutoConstraints *farkas_cst, int num_elim)
{
    int i;
    int best_elim;

    if (options->moredebug) {
        printf("To start with: %d constraints, %d to be eliminated out of %d vars\n", 
                farkas_cst->nrows, num_elim, farkas_cst->ncols-1);
    }

    for (i=0; i<num_elim; i++)  {
        best_elim = pluto_constraints_best_elim_candidate(farkas_cst, num_elim-i);
        fourier_motzkin_eliminate(farkas_cst, best_elim);
        if (options->moredebug) {
            printf("After elimination of %d variable: %d constraints\n", 
                    num_elim-i, farkas_cst->nrows); 
        }
        // pluto_constraints_print(stdout, farkas_cst);
    }

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
        IF_DEBUG2(pluto_constraints_print(stdout, orthcst[i]));
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
    int64 *sol;
    bool retval;

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
    sol = pluto_constraints_solve(cst,DO_NOT_ALLOW_NEGATIVE_COEFF);
    pluto_constraints_free(cst);

    retval = (sol)? false:true;
    free(sol);

    return retval;
}




/* Direction vector component at level 'level'
 * TODO: assumes no parametric shifts 
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

    cst = pluto_constraints_alloc(2*(2+dep->dpolytope->nrows), 
            (src_dim+dest_dim)+npar+1);

    /*
     * Check for zero
     *
     * To test \phi (dest) - \phi(src) = 0, we try 
     *
     * \phi(dest) - \phi(src) >= 1
     */
    cst->is_eq[0] = 0;
    for (j=0; j<src_dim; j++)    {
        cst->val[0][j] = -stmts[src]->trans->val[level][j];
    }
    for (j=src_dim; j<src_dim+dest_dim; j++)    {
        cst->val[0][j] = stmts[dest]->trans->val[level][j-src_dim];
    }
    cst->val[0][src_dim+dest_dim+npar] = 
        -stmts[src]->trans->val[level][src_dim+npar] + stmts[dest]->trans->val[level][dest_dim+npar]-1;
    cst->nrows = 1;

    pluto_constraints_add(cst, dep->dpolytope);

    int64 *sol = pluto_constraints_solve(cst,DO_NOT_ALLOW_NEGATIVE_COEFF);

    if (!sol)   {
        for (j=0; j<src_dim; j++)    {
            cst->val[0][j] = stmts[src]->trans->val[level][j];
        }
        for (j=src_dim; j<src_dim+dest_dim; j++)    {
            cst->val[0][j] = -stmts[dest]->trans->val[level][j-src_dim];
        }
        cst->val[0][src_dim+dest_dim+npar] = 
            stmts[src]->trans->val[level][src_dim+npar] - stmts[dest]->trans->val[level][dest_dim+npar]-1;
        cst->nrows=1;

        pluto_constraints_add(cst, dep->dpolytope);

        sol = pluto_constraints_solve(cst,DO_NOT_ALLOW_NEGATIVE_COEFF);

        /* If no solution exists, all points satisfy \phi (dest) - \phi (src) = 0 */
        if (!sol)   {
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

    for (j=0; j<src_dim; j++)    {
        cst->val[0][j] = stmts[src]->trans->val[level][j];
    }
    for (j=src_dim; j<src_dim+dest_dim; j++)    {
        cst->val[0][j] = -stmts[dest]->trans->val[level][j-src_dim];
    }
    cst->val[0][src_dim+dest_dim+npar] = 
        stmts[src]->trans->val[level][src_dim+npar] - stmts[dest]->trans->val[level][dest_dim+npar] -1;

    cst->nrows=1;

    pluto_constraints_add(cst, dep->dpolytope);

    free(sol);
    sol = pluto_constraints_solve(cst,DO_NOT_ALLOW_NEGATIVE_COEFF);

    if (!sol)   {
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

    for (j=0; j<src_dim; j++)    {
        cst->val[0][j] = -stmts[src]->trans->val[level][j];
    }
    for (j=src_dim; j<src_dim+dest_dim; j++)    {
        cst->val[0][j] = stmts[dest]->trans->val[level][j-src_dim];
    }
    cst->val[0][src_dim+dest_dim+npar] = 
        -stmts[src]->trans->val[level][src_dim+npar] + stmts[dest]->trans->val[level][dest_dim+npar] -1;
    cst->nrows=1;

    pluto_constraints_add(cst, dep->dpolytope);

    free(sol);
    sol = pluto_constraints_solve(cst,DO_NOT_ALLOW_NEGATIVE_COEFF);
    pluto_constraints_free(cst);

    if (!sol)   {   
        return DEP_MINUS;
    }

    free(sol);

    /* Neither ZERO, nor PLUS, nor MINUS, has to be STAR */
    return DEP_STAR;
}
