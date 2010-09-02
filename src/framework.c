/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007--2008 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
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
#include <math.h>
#include <assert.h>

#include "math_support.h"
#include "constraints.h"
#include "pluto.h"

#include "polylib/polylib64.h"
#include "candl/candl.h"

static void eliminate_farkas_multipliers(PlutoConstraints *farkas_cst, int num_elim);
static PlutoMatrix *get_orthogonal_subspace(Matrix *h);

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


/* Builds legality constraints for a non-uniform dependence */
static PlutoConstraints *get_permutability_constraints_nonuniform_dep(Dep *dep, PlutoProg *prog)
{
    PlutoConstraints *farkas_cst, *comm_farkas_cst, *cst;
    int src_stmt, dest_stmt, j, k;
    int src_offset, dest_offset;

    int nvar = prog->nvar;
    int npar = prog->npar;
    Stmt *stmts = prog->stmts;
    int nstmts = prog->nstmts;

    dest_stmt = dep->dest;
    src_stmt = dep->src;

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
        farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, 2*nvar+2+dep->dpolytope->nrows+2);
        comm_farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, npar+1+2*nvar+2+dep->dpolytope->nrows+2);

        farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        farkas_cst->ncols = 2*(nvar+1)+dep->dpolytope->nrows+2;

        comm_farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        comm_farkas_cst->ncols = npar+1+2*(nvar+1)+dep->dpolytope->nrows+2;
    }else{
        /* Intra-statement non-uniform dependence */
        farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, nvar+1+dep->dpolytope->nrows+2);
        comm_farkas_cst = pluto_constraints_alloc(MAX_FARKAS_CST, npar+1+nvar+1+dep->dpolytope->nrows+2);

        farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        farkas_cst->ncols = (nvar+1)+dep->dpolytope->nrows+2;

        comm_farkas_cst->nrows = (2*nvar+npar+1)+1+dep->dpolytope->nrows+1;
        comm_farkas_cst->ncols = npar+1+(nvar+1)+dep->dpolytope->nrows+2;
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
            for (k=2*nvar+2; k<2*nvar+2+dep->dpolytope->nrows; k++)  {
                /* Note that dep polytope is dpolytope->nrows x (2*nvar+npar+1) */
                farkas_cst->val[j][k] = -dep->dpolytope->val[k-2*nvar-2][j];
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
        for (j=0; j<dep->dpolytope->nrows+1; j++)  {
            for (k=0; k<dep->dpolytope->nrows+1; k++)  {
                farkas_cst->val[2*nvar+npar+2+j][2*nvar+2+k] = ((j==k)?1:0);
            }
        }

        /* Bounding function constraints */
        for (k=0; k<npar; k++)  {
            comm_farkas_cst->val[2*nvar+k][k] = 1;
        }

        comm_farkas_cst->val[2*nvar+npar][npar] = 1;

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=0; k<farkas_cst->ncols-dep->dpolytope->nrows-2; k++)
                comm_farkas_cst->val[j][npar+1+k] = -farkas_cst->val[j][k];

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=farkas_cst->ncols-dep->dpolytope->nrows-2; k<farkas_cst->ncols; k++)
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
            for (k=nvar+1; k<nvar+1+dep->dpolytope->nrows; k++)  {
                farkas_cst->val[j][k] = -dep->dpolytope->val[k-nvar-1][j];
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
        for (j=0; j<dep->dpolytope->nrows+1; j++)  {
            for (k=0; k<dep->dpolytope->nrows+1; k++)  {
                farkas_cst->val[2*nvar+npar+2+j][nvar+1+k] = ((j==k)?1:0);
            }
        }

        /* Bounding function constraints */
        for (k=0; k<npar; k++)  {
            comm_farkas_cst->val[2*nvar+k][k] = 1;
        }

        comm_farkas_cst->val[2*nvar+npar][npar] = 1;

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=0; k<farkas_cst->ncols-dep->dpolytope->nrows-2; k++)
                comm_farkas_cst->val[j][npar+1+k] = -farkas_cst->val[j][k];

        for (j=0; j<2*nvar+npar+1; j++)
            for (k=farkas_cst->ncols-dep->dpolytope->nrows-2; k<farkas_cst->ncols; k++)
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
        if (!stmts[src_stmt].is_outer_loop[k])  {
            for (j=0; j < cst->nrows; j++)   {
                cst->val[j][src_offset+k] = 0;
            }
        }
        if (src_stmt != dest_offset && !stmts[dest_stmt].is_outer_loop[k])  {
            for (j=0; j < farkas_cst->nrows+comm_farkas_cst->nrows; j++)   {
                cst->val[j][dest_offset+k] = 0;
            }
        }
    }

    pluto_constraints_free(farkas_cst);
    pluto_constraints_free(comm_farkas_cst);

    return cst;
}


PlutoConstraints *get_permutability_constraints(Dep *deps, int ndeps, 
        PlutoProg *prog)
{
    int i, dest_stmt, src_stmt;
    Dep *dep;
    static PlutoConstraints *globcst = NULL;
    static PlutoConstraints **depcst = NULL;

    int nstmts = prog->nstmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    if (!depcst)   {
        depcst = (PlutoConstraints **) malloc(ndeps*sizeof(PlutoConstraints *));
        for (i=0; i<ndeps; i++) {
            depcst[i] = NULL;
        }
    }

    int total_cst_rows = 0;

// #pragma omp parallel for private(i,dep,dest_stmt,src_stmt) reduction(+:total_cst)
    for (i=0; i<ndeps; i++) {
        dep = &deps[i];

        dest_stmt = dep->dest;
        src_stmt = dep->src;

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

        if (!depcst[i]) {
            /* First time, get the constraints */

            // Candl doesn't separate out uniform depedences and
            // h-transformation
            // if (src_stmt == dest_stmt && IS_UNIFORM(deps[i].type)) {
                /* Uniform self-edge */
                // depcst[i] = get_permutability_constraints_uniform_dep(dep);
            // }else{
                /* Non-uniform dependences */
            depcst[i] = get_permutability_constraints_nonuniform_dep(dep, prog);
            // }

            IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", i+1, depcst[i]->nrows));
            total_cst_rows += depcst[i]->nrows;
        }
    }

    if (!globcst) globcst = pluto_constraints_alloc(total_cst_rows, CST_WIDTH);
    globcst->ncols = CST_WIDTH;
    globcst->nrows = 0;

    for (i=0; i<ndeps; i++) {
        dep = &deps[i];

        if (options->rar == 0 && IS_RAR(dep->type))  {
            continue;
        }

        /* Note that dependences would be marked satisfied (in
         * pluto_auto_transform) only after all possible independent solutions 
         * are found to the formulation
         */ 
        if (dep_is_satisfied(dep)) continue;

        /* Subsequent calls can just use the old ones */
        pluto_constraints_add(globcst, depcst[i]);

        IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints: %d\n", i+1, globcst->nrows));
        if (globcst->nrows >= 0.7*MAX_CONSTRAINTS)  {
            IF_DEBUG(fprintf(stdout, "After dep: %d; num_constraints_simplified: %d\n", i+1, globcst->nrows));
        }
        pluto_constraints_simplify(globcst);
        IF_DEBUG2(pluto_constraints_print(stdout, globcst));
    }

    pluto_constraints_simplify(globcst);

    IF_DEBUG(fprintf(stdout, "After all dependences: num constraints: %d\n", globcst->nrows));

    return globcst;
}


/* PlutoConstraints to avoid trivial solutions (all zeros) */
PlutoConstraints *get_non_trivial_sol_constraints(PlutoProg *prog)
{
    PlutoConstraints *nzcst;
    int i, j, stmt_offset;

    Stmt *stmts = prog->stmts;
    int nstmts = prog->nstmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    nzcst = pluto_constraints_alloc(nstmts, CST_WIDTH);
    nzcst->ncols = CST_WIDTH;

    for (i=0; i<nstmts; i++) {
        /* Don't add the constraint if enough solutions have been found */
        if (stmts[i].num_ind_sols >= stmts[i].dim)   {
            IF_DEBUG2(fprintf(stdout, "non-zero cst: skipping stmt %d\n", i));
            continue;
        }
        stmt_offset = npar+1+i*(nvar+1);
        for (j=0; j<nvar; j++)  {
            if (stmts[i].is_outer_loop[j] == 1)
                nzcst->val[nzcst->nrows][stmt_offset+j] = 1;
        }
        nzcst->val[nzcst->nrows][CST_WIDTH-1] = -1;
        nzcst->nrows++;
    }

    return nzcst;
}


/*
 * Eliminates the last num_elim variables from farkas_cst -- these are the
 * farkas multipliers
 */
static void eliminate_farkas_multipliers(PlutoConstraints *farkas_cst, int num_elim)
{
    int i;
    int best_elim;

    /* printf("To start with: %d constraints, %d to be eliminated out of %d\n", 
            farkas_cst->nrows, num_elim, farkas_cst->ncols-1); */

    for (i=0; i<num_elim; i++)  {
        best_elim = best_elim_candidate(farkas_cst, num_elim-i);
        fourier_motzkin_eliminate(farkas_cst, best_elim);
        /* printf("After elimination of %d variable: %d constraints\n", 
                num_elim-i, farkas_cst->nrows); */
        /* constraints_print(stdout, farkas_cst); */
    }

}


/* Returns linear independence constraints for a single statement */
PlutoConstraints **get_stmt_ortho_constraints(Stmt *stmt, PlutoProg *prog,
        HyperplaneProperties *hProps,  int *orthonum)
{
    int i, j, k, p, q;
    Matrix *h;
    PlutoConstraints **orthcst;

    int nvar = prog->nvar;
    int npar = prog->npar;
    int nstmts = prog->nstmts;

    orthcst = (PlutoConstraints **) malloc(nvar*sizeof(PlutoConstraints *)); 

    for (i=0; i<nvar; i++)  {
        orthcst[i] = pluto_constraints_alloc(1, CST_WIDTH);
        orthcst[i]->ncols = CST_WIDTH;
    }

    if (stmt->num_ind_sols >= stmt->dim) {
        *orthonum = 0;
        return orthcst;
    }

    h = Matrix_Alloc(stmt->trans->nrows, nvar);

    /* Get rid of the variables that don't appear in the domain of this
     * statement and also beta rows
     */
    p=0; 
    q=0;
    for (i=0; i<nvar; i++) {
        if (stmt->is_outer_loop[i])    {
            q=0;
            for (j=0; j<stmt->trans->nrows; j++) {
                /* Skip rows of h that are zero */
                if (hProps[j].type != H_SCALAR)   {
                    h->p[q][p] = stmt->trans->val[j][i];
                    q++;
                }
            }
            p++;
        }
    }

    h->NbRows = q;
    h->NbColumns = p;

    if (h->NbRows == 0) {
        /* no need to add any orthogonality constraints */
        *orthonum = 0;
        Matrix_Free(h);
        return orthcst;
    }

    PlutoMatrix *ortho = get_orthogonal_subspace(h);

    /* Initialize to zero */
    for (k=0; k<nvar; k++) {
        for (i=0; i<1; i++) {
            for (j=0; j<orthcst[k]->ncols; j++) {
                orthcst[k]->val[i][j] = 0;
            }
        }
    }

    /* Positive orthant only */
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
    // pluto_matrix_print(stdout, ortho); 

    p=0;
	assert(h->NbColumns == ortho->nrows);
	assert(h->NbColumns == ortho->ncols);
    for (i=0; i<ortho->ncols; i++) {
        for (j=0; j<ortho->nrows; j++) {
            if (ortho->val[j][i] != 0) break;
        }
        /* Ignore all zero cols */
        if (j==ortho->nrows) continue;

        /* Ignore cols that are -ves of previous ones */
		for (k=0; k<i; k++)	{
            for (j=0; j<ortho->nrows; j++) {
                if (ortho->val[j][i] != -ortho->val[j][k]) break;
            }
            if (j==ortho->nrows) break;
        }
		if (k < i)	continue;

        /* We have a non-zero col */
        j=0;
        for (q=0; q<nvar; q++) {
            if (stmt->is_outer_loop[q])    {
                orthcst[p]->val[0][npar+1+(stmt->id)*(nvar+1)+q] = ortho->val[j][i];
                j++;
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = 0;
        p++;
        assert(p<=nvar-1);
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

#if 0
    /* Since each of the ortho constraints is tried and the
     * best of the solutions will be kept; give all constraints for the 
     * statement
     * */

    p=0;
    for (i=0; i<ncols; i++) {
        for (j=0; j<ncols; j++) {
            if (ortho->val[j][i] != 0) break;
        }
        /* Ignore all zero cols */
        if (j==ncols) continue;

        /* We have a non-zero col */
        j=0;
        for (q=0; q<nvar; q++) {
            if (stmt->is_outer_loop[q])    {
                orthcst[p]->val[0][npar+1+(stmt->id)*(nvar+1)+q] = ortho->val[j][i];
                j++;
            }
        }
        orthcst[p]->nrows = 1;
        orthcst[p]->val[0][CST_WIDTH-1] = -1;
        p++;
    }
#endif

    *orthonum = p;

    /* Free the unnecessary ones */
    for (i=p; i<nvar; i++)    {
        pluto_constraints_free(orthcst[i]);
    }

    /* printf("Ortho constraints: %d\n", *orthonum); */
    // for (i=0; i<*orthonum; i++) {
        // IF_DEBUG2(constraints_print(stdout, orthcst[i]));
    // }

    Matrix_Free(h);
    pluto_matrix_free(ortho);

    return orthcst;
}


/* Given H, Returns I - H^T(HH^T)H which is the subspace orthogonal to 
 * H: i.e., the null space of H */
static PlutoMatrix *get_orthogonal_subspace(Matrix *h)
{
    int nrows, ncols;
    PlutoMatrix *ortho;
    int i, j;
    int scale;
    Matrix *htrans,  *inv, *mat1, *mat2, *mat3, *mat4, *save;

    nrows = h->NbRows;
    ncols = h->NbColumns;

    // print_matrix(stdout, stmt->trans, stmt->trans->nrows, nvar+1);

    htrans = Matrix_Alloc(ncols, nrows);
    mat1 = Matrix_Alloc(nrows, nrows);
    inv = Matrix_Alloc(nrows, nrows);
    mat2 = Matrix_Alloc(nrows, nrows);
    mat3 = Matrix_Alloc(nrows, ncols);
    mat4 = Matrix_Alloc(ncols, ncols);
    save = Matrix_Alloc(nrows, nrows);
    ortho = pluto_matrix_alloc(ncols, ncols);

    Matrix_Free(htrans);
    htrans = Transpose(h);

    // Matrix_Print(stdout, P_VALUE_FMT, h);
    // Matrix_Print(stdout, P_VALUE_FMT, htrans);
    /* compute H.H^T */
    Matrix_Product(h, htrans, mat1);

    /* HACK: polylib seems to be modifying its first argument to inverse;
     * hence saving it */
    for(i=0; i<nrows; i++)
        for(j=0; j<nrows; j++)  {
            save->p[i][j] = mat1->p[i][j];
        }

    Matrix_Inverse(mat1, inv);

    Matrix_Product(inv, save, mat2);
    scale = mat2->p[0][0];

    Matrix_Product(inv, h, mat3);
    Matrix_Product(htrans, mat3, mat4);


    for (i=0; i<ncols; i++) {
        for (j=0; j<ncols; j++)
            ortho->val[i][j] = ((i==j)?scale: 0) - mat4->p[i][j];
    }

    Matrix_Free(htrans);
    Matrix_Free(mat1);
    Matrix_Free(mat2);
    Matrix_Free(mat3);
    Matrix_Free(mat4);
    Matrix_Free(save);
    Matrix_Free(inv);

    return ortho;
}

/*
 * Check whether the dependence is carried at level 'level'
 * (works whether the dep is const or non-const, inter-stmt or
 * self edge
 */
bool dep_satisfaction_test(Dep *dep, PlutoProg *prog, int level)
{
    static PlutoConstraints *cst = NULL;
    int i, j, src, dest, *sol;

    int nvar = prog->nvar;
    int npar = prog->npar;

    Stmt *stmts = prog->stmts;

    src = dep->src;
    dest = dep->dest;

    assert(level < stmts[src].trans->nrows);
    assert(level < stmts[dest].trans->nrows);

    if (!cst || cst->alloc_nrows < 1+dep->dpolytope->nrows)   {
        if (cst) pluto_constraints_free(cst);
        /* rougly allocate twice to prevent frequent increase */
        cst = pluto_constraints_alloc(2*(1+dep->dpolytope->nrows), 2*nvar+npar+1);
    }
    cst->ncols = 2*nvar+npar+1;

    /*
     * constraint format 
     * \phi(src) - \phi (dest) >= 0
     * (reverse of satisfaction)
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = -stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        stmts[src].trans->val[level][nvar] - stmts[dest].trans->val[level][nvar];

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    /* if no solution exists, the dependence is carried, i.e., no points
     * satisfy \geq 0 */ 
    sol = pluto_constraints_solve(cst);

    bool retval = (sol)? false:true;
    free(sol);

    return retval;
}


int get_dep_direction(Dep *dep, PlutoProg *prog, int level)
{
    static PlutoConstraints *cst = NULL;
    int i, j, src, dest;

    int nvar = prog->nvar;
    int npar = prog->npar;
    Stmt *stmts = prog->stmts;

    src = dep->src;
    dest = dep->dest;

    assert(level < stmts[src].trans->nrows);
    assert(level < stmts[dest].trans->nrows);

    if (!cst || cst->alloc_nrows < 2+dep->dpolytope->nrows)   {
        if (cst) pluto_constraints_free(cst);
        /* Rougly allocate twice to prevent frequent increase */
        cst = pluto_constraints_alloc(2*(2+dep->dpolytope->nrows), 2*nvar+npar+1);
    }
    cst->ncols = 2*nvar+npar+1;

    /*
     * Check for zero
     *
     * To test \phi (dest) - \phi(src) = 0, we try 
     *
     * \phi(dest) - \phi(src) >= 1
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = -stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        -stmts[src].trans->val[level][nvar] + stmts[dest].trans->val[level][nvar]-1;

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    int *sol = pluto_constraints_solve(cst);

    if (!sol)   {
        free(sol);

        for (j=0; j<nvar; j++)    {
            cst->val[0][j] = stmts[src].trans->val[level][j];
        }
        for (j=nvar; j<2*nvar; j++)    {
            cst->val[0][j] = -stmts[dest].trans->val[level][j-nvar];
        }
        cst->val[0][2*nvar+npar] = 
            stmts[src].trans->val[level][nvar] - stmts[dest].trans->val[level][nvar]-1;

        for (i=0; i<dep->dpolytope->nrows; i++)  {
            for (j=0; j<2*nvar+npar+1; j++)  {
                cst->val[1+i][j] = dep->dpolytope->val[i][j];
            }
        }

        cst->nrows = 1+dep->dpolytope->nrows;

        sol = pluto_constraints_solve(cst);

        /* If no solution exists, all points satisfy \phi (dest) - \phi (src) = 0 */
        if (!sol)   {
            free(sol);
            return DEP_ZERO;
        }
    }


    /*
     * Check for PLUS
     * Constraint format 
     * \phi(dest) - \phi (src) <= -1
     * (reverse of plus)
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = -stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        stmts[src].trans->val[level][nvar] - stmts[dest].trans->val[level][nvar] -1;

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    sol = pluto_constraints_solve(cst);

    if (!sol)   {
        free(sol);
        return DEP_PLUS;
    }

    /*
     * Check for MINUS
     *
     * Constraint format 
     * \phi(dest) - \phi (src) >= 1
     * reverse of minus, we alraedy know that it's not zero
     */

    for (j=0; j<nvar; j++)    {
        cst->val[0][j] = -stmts[src].trans->val[level][j];
    }
    for (j=nvar; j<2*nvar; j++)    {
        cst->val[0][j] = stmts[dest].trans->val[level][j-nvar];
    }
    cst->val[0][2*nvar+npar] = 
        -stmts[src].trans->val[level][nvar] + stmts[dest].trans->val[level][nvar] -1;

    for (i=0; i<dep->dpolytope->nrows; i++)  {
        for (j=0; j<2*nvar+npar+1; j++)  {
            cst->val[1+i][j] = dep->dpolytope->val[i][j];
        }
    }

    cst->nrows = 1+dep->dpolytope->nrows;

    sol = pluto_constraints_solve(cst);

    if (!sol)   {   
        free(sol);
        return DEP_MINUS;
    }

    /* Neither ZERO, nor PLUS, nor MINUS, has to be STAR */
    return DEP_STAR;
}
