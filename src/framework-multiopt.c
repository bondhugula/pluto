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

#define CONSTRAINTS_SIMPLIFY_THRESHOLD 10000
#define MAX_FARKAS_CST  2000

static void multiopt_compute_permutability_constraints_dep(Dep *dep, PlutoProg *prog)
{
    PlutoConstraints *cst, *tiling_valid_cst, *bounding_func_cst;
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
         *          * so nothing for 2*nvar+npar */
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
            /* No parametric shift coefficients */
        }
        /* Statement's translation coefficients cancel out */

        /* for w */
        phi->val[2*nvar+npar][npar] = 1;
    }

    /* Apply Farkas lemma for bounding function constraints */
    bounding_func_cst = farkas_lemma_affine(dep->bounding_poly, phi);

    pluto_matrix_free(phi);
    pluto_constraints_free(dpoly);
    /* Aggregate permutability and bounding function constraints together in
     *      * global format; note that tiling_valid_cst and bounding_func_cst are 
     *           * local to a  dependence/statements pertaining to it) */

    /* Everything initialized to zero during allocation */
    cst = pluto_constraints_alloc(tiling_valid_cst->nrows + bounding_func_cst->nrows, CST_WIDTH);
    cst->nrows = 0;
    cst->ncols = CST_WIDTH;

    src_offset = npar+1+src_stmt*(nvar+1);
    dest_offset = npar+1+dest_stmt*(nvar+1);

    /* Permutability constraints */
    if (!IS_RAR(dep->type)) {
        /* Permutability constraints only for non-RAR deps */
        for (k=0; k<tiling_valid_cst->nrows; k++) {
            pluto_constraints_add_constraint(cst, tiling_valid_cst->is_eq[k]);
            for (j=0; j<nvar+1; j++)  {
                cst->val[cst->nrows-1][src_offset+j] = tiling_valid_cst->val[k][j];
                if (src_stmt != dest_stmt) {
                    cst->val[cst->nrows-1][dest_offset+j] = tiling_valid_cst->val[k][nvar+1+j];
                }
            }
            /* constant part */
            if (src_stmt == dest_stmt) {
                cst->val[cst->nrows-1][cst->ncols-1] = tiling_valid_cst->val[k][nvar+1];
            }else{
                cst->val[cst->nrows-1][cst->ncols-1] = tiling_valid_cst->val[k][2*nvar+2];
            }
        }
    }

    /* Add bounding function constraints */
    if (!options->nodepbound)   {
        /* Bounding function constraints in global format */
        PlutoConstraints *bcst_g;

        src_offset = npar+1+src_stmt*(nvar+1);
        dest_offset = npar+1+dest_stmt*(nvar+1);

        bcst_g = pluto_constraints_alloc(bounding_func_cst->nrows, CST_WIDTH);

        for (k=0; k<bounding_func_cst->nrows; k++)   {
            pluto_constraints_add_constraint(bcst_g, bounding_func_cst->is_eq[k]);
            for (j=0; j<npar+1; j++)  {
                bcst_g->val[bcst_g->nrows-1][j] = bounding_func_cst->val[k][j];
            }
            for (j=0; j<nvar+1; j++)  {
                bcst_g->val[bcst_g->nrows-1][src_offset+j] = bounding_func_cst->val[k][npar+1+j];
                if (src_stmt != dest_stmt) {
                    bcst_g->val[bcst_g->nrows-1][dest_offset+j] = bounding_func_cst->val[k][npar+1+nvar+1+j];
                }
            }
            /* constant part */
            if (src_stmt == dest_stmt) {
                bcst_g->val[bcst_g->nrows-1][bcst_g->ncols-1] = bounding_func_cst->val[k][npar+1+nvar+1];
            }else{
                bcst_g->val[bcst_g->nrows-1][bcst_g->ncols-1] = bounding_func_cst->val[k][npar+1+2*nvar+2];
            }
        }
        pluto_constraints_add(cst, bcst_g);

        pluto_constraints_free(dep->bounding_cst);
        dep->bounding_cst = bcst_g;
    }

    /* Coefficients of those dimensions that were added for padding
     *      * are of no utility */
    for (k=0; k<nvar; k++)    {
        if (!stmts[src_stmt]->is_orig_loop[k]) {
            for (j=0; j < cst->nrows; j++)   {
                cst->val[j][src_offset+k] = 0;
            }
        }
        if (src_stmt != dest_offset && !stmts[dest_stmt]->is_orig_loop[k])  {
            for (j=0; j < cst->nrows; j++)   {
                cst->val[j][dest_offset+k] = 0;
            }
        }
    }

    pluto_constraints_free(dep->cst);
    dep->cst = cst;

    pluto_constraints_free(tiling_valid_cst);
    pluto_constraints_free(bounding_func_cst);
}

/* Each connected component will have a separate value of u and w. This is done by having replacing the current u with a new one corresponding to the CC to which the dep edge belongs */
PlutoConstraints* multiopt_get_permutability_constraints(PlutoProg *prog)
{
    int i, j, k, ndeps, num_ccs, offset, src_id, cc_id;
    int nvar, npar, nstmts;
    Graph *ddg;
    PlutoConstraints* cc_permute_cst, *new_dep_cst;
    Dep **deps;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;

    deps = prog->deps;
    ddg = prog->ddg;
    num_ccs = prog->ddg->num_ccs;
    ndeps = prog->ndeps;
    cc_permute_cst = NULL;
    new_dep_cst = NULL;

    if (prog->globcst == NULL) {
        get_permutability_constraints(prog);
    }
    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];
        /* Note that you need to add code to check for skipdeps here */
        if (options->rar == 0 && IS_RAR(dep->type)) continue;
        if (dep_is_satisfied(dep)) {
            continue;
        }

        /* This is a redundant condition */
        if (dep->cst == NULL) {
            /* First time, compute the constraints */
            multiopt_compute_permutability_constraints_dep(dep, prog);

            IF_DEBUG(fprintf(stdout, "\tFor dep %d; num_constraints: %d\n",
                        i+1, dep->cst->nrows));
        }
        /* This can be ensured by a call to get_permutability_cst. This is done to enusre that we allocate momory for constraints only once. */
        assert(prog->globcst != NULL);
        if(cc_permute_cst == NULL) {
            cc_permute_cst = pluto_constraints_alloc(prog->globcst->nrows+1, (num_ccs-1)*(npar+1) + CST_WIDTH +1);
            cc_permute_cst->nrows = 1;
            cc_permute_cst->ncols = (num_ccs-1) * (npar+1) + CST_WIDTH + 1;

            cc_permute_cst->val[0][0] = 1;
            for (j=0; j<num_ccs; j++){
                for (k=0;k<npar+1; k++){
                    cc_permute_cst->val[0][j*(npar+1)+k+1] = -1;
                }
            }
        }

        src_id = dep->src;
        cc_id = ddg->vertices[src_id].cc_id;
        if(new_dep_cst == NULL) {
            new_dep_cst = pluto_constraints_alloc (dep->cst->nrows, dep->cst->ncols);
        }
        pluto_constraints_copy(new_dep_cst, dep->cst);
        /* new_dep_cst = pluto_constraints_dup(dep->cst); */

        printf("Dependence Constraints before resize\n");
        pluto_constraints_cplex_print(stdout, new_dep_cst);

        pluto_constraints_resize_single(new_dep_cst, dep->cst->nrows,cc_permute_cst->ncols);

        printf("Dependence Constraints after resize\n");
        pluto_constraints_cplex_print(stdout, new_dep_cst);

        offset = (num_ccs-1)*(npar+1)+1;
        for (j=0; j<new_dep_cst->nrows; j++) {
            for(k=CST_WIDTH-1; k>=npar+1; k--){
                new_dep_cst->val[j][offset+k] = new_dep_cst->val[j][k];
                new_dep_cst->val[j][k] = 0;
            }
        }

        offset = cc_id*(npar+1)+1;
        for (j=0; j<new_dep_cst->nrows; j++){
            for(k=0; k<npar+1; k++){
                new_dep_cst->val[j][offset+k] = new_dep_cst->val[j][k];
                new_dep_cst->val[j][k]=0;
            }
        }

        pluto_constraints_add(cc_permute_cst, new_dep_cst);
        printf("Dep cst \n");
        pluto_constraints_cplex_print(stdout, dep->cst);
        printf(" Resized Dep cst \n");
        pluto_constraints_cplex_print(stdout, new_dep_cst);
        pluto_constraints_free(new_dep_cst);
    }

    /* There are no unsatisfied deps */
    if (cc_permute_cst == NULL) {
        cc_permute_cst = pluto_constraints_alloc(0, (num_ccs-1) * (npar+1) + CST_WIDTH +1);
    }
    return cc_permute_cst;
}

void multiopt_add_stmt_hyperplane_from_ilp_solutions (int *bestsol, PlutoProg *prog) {
    int j, k, nstmts, nvar, npar, num_ccs;
    Stmt **stmts;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    num_ccs = prog->ddg->num_ccs;

    stmts = prog->stmts;

    for (j=0; j<nstmts; j++) {
        Stmt *stmt = stmts[j];
        pluto_stmt_add_hyperplane(stmt, H_UNKNOWN, stmt->trans->nrows);
        for (k=0; k<nvar; k++)    {
            stmt->trans->val[stmt->trans->nrows-1][k] =
                bestsol[num_ccs*(npar+1)+1+j*(nvar+1)+k];
        }
        /* No parameteric shifts */
        for (k=nvar; k<nvar+npar; k++) {
            stmt->trans->val[stmt->trans->nrows-1][k] = 0;
        }
        stmt->trans->val[stmt->trans->nrows-1][nvar+npar] =
            bestsol[num_ccs*(npar+1)+1+j*(nvar+1)+nvar];

        stmt->hyp_types[stmt->trans->nrows-1] =
            pluto_is_hyperplane_scalar(stmt, stmt->trans->nrows-1)?
            H_SCALAR: H_LOOP;

    }
    return;

}

void resize_cst_multiopt(PlutoConstraints *cst, PlutoProg* prog)
{
    int nvar, npar, nstmts, num_ccs, offset;
    int i,j;
    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    num_ccs = prog->ddg->num_ccs;

    assert (cst->ncols == CST_WIDTH);
    pluto_constraints_resize_single(cst, cst->nrows, (num_ccs-1)*(npar+1)+CST_WIDTH+1);
    offset = (num_ccs-1)*(npar+1)+1;
    for (i=0; i<cst->nrows; i++){
        for(j=CST_WIDTH-1; j>=npar+1; j--){
            cst->val[i][offset+j] = cst->val[i][j];
            cst->val[i][j] = 0;
        }
        for (j=0; j<num_ccs*(npar+1)+1; j++){
            cst->val[i][j] = 0;
        }
    }
}


