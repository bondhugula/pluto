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
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "pluto.h"
#include "math_support.h"
#include "constraints.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"
#include "ddg.h"
#include "version.h"

/* Iterative search modes */
#define EAGER 0
#define LAZY 1

int dep_satisfaction_update(PlutoProg *prog, int level);
bool dep_satisfaction_test(Dep *dep, PlutoProg *prog, int level);

int get_num_unsatisfied_deps(Dep **deps, int ndeps);
int get_num_unsatisfied_inter_stmt_deps(Dep **deps, int ndeps);

int pluto_diamond_tile(PlutoProg *prog);

/*
 * Returns the number of (new) satisfied dependences at this level
 */
int dep_satisfaction_update(PlutoProg *prog, int level)
{
    int i;
    int num_new_carried;

    int ndeps = prog->ndeps;
    Dep **deps = prog->deps;

    num_new_carried = 0;

    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];
        if (!dep_is_satisfied(dep))   {
            dep->satisfied = dep_satisfaction_test(dep, prog, level);
            if (dep->satisfied)    { 
                if (!IS_RAR(dep->type)) num_new_carried++;
                dep->satisfaction_level = level;
            }
        }
    }

    return num_new_carried;
}


/* Check whether all deps are satisfied */
int deps_satisfaction_check(Dep **deps, int ndeps)
{
    int i;

    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i]->type)) continue;
        if (!dep_is_satisfied(deps[i]))    {
            return false;
        }
    }
    return true;
}


void pluto_compute_dep_satisfaction(PlutoProg *prog)
{
    int i;

    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i]->satisfied =  false;
        prog->deps[i]->satisfaction_level =  prog->num_hyperplanes-1;
    }

    for (i=0; i<prog->num_hyperplanes; i++) {
        dep_satisfaction_update(prog, i);
    }

    /* Create and set satisfaction vectors */
    for (i=0; i<prog->ndeps; i++) {
        int level;
        Dep *dep = prog->deps[i];
        if (IS_RAR(dep->type)) continue;

        /* Dep satisfaction level should have been set */
        assert(dep->satisfaction_level >= 0);
        if (dep->satvec != NULL) free(dep->satvec);
        dep->satvec = (int *) malloc(prog->num_hyperplanes*sizeof(int));

        assert(dep->dirvec != NULL);

        for (level=0; level<prog->num_hyperplanes; level++) {
            if (dep->dirvec[level] != DEP_ZERO && dep->satisfaction_level >= level) {
                dep->satvec[level] = 1;
            }else{
                dep->satvec[level] = 0;
            }
        }

        dep->satvec[dep->satisfaction_level] = 1;
    }
}



bool dep_is_satisfied(Dep *dep)
{
    return dep->satisfied;
}


int num_satisfied_deps (Dep *deps, int ndeps)
{
    int i;

    int num_satisfied = 0;
    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i].type)) continue;
        if (dep_is_satisfied(&deps[i])) num_satisfied++;
    }

    return num_satisfied;
}


int num_inter_stmt_deps (Dep *deps, int ndeps)
{
    int i;
    int count;

    count=0;
    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i].type)) continue;
        if (deps[i].src != deps[i].dest)    {
            count++;
        }
    }
    return count;
}

int num_inter_scc_deps (Stmt *stmts, Dep *deps, int ndeps)
{
    int i, count;

    count=0;
    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i].type)) continue;
        if (dep_is_satisfied(&deps[i])) continue;
        if (stmts[deps[i].src].scc_id != stmts[deps[i].dest].scc_id)
            count++;
    }
    return count;
}


/* PlutoConstraints to avoid trivial solutions (all zeros)
 *
 * loop_search_mode = EAGER: If a statement's transformation is not full-ranked, 
 * a hyperplane, if found, will be a loop hyperplane.
 *                 = LAZY: at least one of the hyperplanes for non-full statements 
 *  should be a loop hyperplane as opposed to all 
 */
PlutoConstraints *get_non_trivial_sol_constraints(const PlutoProg *prog,
        bool loop_search_mode)
{
    PlutoConstraints *nzcst;
    int i, j, stmt_offset, nvar, npar, nstmts;

    Stmt **stmts = prog->stmts;
    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;

    nzcst = pluto_constraints_alloc(nstmts, CST_WIDTH);
    nzcst->ncols = CST_WIDTH;

    if (loop_search_mode == EAGER) {
        for (i=0; i<nstmts; i++) {
            /* Don't add the constraint if enough solutions have been found */
            if (pluto_stmt_get_num_ind_hyps(stmts[i]) >= stmts[i]->dim_orig)   {
                IF_DEBUG2(fprintf(stdout, "non-zero cst: skipping stmt %d\n", i));
                continue;
            }
            stmt_offset = npar+1+i*(nvar+1);
            for (j=0; j<nvar; j++)  {
                if (stmts[i]->is_orig_loop[j] == 1) {
                    nzcst->val[nzcst->nrows][stmt_offset+j] = 1;
                }
            }
            nzcst->val[nzcst->nrows][CST_WIDTH-1] = -1;
            nzcst->nrows++;
        }
    }else{
        assert(loop_search_mode == LAZY);
        for (i=0; i<nstmts; i++) {
            /* Don't add the constraint if enough solutions have been found */
            if (pluto_stmt_get_num_ind_hyps(stmts[i]) >= stmts[i]->dim_orig)   {
                IF_DEBUG2(fprintf(stdout, "non-zero cst: skipping stmt %d\n", i));
                continue;
            }
            stmt_offset = npar+1+i*(nvar+1);
            for (j=0; j<nvar; j++)  {
                if (stmts[i]->is_orig_loop[j] == 1) {
                    nzcst->val[0][stmt_offset+j] = 1;
                }
            }
            nzcst->val[0][CST_WIDTH-1] = -1;
        }
        nzcst->nrows = 1;
    }

    return nzcst;
}

/**
 * Bounds for Pluto ILP variables
 */
static PlutoConstraints *get_coeff_bounding_constraints(PlutoProg *prog)
{
    int i, npar, nstmts, nvar;
    PlutoConstraints *cst;

    npar = prog->npar;
    nstmts = prog->nstmts;
    nvar = prog->nvar;

    cst = pluto_constraints_alloc(1, CST_WIDTH);

    /* Lower bound for bounding coefficients */
    for (i=0; i<npar+1; i++)  {
        pluto_constraints_add_lb(cst, i, 0);
    }
    /* Lower bound for transformation coefficients */
    for (i=0; i<cst->ncols-npar-1-1; i++)  {
        IF_DEBUG2(printf("Adding lower bound %d for transformation coefficients\n", 0););
        pluto_constraints_add_lb(cst, npar+1+i, 0);
    }

    if (options->coeff_bound != -1) {
        for (i=0; i<cst->ncols-npar-1-1; i++)  {
            IF_DEBUG2(printf("Adding upper bound %d for transformation coefficients\n", options->coeff_bound););
            pluto_constraints_add_ub(cst, npar+1+i, options->coeff_bound);
        }
    }else{
        /* Add upper bounds for transformation coefficients */
        int ub = pluto_prog_get_largest_const_in_domains(prog);

        /* Putting too small an upper bound can prevent useful transformations;
         * also, note that an upper bound is added for all statements globally due
         * to the lack of an easy way to determine bounds for each coefficient to
         * prevent spurious transformations that involve shifts proportional to
         * loop bounds
         */
        if (ub >= 10)   {
            for (i=0; i<cst->ncols-npar-1-1; i++)  {
                IF_DEBUG2(printf("Adding upper bound %d for transformation coefficients\n", ub););
                pluto_constraints_add_ub(cst, npar+1+i, ub);
            }
        }
    }

    return cst;
}


/**
 * Coefficient bounds when finding the cone complement; the cone complement
 * could have (and always has in the case of Pluto as opposed to Pluto+)
 * negative coefficients. So, we can't assume non-negative coefficients as in
 * the remaining Pluto hyperplanes
 */

static PlutoConstraints *get_coeff_bounding_constraints_for_cone_complement(PlutoProg *prog)
{
    int i, npar, nstmts, nvar, s;
    PlutoConstraints *cst;

    npar = prog->npar;
    nstmts = prog->nstmts;
    nvar = prog->nvar;

    cst = pluto_constraints_alloc(1, CST_WIDTH);

    /* Lower bound for bounding coefficients */
    for (i=0; i<npar+1; i++)  {
        pluto_constraints_add_lb(cst, i, 0);
    }
    /* Lower bound for transformation coefficients */
    for (s=0; s<nstmts; s++)  {
        for (i=0; i<nvar; i++)  {
            /* Set this to -4 (is enough) */
            IF_DEBUG2(printf("Adding lower bound %d for stmt dim coefficients\n", -4););
            pluto_constraints_add_lb(cst, npar+1+s*(nvar+1)+i, -4);
        }
        IF_DEBUG2(printf("Adding lower bound %d for stmt translation coefficient\n", 0););
        pluto_constraints_add_lb(cst, npar+1+s*(nvar+1)+nvar, 0);
    }
    return cst;
}


/*
 * This calls pluto_constraints_lexmin, but before doing that does some preprocessing
 * - removes variables that we know will be assigned 0 - also do some
 *   permutation/substitution of variables
 */
int64 *pluto_prog_constraints_lexmin(PlutoConstraints *cst, PlutoProg *prog)
{
    Stmt **stmts;
    int nstmts, nvar, npar;

    stmts = prog->stmts;
    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;

    /* Remove redundant variables - that don't appear in your outer loops */
    int redun[npar+1+nstmts*(nvar+1)+1];
    int i, j, k, q;
    int64 *sol, *fsol;
    PlutoConstraints *newcst;

    assert(cst->ncols - 1 == CST_WIDTH - 1);

    for (i=0; i<npar+1; i++)    {
        redun[i] = 0;
    }

    for (i=0; i<nstmts; i++)    {
        for (j=0; j<nvar; j++)    {
            redun[npar+1+i*(nvar+1)+j] = !stmts[i]->is_orig_loop[j];
        }
        redun[npar+1+i*(nvar+1)+nvar] = 0;
    }
    redun[npar+1+nstmts*(nvar+1)] = 0;

    int del_count = 0;
    newcst = pluto_constraints_dup(cst);
    for (j = 0; j < cst->ncols-1; j++) {
        if (redun[j]) {
            pluto_constraints_remove_dim(newcst, j-del_count);
            del_count++;
        }
    }
    IF_DEBUG2(printf("Constraints after reductions\n"));
    IF_DEBUG2(pluto_constraints_compact_print(stdout,newcst));

        
    /* Reverse the variable order for stmts */
    PlutoMatrix *perm_mat = pluto_matrix_alloc(newcst->ncols, newcst->ncols);
    PlutoMatrix *newcstmat = pluto_matrix_alloc(newcst->nrows, newcst->ncols);

    for (i=0; i<newcst->ncols; i++) {
        bzero(perm_mat->val[i], sizeof(int64)*newcst->ncols);
    }

    for (i=0; i<npar+1; i++) {
        perm_mat->val[i][i] = 1;
    }

    j=npar+1;
    for (i=0; i<nstmts; i++)    {
        for (k=j; k<j+stmts[i]->dim_orig; k++) {
            perm_mat->val[k][2*j+stmts[i]->dim_orig-k-1] = 1;
        }
        perm_mat->val[k][k] = 1;
        j += stmts[i]->dim_orig+1;
    }
    perm_mat->val[j][j] = 1;

    for (i=0; i<newcst->nrows; i++) {
        for (j=0; j<newcst->ncols; j++) {
            newcstmat->val[i][j] = 0;
            for (k=0; k<newcst->ncols; k++) {
                newcstmat->val[i][j] += newcst->val[i][k]*perm_mat->val[k][j];
            }
        }
    }
    /* pluto_matrix_print(stdout, newcst->val, newcst->nrows, newcst->ncols); */
    /* pluto_matrix_print(stdout, newcstmat, newcst->nrows, newcst->ncols); */

    PlutoConstraints *newcst_permuted;
    newcst_permuted = pluto_constraints_from_inequalities(newcstmat);
    pluto_matrix_free(newcstmat);

    IF_DEBUG(printf("[Pluto] pluto_prog_constraints_lexmin (%d variables, %d constraints)\n",
                cst->ncols-1, cst->nrows););

    /* Solve the constraints */
    if (options->islsolve) {
        sol = pluto_constraints_lexmin_isl(newcst_permuted, ALLOW_NEGATIVE_COEFF);
        /* print_polylib_visual_sets("csts", newcst); */
    }else{
        sol = pluto_constraints_lexmin_pip(newcst_permuted, ALLOW_NEGATIVE_COEFF);
        /* print_polylib_visual_sets("csts", newcst); */
    }

    pluto_constraints_free(newcst_permuted);

    fsol = NULL;
    if (sol != NULL)    {

        PlutoMatrix *actual_sol = pluto_matrix_alloc(1, newcst->ncols-1);
        for (j=0; j<newcst->ncols-1; j++) {
            actual_sol->val[0][j] = 0;
            for (k=0; k<newcst->ncols-1; k++) {
                actual_sol->val[0][j] += sol[k]*perm_mat->val[k][j];
            }
        }
        free(sol);

        fsol = (int64 *)malloc(cst->ncols*sizeof(int64));
        /* Fill the soln with zeros for the redundant variables */
        q = 0;
        for (j=0; j<cst->ncols-1; j++) {
            if (redun[j])  {
                fsol[j] = 0;
            }else{
                fsol[j] = actual_sol->val[0][q++];
            }
        }
        pluto_matrix_free(actual_sol);
    }

    pluto_matrix_free(perm_mat);
    pluto_constraints_free(newcst);

    return fsol;
}


/* Is there an edge between some vertex of SCC1 and some vertex of SCC2? */
int ddg_sccs_direct_connected(Graph *g, PlutoProg *prog, int scc1, int scc2)
{
    int i, j;

    for (i=0; i<prog->nstmts; i++)  {
        if (prog->stmts[i]->scc_id == scc1)  {
            for (j=0; j<prog->nstmts; j++)  {
                if (prog->stmts[j]->scc_id == scc2)  {
                    if (g->adj->val[i][j] > 0)  {
                        return 1;
                    }
                }
            }
        }
    }

    return 0;
}

/* Cut dependences between two SCCs 
 * Returns: number of dependences cut  */
int cut_between_sccs(PlutoProg *prog, Graph *ddg, int scc1, int scc2)
{
    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    int nvar = prog->nvar;
    int npar = prog->npar;

    int i, j, num_satisfied;

    if (!ddg_sccs_direct_connected(ddg, prog, scc1, scc2))    {
        return 0;
    }

    IF_DEBUG(printf("[pluto] Cutting between SCC id %d and id %d\n", scc1, scc2));

    pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

    for (i=0; i<nstmts; i++) {
        pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
        for (j=0; j<nvar+npar; j++)  {
            stmts[i]->trans->val[stmts[i]->trans->nrows-1][j] = 0;
        }
        if (stmts[i]->scc_id < scc2)   {
            stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar+npar] = 0;
        }else{
            stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar+npar] = 1;
        }

    }
    num_satisfied =  dep_satisfaction_update(prog, stmts[0]->trans->nrows-1);
    if (num_satisfied >= 1) {
        ddg_update(ddg, prog);
    }else{
        for (i=0; i<nstmts; i++) {
            stmts[i]->trans->nrows--;
        }
        prog->num_hyperplanes--;
    }

    return num_satisfied;
}


/*
 * Cut dependences between all SCCs 
 */
int cut_all_sccs(PlutoProg *prog, Graph *ddg)
{
    int i, j, num_satisfied;
    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    IF_DEBUG(printf("[pluto] Cutting between all SCCs\n"));

    if (ddg->num_sccs == 1) {
        IF_DEBUG(printf("\t only one SCC\n"));
        return 0;
    }

    pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

    for (i=0; i<nstmts; i++)    {
        pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
        for (j=0; j<nvar+npar; j++)  {
            stmts[i]->trans->val[stmts[i]->trans->nrows-1][j] = 0;
        }
        stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar+npar] = stmts[i]->scc_id;

    }
    num_satisfied = dep_satisfaction_update(prog, stmts[0]->trans->nrows-1);
    ddg_update(ddg, prog);

    return num_satisfied;
}


/* 
 * Cut based on dimensionalities of SCCs; if two SCCs are of different 
 * dimensionalities; separate them 
 * SCC1 -> SCC2 -> SCC3 ... ->SCCn 
 * Two neighboring SCCs won't be cut if they are of the same
 * dimensionality
 */
int cut_scc_dim_based(PlutoProg *prog, Graph *ddg)
{
    int i, j, k, count;
    Stmt **stmts = prog->stmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    if (ddg->num_sccs == 1) return 0;

    IF_DEBUG(printf("Cutting based on SCC dimensionalities\n"));

    count = 0;

    int cur_max_dim = ddg->sccs[0].max_dim;

    pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

    for (k=0; k<ddg->num_sccs; k++) {
        if (cur_max_dim != ddg->sccs[k].max_dim)   {
            cur_max_dim = ddg->sccs[k].max_dim;
            count++;
        }

        for (i=0; i<prog->nstmts; i++) {
            if (stmts[i]->scc_id == k)  {
                pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
                for (j=0; j<nvar; j++)  {
                    stmts[i]->trans->val[stmts[i]->trans->nrows-1][j] = 0;
                }
                stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar+npar] = count;
            }
        }
    }

    int num_new_carried = dep_satisfaction_update(prog, stmts[0]->trans->nrows-1);

    if (num_new_carried >= 1)   {
        ddg_update(ddg, prog);
    }else{
        for (i=0; i<prog->nstmts; i++) {
            stmts[i]->trans->nrows--;
        }
        prog->num_hyperplanes--;
    }

    return num_new_carried;
}

/* Heuristic cut */
void cut_smart(PlutoProg *prog, Graph *ddg)
{
    if (ddg->num_sccs == 0) return;

    if (pluto_transformations_full_ranked(prog)) {
        /* Enough linearly independent solutions have been found */
        cut_all_sccs(prog, ddg);
        return;
    }

    int i, j;

    int num_new_carried = 0;

    /* First time, cut between SCCs of different dimensionalities */
    if (cut_scc_dim_based(prog, ddg)) {
        return;
    }

    /* Cut in the center */
    if (cut_between_sccs(prog,ddg,ceil(ddg->num_sccs/2.0)-1, 
                ceil(ddg->num_sccs/2.0))) {
        return;
    }

    /* Cut between SCCs that are far away */
    for (i=0; i<ddg->num_sccs-1; i++) {
        for (j=ddg->num_sccs-1; j>=i+1; j--) {
            if (prog->stmts[0]->trans->nrows <= 4*prog->nvar+2)   {
                if (ddg_sccs_direct_connected(ddg, prog, i, j))    {
                    // if (ddg->sccs[i].max_dim == ddg->sccs[j].max_dim) {
                    num_new_carried += cut_between_sccs(prog,ddg,i,j);
                    // }
                }
            }else{
                cut_all_sccs(prog, ddg);
                return;
            }
        }
    }
}


/* Distribute conservatively to maximize (rather random) fusion chance */
void cut_conservative(PlutoProg *prog, Graph *ddg)
{
    int i, j;

    if (cut_scc_dim_based(prog,ddg))   {
        return;
    }

    /* Cut in the center */
    if (cut_between_sccs(prog,ddg,ceil(ddg->num_sccs/2.0)-1,
                ceil(ddg->num_sccs/2.0)))  {
        return;
    }

    /* Cut between SCCs that are far away */
    for (i=0; i<ddg->num_sccs-1; i++) {
        for (j=ddg->num_sccs-1; j>=i+1; j--) {
            if (prog->stmts[0]->trans->nrows <= 4*prog->nvar+2)   {
                if (cut_between_sccs(prog,ddg,i,j)) {
                    return;
                }
            }else{
                cut_all_sccs(prog, ddg);
                return;
            }
        }
    }
}
/*
 * Determine constraints to ensure linear independence of data hyperplanes
 *
 * lin_ind_mode = EAGER: all array hyperplanes have to be linearly independent
 * w.r.t existing ones (ignoring stmts that already have enough lin ind solns)
 *              = LAZY: at least one arry that does not have enough
 * linearly independent solutions will get a new linearly independent
 * hyperplane (this is enough to make progress)
 */
PlutoConstraints *pluto_dist_get_linear_ind_constraints(const PlutoProg *prog,
        const PlutoConstraints *cst, bool lin_ind_mode)
{
    int narrays, i, j, k, orthosum;
    int orthonum[prog->narrays];
    PlutoConstraints ***orthcst;
    Array **arrays = prog->arrays;

    narrays = prog->narrays;


    orthcst = (PlutoConstraints ***) malloc(narrays*sizeof(PlutoConstraints **));

    orthosum = 0;

    /* Get orthogonality constraints for each array*/
    for (j=0; j<narrays; j++)    {
        orthcst[j] = get_array_ortho_constraints(arrays[j],
                prog, cst, &orthonum[j]);
        orthosum += orthonum[j];
    }


    PlutoConstraints *indcst = pluto_constraints_alloc(1, cst->ncols);

    if (orthosum >= 1) {
        if (lin_ind_mode == EAGER) {
            /* Look for linearly independent hyperplanes for all arrays*/
            for (j=0; j<narrays; j++)    {
                if (orthonum[j] >= 1)   {
                    IF_DEBUG2(printf("Added ortho constraints for S%d\n", j+1););
                    pluto_constraints_add(indcst, orthcst[j][orthonum[j]-1]);
                }
            }
        }else{
            assert(lin_ind_mode == LAZY);
            /* At least one array should have a linearly independent hyperplane */

            for (i=0; i<prog->narrays; i++) {
                /* Everything was initialized to zero */
                if (orthonum[i] >= 1) {
                    for (j=0; j<cst->ncols-1; j++) {
                        indcst->val[0][j] += orthcst[i][orthonum[i]-1]->val[0][j];
                    }
                }
            }
            indcst->val[0][cst->ncols-1] = -1;
            indcst->nrows = 1;
            IF_DEBUG2(printf("Added \"at least one\" linear ind constraints\n"););
            IF_DEBUG2(pluto_constraints_pretty_print(stdout, indcst););
        }
    }

    for (j=0; j<narrays; j++)    {
        for (k=0; k<orthonum[j]; k++)   {
            pluto_constraints_free(orthcst[j][k]);
        }
        free(orthcst[j]);
    }
    free(orthcst);

    return indcst;
}

/* 
 * Determine constraints to ensure linear independence of hyperplanes 
 *
 * lin_ind_mode = EAGER: all statement hyperplanes have to be linearly independent
 * w.r.t existing ones (ignoring stmts that already have enough lin ind solns)
 *              = LAZY: at least one statement that does not have enough
 * linearly independent solutions will get a new linearly independent
 * hyperplane (this is enough to make progress)
 */
PlutoConstraints *get_linear_ind_constraints(const PlutoProg *prog, 
        const PlutoConstraints *cst, bool lin_ind_mode)
{
    int npar, nvar, nstmts, i, j, k, orthosum;
    int orthonum[prog->nstmts];
    PlutoConstraints ***orthcst;
    Stmt **stmts;

    IF_DEBUG(printf("[Pluto] get_linear_ind_constraints\n"););

    npar = prog->npar;
    nvar = prog->nvar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;

    orthcst = (PlutoConstraints ***) malloc(nstmts*sizeof(PlutoConstraints **));

    orthosum = 0;

    /* Get orthogonality constraints for each statement */
    for (j=0; j<nstmts; j++)    {
        orthcst[j] = get_stmt_ortho_constraints(stmts[j], 
                prog, cst, &orthonum[j]);
        orthosum += orthonum[j];
    }

    PlutoConstraints *indcst = pluto_constraints_alloc(1, CST_WIDTH);

    if (orthosum >= 1) {
        if (lin_ind_mode == EAGER) {
            /* Look for linearly independent hyperplanes for all stmts */
            for (j=0; j<nstmts; j++)    {
                if (orthonum[j] >= 1)   {
                    IF_DEBUG2(printf("Added ortho constraints for S%d\n", j+1););
                    pluto_constraints_add(indcst, orthcst[j][orthonum[j]-1]);
                }
            }
        }else{
            assert(lin_ind_mode == LAZY);
            /* At least one stmt should have a linearly independent hyperplane */
            for (i=0; i<prog->nstmts; i++) {
                /* Everything was initialized to zero */
                if (orthonum[i] >= 1) {
                    for (j=0; j<CST_WIDTH-1; j++) {
                        indcst->val[0][j] += orthcst[i][orthonum[i]-1]->val[0][j];
                    }
                }
            }
            indcst->val[0][CST_WIDTH-1] = -1;
            indcst->nrows = 1;
            IF_DEBUG2(printf("Added \"at least one\" linear ind constraints\n"););
            IF_DEBUG2(pluto_constraints_pretty_print(stdout, indcst););
        }
    }

    for (j=0; j<nstmts; j++)    {
        for (k=0; k<orthonum[j]; k++)   {
            pluto_constraints_free(orthcst[j][k]);
        }
        free(orthcst[j]);
    }
    free(orthcst);

    return indcst;
}


/* Find all linearly independent permutable band of hyperplanes at a level. 
 *
 * See sub-functions for loop_search_mode and lin_ind_mode
 *
 * If all statements already have enough linearly independent solutions, no
 * independence constraints will be generated, and since no non-trivial
 * solution constraints are added in such a case, the trivial zero solution will 
 * end up being returned.
 * */
int find_permutable_hyperplanes(PlutoProg *prog, bool lin_ind_mode, 
        bool loop_search_mode, int max_sols)
{
    int num_sols_found, j, k;
    int64 *bestsol;
    PlutoConstraints *basecst, *nzcst, *boundcst;
    PlutoConstraints *currcst;

    int nstmts = prog->nstmts;
    Stmt **stmts = prog->stmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    IF_DEBUG(fprintf(stdout, "Finding hyperplanes: max %d\n", max_sols));

    assert(max_sols >= 0);

    if (max_sols == 0) return 0;

    /* Don't free basecst */
    basecst = get_permutability_constraints(prog);
    boundcst = get_coeff_bounding_constraints(prog);
    pluto_constraints_add(basecst, boundcst);
    pluto_constraints_free(boundcst);
    // print_polylib_visual_sets("pluto", basecst);

    num_sols_found = 0;
    /* We don't expect to add a lot to basecst - just ortho constraints
     * and trivial soln avoidance constraints; instead of duplicating basecst,
     * we will just allocate once and copy each time */
    currcst = pluto_constraints_alloc(basecst->nrows+nstmts+nvar*nstmts, 
            CST_WIDTH);

    do{
        pluto_constraints_copy(currcst, basecst);
        nzcst = get_non_trivial_sol_constraints(prog, loop_search_mode);
        pluto_constraints_add(currcst, nzcst);
        pluto_constraints_free(nzcst);
        // print_polylib_visual_sets("curr", currcst);

        PlutoConstraints *indcst = get_linear_ind_constraints(prog, currcst, lin_ind_mode);
        // print_polylib_visual_sets("ind", indcst);

        if (indcst->nrows == 0) {
            /* If you don't have any independence constraints, we would end up finding
             * the same solution that was found earlier; so we won't find anything
             * new */
            bestsol = NULL;
        }else{
            pluto_constraints_add(currcst, indcst);
            IF_DEBUG2(printf("Solving for %d solution\n", num_sols_found+1));
            IF_DEBUG2(pluto_constraints_pretty_print(stdout, currcst));
            bestsol = pluto_prog_constraints_lexmin(currcst, prog);
        }
        pluto_constraints_free(indcst);

        if (bestsol != NULL)    {
            IF_DEBUG(fprintf(stdout, "Found a hyperplane\n"));
            num_sols_found++;

            pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_LOOP);

            for (j=0; j<nstmts; j++)    {
                Stmt *stmt = stmts[j];
                pluto_stmt_add_hyperplane(stmt, H_UNKNOWN, stmt->trans->nrows);
                for (k=0; k<nvar; k++)    {
                    stmt->trans->val[stmt->trans->nrows-1][k] = 
                        bestsol[npar+1+j*(nvar+1)+k];
                }
                /* No parameteric shifts */
                for (k=nvar; k<nvar+npar; k++)    {
                    stmt->trans->val[stmt->trans->nrows-1][k] = 0;
                }
                stmt->trans->val[stmt->trans->nrows-1][nvar+npar] = 
                    bestsol[npar+1+j*(nvar+1)+nvar];

                stmt->hyp_types[stmt->trans->nrows-1] =  
                    pluto_is_hyperplane_scalar(stmt, stmt->trans->nrows-1)?
                    H_SCALAR: H_LOOP;

            }
            free(bestsol);
        }
    }while (num_sols_found < max_sols && bestsol != NULL);

    pluto_constraints_free(currcst);

    /* Same number of solutions are found for each stmt */
    return num_sols_found;
}

/*
 * Returns H_LOOP if this hyperplane is a real loop or H_SCALAR if it's a scalar
 * dimension (beta row or node splitter)
 */
int get_loop_type(Stmt *stmt, int level)
{
    int j;

    for (j=0; j<stmt->trans->ncols-1; j++)    {
        if (stmt->trans->val[level][j] > 0)  {
            return H_LOOP;
        }
    }

    return H_SCALAR;
}


/* Cut based on the .fst file; returns 0 if it fails  */
bool precut(PlutoProg *prog, Graph *ddg, int depth)
{
    int ncomps;

    int nstmts = prog->nstmts;

    int stmtGrp[nstmts][nstmts];
    int grpCount[nstmts];

    int i, j, k;

    if (depth != 0) return false;

    Stmt **stmts = prog->stmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    FILE *cutFp = fopen(".fst", "r");

    if (cutFp)  {
        int tile;

        fscanf(cutFp, "%d", &ncomps);

        if (ncomps > nstmts)   {
            printf("You have an .fst in your directory that is invalid for this source\n");
            printf("No fusion/distribution forced\n");
            return false;
        }

        for (i=0; i<ncomps; i++)    {
            fscanf(cutFp, "%d", &grpCount[i]);
            assert(grpCount[i] <= nstmts);
            for (j=0; j<grpCount[i]; j++)    {
                fscanf(cutFp, "%d", &stmtGrp[i][j]);
                assert(stmtGrp[i][j] <= nstmts-1);
            }
            fscanf(cutFp, "%d", &tile);
            for (j=0; j<grpCount[i]; j++)    
                for (k=0; k<stmts[stmtGrp[i][j]]->dim_orig; k++)
                    stmts[stmtGrp[i][j]]->tile = tile;
        }

        fclose(cutFp);

        /* Update transformation matrices */
        for (i=0; i<nstmts; i++)    {
            pluto_stmt_add_hyperplane(stmts[i], H_SCALAR, stmts[i]->trans->nrows);
        }

        for (i=0; i<ncomps; i++)    {
            for (j=0; j<grpCount[i]; j++)    {
                int id = stmtGrp[i][j];
                for (k=0; k<nvar+npar; k++)  {
                    stmts[id]->trans->val[stmts[id]->trans->nrows-1][k] = 0;
                }
                stmts[id]->trans->val[stmts[id]->trans->nrows-1][nvar+npar] = i;
            }
        }

        pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_SCALAR);

        dep_satisfaction_update(prog, prog->num_hyperplanes-1);
        ddg_update(ddg, prog);

        return true;
    }else{
        FILE *precut = fopen(".precut", "r");
        int ignore, rows, cols, tile, tiling_depth;

        if (precut) {
            /* Num of statements */
            fscanf(precut, "%d", &ignore);

            assert (ignore == prog->nstmts);

            /* Tiling depth */
            fscanf(precut, "%d", &tiling_depth);

            for (i=0; i<prog->nstmts; i++)  {
                /* Read scatterings */
                fscanf(precut, "%d", &rows);
                fscanf(precut, "%d", &cols);

                for (k=0; k<rows; k++)  {
                    /* Transformation is in polylib format
                     * <stmt_orig_dim>+<npar>+1 (first column for equality) */
                    assert(cols == 1+stmts[i]->dim_orig+npar+1);

                    /* For equality - ignore the zero */
                    fscanf(precut, "%d", &ignore);
                    assert(ignore == 0);

                    pluto_matrix_add_row(stmts[i]->trans, stmts[i]->trans->nrows);

                    for (j=0; j<nvar; j++)    {
                        if (stmts[i]->is_orig_loop[j])  {
                            fscanf(precut, "%lld", &stmts[i]->trans->val[stmts[i]->trans->nrows-1][j]);
                        }else{
                            stmts[i]->trans->val[stmts[i]->trans->nrows-1][j] = 0;
                        }
                    }
                    for (j=0; j<npar; j++)    {
                        fscanf(precut, "%d", &ignore);
                        stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar] = 0;
                    }
                    /* Constant part */
                    fscanf(precut, "%lld", &stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar]);

                    // stmts[i]->trans_loop_type[stmts[i]->trans->nrows] = 
                    // (get_loop_type(stmts[i], stmts[i]->trans->nrows) 
                    // == H_SCALAR)? SCALAR:LOOP;
                }

                /* Number of levels */
                fscanf(precut, "%d", &ignore);

                /* FIX this: to tile or not is specified depth-wise, why? Just
                 * specify once */
                for (j=0; j<tiling_depth; j++)    {
                    fscanf(precut, "%d", &tile);
                }
                stmts[i]->tile = tile;
            }

            /* Set hProps correctly and update satisfied dependences */
            for (k=0; k<rows; k++)  {
                pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_UNKNOWN);
                for (i=0; i<nstmts; i++)    {
                    if (get_loop_type(stmts[i], stmts[0]->trans->nrows-rows+k)
                            == H_LOOP)  {
                        stmts[i]->hyp_types[prog->num_hyperplanes-1] = H_LOOP;
                        prog->hProps[prog->num_hyperplanes-1].type = H_LOOP;
                    }else{
                        stmts[i]->hyp_types[prog->num_hyperplanes-1] = H_SCALAR;
                        prog->hProps[prog->num_hyperplanes-1].type = H_SCALAR;
                    }
                }

                dep_satisfaction_update(prog, prog->num_hyperplanes-1);
                ddg_update(ddg, prog);
            }
            return true;
        }

        return false;
    }
}


void pluto_compute_dep_directions(PlutoProg *prog)
{
    int i, level;

    Dep **deps = prog->deps;

    for (i=0; i<prog->ndeps; i++)   {
        if (deps[i]->dirvec != NULL)  {
            free(deps[i]->dirvec);
        }
        deps[i]->dirvec = (DepDir *)malloc(prog->num_hyperplanes*sizeof(DepDir));
        for (level=0; level < prog->num_hyperplanes; level++)  {
            deps[i]->dirvec[level] = get_dep_direction(deps[i], prog, level);
        }
    }
}

void pluto_detect_hyperplane_types_stmtwise(PlutoProg *prog)
{
    int s, i;

    for (s=0; s<prog->nstmts; s++) {
        Stmt *stmt = prog->stmts[s];
        for (i=0; i<stmt->trans->nrows; i++) {
            stmt->hyp_types[i] = pluto_is_hyperplane_loop(stmt, i)? H_LOOP:H_SCALAR;
        }
    }
}

/* Detect H_LOOP or H_SCALAR from scratch */
void pluto_detect_hyperplane_types(PlutoProg *prog)
{
    int i, depth;
    int nstmts = prog->nstmts;

    for (depth=0; depth<prog->num_hyperplanes; depth++) {
        for (i=0; i<nstmts; i++) {
            if (pluto_is_hyperplane_loop(prog->stmts[i], depth)) break;
        }
        prog->hProps[depth].type = (i<nstmts)? H_LOOP: H_SCALAR;
    }
}


/* Detect tilable bands; calculate dependence components (in transformed 
 * space); calls simple dep satisfaction checks */
void pluto_detect_transformation_properties(PlutoProg *prog)
{
    int level, i, j;
    Stmt **stmts = prog->stmts;
    Dep **deps = prog->deps;
    int band, num_loops_in_band;

    IF_DEBUG(printf("[Pluto] pluto_detect_transformation_properties\n"););

    if (prog->nstmts == 0) return;

    HyperplaneProperties *hProps = prog->hProps;

    assert(prog->num_hyperplanes == stmts[0]->trans->nrows);

    // pluto_deps_print(stdout, prog);

    /* First compute satisfaction levels */
    pluto_compute_dep_directions(prog);
    pluto_compute_dep_satisfaction(prog);

    band = 0;
    level = 0;
    num_loops_in_band = 0;
    int bandStart = 0;

    do{
        for (i=0; i<prog->ndeps; i++)   {
            if (IS_RAR(deps[i]->type)) continue;
            if (deps[i]->satisfaction_level < level && 
                    hProps[deps[i]->satisfaction_level].type == H_SCALAR) continue;
            if (deps[i]->satisfaction_level >= bandStart 
                    && deps[i]->dirvec[level] != DEP_ZERO) 
                break;
        }

        if (i==prog->ndeps) {
            /* This band information is not used later; since band detection
             * is done again more accurately based on the scattering tree as
             * opposed to global depths; this is only used to output band
             * numbers conservatively when printing transformation properties */
            hProps[level].dep_prop = PARALLEL;
            hProps[level].band_num = band;
            if (hProps[level].type != H_SCALAR) num_loops_in_band++;
            level++;

        }else{

            for (i=0; i<prog->ndeps; i++)   {
                if (IS_RAR(deps[i]->type)) continue;
                if (deps[i]->satisfaction_level < level && 
                        hProps[deps[i]->satisfaction_level].type == H_SCALAR) continue;
                if (deps[i]->satisfaction_level >= bandStart 
                        && (deps[i]->dirvec[level] == DEP_MINUS 
                            || deps[i]->dirvec[level] == DEP_STAR))
                    break;
            }
            if (i==prog->ndeps) {
                hProps[level].dep_prop = PIPE_PARALLEL;
                hProps[level].band_num = band;
                if (hProps[level].type != H_SCALAR) num_loops_in_band++;

                level++;
            }else{
                /* Dependence violation if assertion fails: 
                 * basically, the current level has negative
                 * components for some unsatisfied dependence
                 */
                if (num_loops_in_band == 0) {
                    fprintf(stderr, "[Pluto] Unfortunately, the transformation computed has violated a dependence.\n");
                    fprintf(stderr, "\tPlease make sure there is no inconsistent/illegal .fst file in your working directory.\n");
                    fprintf(stderr, "\tIf not, this usually is a result of a bug in the dependence tester,\n");
                    fprintf(stderr, "\tor a bug in Pluto's auto transformation.\n");
                    fprintf(stderr, "\tPlease send this input file to the author if possible.\n");
                    IF_DEBUG(pluto_stmts_print(stdout, prog->stmts, prog->nstmts););
                    pluto_transformations_pretty_print(prog);
                    pluto_compute_dep_directions(prog);
                    pluto_print_dep_directions(prog);
                    assert(0);
                }

                band++;
                bandStart = level;
                if (num_loops_in_band == 1) {
                    if (hProps[level-1].dep_prop == PIPE_PARALLEL)
                        hProps[level-1].dep_prop = SEQ;
                }
                num_loops_in_band = 0;
            }

        }
    }while (level < prog->num_hyperplanes);

    if (num_loops_in_band == 1) {
        if (hProps[level-1].dep_prop == PIPE_PARALLEL)
            hProps[level-1].dep_prop = SEQ;
    }

    /* Permutable bands of loops could have inner parallel loops; they 
     * all have been detected as fwd_dep (except the outer parallel one of a band); 
     * we just modify those to parallel */
    for (i=0; i<prog->num_hyperplanes; i++)  {
        for (j=0; j<prog->ndeps; j++) {
            if (IS_RAR(deps[j]->type)) continue;
            if (deps[j]->satisfaction_level >= i && deps[j]->dirvec[i] != DEP_ZERO) break;
        }

        if (j==prog->ndeps)   {
            // couldn't have been marked sequential
            assert(hProps[i].dep_prop != SEQ);
            if (hProps[i].dep_prop == PIPE_PARALLEL)    {
                hProps[i].dep_prop = PARALLEL;
            }
        }
    }

    pluto_detect_hyperplane_types_stmtwise(prog);
}


void pluto_print_depsat_vectors(PlutoProg *prog, int levels)
{
    int i, j;
    Dep **deps;

    deps = prog->deps;

    printf("\nSatisfaction vectors for transformed program\n");

    for (i=0; i<prog->ndeps; i++) {
        assert(deps[i]->satvec != NULL);
        printf("Dep %d: S%d to S%d: ", i+1, deps[i]->src+1, deps[i]->dest+1);
        printf("(");
        for (j=0; j<levels; j++) {
            printf("%d, ", deps[i]->satvec[j]);
        }
        printf(")\n");
    }
}

void pluto_print_dep_directions(PlutoProg *prog)
{
    int i, j;

    Dep **deps = prog->deps;
    int ndeps = prog->ndeps;
    int nlevels = prog->num_hyperplanes;

    printf("\nDirection vectors for transformed program\n");

    for (i=0; i<ndeps; i++) {
        printf("Dep %d: S%d to S%d: ", i+1, deps[i]->src+1, deps[i]->dest+1);
        printf("(");
        for (j=0; j<nlevels; j++) {
            printf("%c, ", deps[i]->dirvec[j]);
        }
        printf(") satisfied: %s, satvec: (", 
                deps[i]->satisfied? "yes":"no");
        for (j=0; j<nlevels; j++) {
            printf("%d, ", deps[i]->satvec[j]);
        }
        printf(")\n");

        for (j=0; j<nlevels; j++) {
            if (deps[i]->dirvec[j] > 0)  {
                break;
            }
            if (deps[i]->dirvec[j] < 0) {
                printf("Dep %d violated: S%d to S%d\n", i, deps[i]->src+1, deps[i]->dest+1);
                printf("%d %d\n", deps[i]->satisfaction_level, deps[i]->satisfied);
            }else if (deps[i]->dirvec[j] < 0) {
                printf("Dep %d violated: S%d to S%d\n", i, deps[i]->src+1, deps[i]->dest+1);
                printf("%d %d\n", deps[i]->satisfaction_level, deps[i]->satisfied);
            }
        }

        printf("satvec: ");
        for (j=0; j<nlevels; j++) {
            printf("%d, ", deps[i]->satvec[j]);
        }
        printf("\n");
    }
}


/* 
 * 1. Pad statement domains to maximum domain depth to make it easier to construct
 * scheduling constraints. These will be removed before the actual ILP
 * runs, and just before pluto_auto_transform returns
 *
 * 2. Pre-process dependence domains to allow a single affine bounding
 * expression to be constructed. Fix this implementation (too brittle).
 */
void normalize_domains(PlutoProg *prog)
{
    int i, j, k;

    int nvar = prog->nvar;
    int npar = prog->npar;

    /* if a dep distance <= N and another <= M, how do you bound it, no
     * way to express max(N,M) as a single affine function, when space is
     * built for each dependence, it doesn't know anything about a parameter
     * that does not appear in its dpolyhedron, and so it will assign the
     * coeff corresponding to that param in the bounding function constraints
     * local to the dependence to zero, and what if some other dep needs that 
     * particular coeff to be >= 1 for bounding? 
     *
     * Solution: global context should be available
     *
     * How to construct? Just put together all constraints on parameters alone in
     * the global context, i.e., eliminate iterators out of each domain and
     * aggregate constraints on the parameters, and add them to each
     * dependence polyhedron
     */
    int count=0;
    if (npar >= 1)	{
        PlutoConstraints *context = pluto_constraints_alloc(prog->nstmts*npar, npar+1);
        pluto_constraints_set_names(context, prog->params);
        for (i=0; i<prog->nstmts; i++)    {
            PlutoConstraints *copy = pluto_constraints_dup(prog->stmts[i]->domain);
            for (j=0; j<prog->stmts[i]->dim_orig; j++)    {
                fourier_motzkin_eliminate(copy, 0);
            }
            assert(copy->ncols == npar+1);
            count += copy->nrows;

            if (count <= prog->nstmts*npar)    {
                pluto_constraints_intersect(context, copy);
                pluto_constraints_free(copy);
            }else{
                pluto_constraints_free(copy);
                break;
            }
        }
        pluto_constraints_simplify(context);
        if (options->debug) {
            printf("[Pluto] Global constraint context\n");
            pluto_constraints_compact_print(stdout, context);
        }

        /* Add context to every dep polyhedron */
        for (i=0; i<prog->ndeps; i++) {
            PlutoConstraints *dpolytope = prog->deps[i]->dpolytope;

            for (k=0; k<context->nrows; k++) {
                pluto_constraints_add_inequality(dpolytope);

                /* Already initialized to zero */

                for (j=0; j<npar+1; j++){
                    dpolytope->val[dpolytope->nrows-1][j+dpolytope->ncols-(npar+1)] = 
                        context->val[k][j];
                }
            }
            /* Update reference, add_row can resize */
            prog->deps[i]->dpolytope = dpolytope;
        }
        pluto_constraints_free(context);
    }else{
        IF_DEBUG(printf("No global context\n"));
    }


    /* Add padding dimensions to statement domains */
    for (i=0; i<prog->nstmts; i++)    {
        Stmt *stmt = prog->stmts[i];
        int orig_depth = stmt->dim_orig;
        assert(orig_depth == stmt->dim);
        for (j=orig_depth; j<nvar; j++)  {
            pluto_sink_statement(stmt, stmt->dim, 0, prog);
        }
    }


    for (i=0; i<prog->ndeps; i++)    {
        Dep *dep = prog->deps[i];
        int src_dim = prog->stmts[dep->src]->dim;
        int target_dim = prog->stmts[dep->dest]->dim;
        assert(dep->dpolytope->ncols == src_dim+target_dim+prog->npar+1);
    }

    /* Normalize rows of dependence polyhedra */
    for (k=0; k<prog->ndeps; k++)   {
        /* Normalize by gcd */
        PlutoConstraints *dpoly = prog->deps[k]->dpolytope;

        for(i=0; i<dpoly->nrows; i++)   {
            pluto_constraints_normalize_row(dpoly, i);
        }
    }

    /* Avoid the need for bounding function coefficients to take negative
     * values (TODO: should do this only for the bounding function constraints) */
    bool *neg = malloc(sizeof(bool)*npar);
    for (k=0; k<prog->ndeps; k++) {
        Dep *dep = prog->deps[k];
        PlutoConstraints *dpoly = dep->dpolytope;

        int j;
        bzero(neg, npar*sizeof(bool));

        for (j=2*nvar; j<2*nvar+npar; j++)  {
            int min = dpoly->val[0][j];
            int max = dpoly->val[0][j];
            for (i=1; i<dpoly->nrows; i++)  {
                min = PLMIN(dpoly->val[i][j], min);
                max = PLMAX(dpoly->val[i][j], max);
            }

            if (min < 0 && max <= 0)    {
                neg[j-2*nvar] = true;
                IF_DEBUG(printf("Dep %d has negative coeff's for parameter %d\n", 
                            dep->id, j-2*nvar+1));
            }
        }

        for (j=0; j<npar; j++)  {
            if (neg[j])   {
                pluto_constraints_add_inequality(dpoly);
                dpoly->val[dpoly->nrows-1][2*nvar+j] = 1;
            }
        }

    }
    free(neg);
}


/* Remove padding dimensions that were added earlier; transformation matrices
 * will have stmt->dim + npar + 1 after this function */
void denormalize_domains(PlutoProg *prog)
{
    int i, j;

    int nvar = prog->nvar;
    int npar = prog->npar;

    for (i=0; i<prog->nstmts; i++)  {
        int del_count;
        Stmt *stmt = prog->stmts[i];
        del_count = 0;
        for (j=0; j<nvar; j++)  {
            if (!stmt->is_orig_loop[j-del_count]) {
                pluto_stmt_remove_dim(stmt, j-del_count, prog);
                if (stmt->evicted_hyp) {
                    pluto_matrix_remove_col(stmt->evicted_hyp, j-del_count);
                }
                del_count++;
            }
        }

        assert(stmt->domain->ncols == stmt->dim+npar+1);
        assert(stmt->trans->ncols == stmt->dim+npar+1);

        for (j=0; j<stmt->dim; j++)  {
            stmt->is_orig_loop[j] = 1;
        }
    }
}

/* VB: Find the face that is allowing concurrent start
 * Used while lbtile option is set
 * Currently, the outermost loop is being assumed as the face with
 * concurrent start */
int *find_face_allowing_con_start(PlutoProg * prog)
{
    int i;
    int *face =  (int *) malloc(sizeof(int) * prog->nvar);
    for (i=0;i<prog->nvar;i++)
      face[i] = (i==0)?1:0;
    return face;
}


/*
 * Find hyperplane inside the cone  of previously found hyperplanes 
 * and the face allowing concurrent start
 *
 * evict_pos: position of the hyperplane to be evicted by the one that will
 * enable concurrent start
 *
 * cone_complement_pos: in case of partial concurrent start, the 
 * hyperplane that will form the cone with the conc start hyperplane
 *
 */
int find_cone_complement_hyperplane(int evict_pos, 
        PlutoMatrix **cone_complement_hyps, int cone_complement_pos, 
        PlutoConstraints *basecst, PlutoProg *prog)
{
    int i, j, k, lambda_k;
    int nstmts = prog->nstmts;
    int nvar = prog->nvar;
    int npar = prog->npar;
    Stmt **stmts = prog->stmts;

    IF_DEBUG(printf("[Pluto] finding cone complement hyperplane\n"););

    int64 *bestsol;
    PlutoConstraints *con_start_cst, *lastcst;

    /* lastcst is the set of additional constraints */
    lastcst = pluto_constraints_alloc(2*nvar*nstmts, CST_WIDTH+nvar*nstmts);

    /* all lambdas >=1 */
    for (i=0; i<nstmts; i++) {
        int stmt_offset = npar+1+nstmts*(nvar+1)+i*nvar;
        for (j=0; j<nvar; j++)  {
            pluto_constraints_add_inequality(lastcst);
            lastcst->val[lastcst->nrows-1][stmt_offset+j] =1;
            lastcst->val[lastcst->nrows-1][lastcst->ncols-1] = -1;
        }
    }

    /* Now, add the constraints for the new hyperplane to be in the cone
     * of the face and the negatives of the hyperplanes already found
     * (excluding the one being evicted: at `evict_pos') */
    for (i=0; i<nstmts; i++) {
        int stmt_offset1= npar+1+i*(nvar+1);
        int stmt_offset2= npar+1+nstmts*(nvar+1)+i*nvar;
        for (j=0; j<nvar; j++)  {
            pluto_constraints_add_equality(lastcst);
            lastcst->val[lastcst->nrows-1][stmt_offset1+j] =1;

            int *face = find_face_allowing_con_start(prog);
            lastcst->val[lastcst->nrows-1][stmt_offset2] = -(face[j]);
            free(face);

            if (options->partlbtile) {
                lastcst->val[lastcst->nrows-1][stmt_offset2+1] = 
                    prog->stmts[i]->trans->val[cone_complement_pos][j];
            }else{
                lambda_k=0;
                for(k=0; k<prog->stmts[i]->trans->nrows; k++){
                    if (k != evict_pos && prog->stmts[i]->hyp_types[k]!= H_SCALAR){
                        lastcst->val[lastcst->nrows-1][stmt_offset2+lambda_k+1] = prog->stmts[i]->trans->val[k][j];
                        lambda_k++;
                    }
                }
            }
            lastcst->val[lastcst->nrows-1][lastcst->ncols-1] = 0;
        }
    }

    /*
     * con_start_cst serves the same purpose as Pluto ILP formulation, but 
     * with expanded constraint-width to incorporate lambdas
     *
     * No need of non-zero solution constraints
     */
    con_start_cst = pluto_constraints_dup(basecst);

    PlutoConstraints *boundcst = get_coeff_bounding_constraints_for_cone_complement(prog);
    pluto_constraints_add(con_start_cst, boundcst);
    pluto_constraints_free(boundcst);

    for (i=0; i<nvar*nstmts; i++) {
        pluto_constraints_add_dim(con_start_cst, basecst->ncols-1, NULL);
    }

    pluto_constraints_add(con_start_cst, lastcst);
    pluto_constraints_free(lastcst);
    // printf("Cone complement constraints\n");
    // pluto_constraints_pretty_print(stdout, con_start_cst);

    /* pluto_constraints_lexmin is being called directly */
    bestsol = pluto_constraints_lexmin(con_start_cst, ALLOW_NEGATIVE_COEFF);
    pluto_constraints_free(con_start_cst);

    /* pluto_constraints_lexmin is being called directly */
    if (bestsol == NULL) {
        printf("[Pluto] No concurrent start possible\n");
    }else{
        IF_DEBUG(printf("[Pluto] Concurrent start possible\n"););
        for (j=0; j<nstmts; j++) {
            cone_complement_hyps[j] =
                pluto_matrix_alloc(1, stmts[j]->dim+npar+1);
        }
        for (j=0; j<nstmts; j++) {
            for (k=0; k<nvar; k++)    {
                cone_complement_hyps[j]->val[0][k] =
                    bestsol[npar+1+j*(nvar+1)+k];
            }
            /* No parameteric shifts */
            for (k=nvar; k<nvar+npar; k++)    {
                cone_complement_hyps[j]->val[0][k] = 0;
            }
            cone_complement_hyps[j]->val[0][nvar+npar] =
                bestsol[npar+1+j*(nvar+1)+nvar];

            IF_DEBUG(printf("S%d: cone complement\n", j+1););
            IF_DEBUG(pluto_matrix_print(stdout, cone_complement_hyps[j]););
        }
        free(bestsol);
    }

    return (cone_complement_hyps[0] == NULL)? 0:1; 
}

/* 
 * Check if the k'th row for any statement is the face 
 * allowing concurrent start
 */
int is_concurrent_start_face(PlutoProg *prog, int k)
{
    int i,j;
    for (i=0; i<prog->nstmts; i++){
        if (prog->stmts[i]->trans->val[k][0] != 1) return 0;
        for(j=1; j<prog->nvar; j++){
            if (prog->stmts[i]->trans->val[k][j] != 0) return 0;
        }
    }
    return 1;
}

/* 
 * Find the hyperplace parallel to the concurrent start face 
 * that will be evicted; if there is none, return
 * the last hyperplane
 */
int find_hyperplane_to_be_evicted(PlutoProg *prog, int first, int num_sols_found)
{
    int j;
    for (j=first; j<first+num_sols_found-1; j++){
        if (is_concurrent_start_face(prog, j)) return j;
    }
    /* Return the last one */
    return first + num_sols_found - 1; 
}



int get_first_non_scalar_hyperplane(PlutoProg *prog, int start, int end)
{
    int i;
    for (i=start; i<=end; i++){
        if (prog->hProps[i].type == H_LOOP) return i;
    }
    return -1;
}

/* Copy h2 into h1 */
static void copy_hyperplane(int64 *h1, int64 *h2, int ncols)
{
    int j;

    for(j=0; j<ncols; j++){
        h1[j] = h2[j];
    }
}

#if 0
/* Swap hyperplanes h1 and h2 */
static void swap_hyperplanes(int64 *h1, int64 *h2, int ncols)
{
    int64 tmp;
    int j;

    for(j=0; j<ncols; j++){
        tmp = h2[j];
        h2[j] = h1[j];
        h1[j] = tmp;
    }
}
#endif

void copy_comp_to_data_hyperplanes(Stmt *stmt, Array *arr,  int row){
    int i;

    //	assert(stmt->trans->ncols == arr->trans->ncols);

    for(i=0;i<arr->trans->ncols;i++)
        arr->trans->val[row][i] = stmt->trans->val[row][i];

    return;
}

void restore_data_hyperplane(Array *arr, int row){
    int i;

    if(options->identity_data_dist) return;

    for(i=0;i<arr->trans->ncols;i++)
        arr->trans->val[row][i] = arr->trans_orig->val[row][i];

    return;
}

PlutoConstraints* pluto_dist_get_nontrivial_solu_constraints(PlutoProg *prog){
    int stmt_cst_width = prog->nvar + prog->npar + 1;

    int arr_cst_width = prog->max_array_dim + prog->npar + 1;
    int cst_width = (prog->npar + 1) * (prog->narrays + 1) + prog->nstmts * (stmt_cst_width) + prog->narrays *(arr_cst_width) + 1;

    int cst_nrows = prog->narrays;
    PlutoConstraints *cst = pluto_constraints_alloc(cst_nrows, cst_width);
    cst->ncols = cst_width;
    cst->nrows = cst_nrows;

    int start_arr = (prog->npar + 1) * (prog->narrays + 1) + prog->nstmts*(stmt_cst_width);
    int i, j, offset, cur_row = 0;

    //sigma ai >= 1 to avoid trivial solu
    Array * arr;
    for(i=0;i<prog->narrays;i++){
        arr = prog->arrays[i];

        //		if(arr->dim_orig < 2) continue;

        if(arr->num_hyperplanes_found == arr->dim)
            continue;

        offset = arr->id * arr_cst_width + start_arr;

        cst->val[cur_row][cst_width-1] = -1;
        for(j=0;j<arr->dim;j++){
            cst->val[cur_row][offset+j] = 1;
        }

        cur_row++;
    }

    cst->nrows = cur_row;

    return cst;
}

void pluto_dist_pretty_print(int64 *sol, PlutoProg *prog){

    int i, j, k;
    int cur = 0;

    if(sol == NULL) {
        printf("No Solution found\n");
        return;
    }

    assert(sol != NULL);

    printf("\nComm volume due to dep = \n\t");
    for(i=0;i<prog->npar + 1; i++)
        printf("%lld\t", sol[i]);
    cur += prog->npar + 1;
    printf("\n");

    for(k=0;k<prog->narrays;k++){
        printf("\n Comm due to data hyperplanes  of Array %s = \n\t", prog->arrays[k]->text);
        for(i=0;i<prog->npar + 1; i++)
            printf("%lld\t", sol[cur + i]);
        cur += prog->npar + 1;
        printf("\n");
    }

    int num_cols_reverse = 0;

    //	printf("\n Comp Hyperplanes = \n");
    for(i=0;i<prog->nstmts;i++){
        printf("stmts_%d\t\t",i);
        num_cols_reverse = prog->stmts[i]->dim_orig;
        for(j=0;j<num_cols_reverse;j++){
            printf("%lld\t", sol[cur+ num_cols_reverse - 1 - j]);
        }
        cur += j;

        for(j=0;j<prog->npar+1;j++)
            printf("%lld\t",sol[cur + j]);

        cur += j;

        printf("\n");
    }

    //	printf("\n Data Hyperplanes = \n");
    for(i=0;i<prog->narrays;i++){
        printf("array_%s\t\t",prog->arrays[i]->text);
        num_cols_reverse = prog->arrays[i]->dim_orig;
        for(j=0;j< num_cols_reverse;j++){
            printf("%lld\t", sol[cur + num_cols_reverse -1 - j]);
        }
        cur += j;

        for(j=0;j<prog->npar+1;j++)
            printf("%lld\t",sol[cur + j]);

        cur += j;

        printf("\n");
    }

    return;
}

int64* pluto_dist_constraints_solve(PlutoConstraints *cst, PlutoProg *prog)
{
    int i, j, k;

    int64 *sol;

    print_polylib_visual_sets("cst", cst);

    //Remove the redundant cols and change order of cols
    int64 **mat = (int64 **)malloc(cst->nrows * sizeof(int64 *));
    for(i=0;i<cst->nrows;i++)
        mat[i] = (int64 *)malloc(cst->ncols * sizeof(int64));

    int curr_col = 0, cst_col = 0;

    //copy npar and const cols
    for(i=0;i< (prog->npar+1) * (prog->narrays + 1);i++)
        for(j=0;j<cst->nrows;j++)
            mat[j][i] = cst->val[j][i];

    curr_col += (prog->npar +1) * (prog->narrays + 1);
    cst_col +=  (prog->npar +1) * (prog->narrays + 1);

    int ncols_to_reverse;
    for(i=0;i<prog->nstmts;i++){
        ncols_to_reverse = prog->stmts[i]->dim_orig ;
        for(j=0;j< ncols_to_reverse;j++){
            for(k=0;k<cst->nrows; k++){
                mat[k][curr_col + j] = cst->val[k][cst_col + ncols_to_reverse - 1 -j];
            }
        }
        curr_col += ncols_to_reverse;
        cst_col += prog->nvar;

        //copy npar and const cols
        for(k=0;k<prog->npar+1;k++)
            for(j=0;j<cst->nrows;j++)
                mat[j][curr_col + k] = cst->val[j][cst_col+k];

        curr_col += prog->npar + 1;
        cst_col += prog->npar + 1;
    }

    for(i=0;i<prog->narrays;i++){
        ncols_to_reverse = prog->arrays[i]->dim_orig;
        for(j=0;j< ncols_to_reverse;j++){
            for(k=0;k<cst->nrows; k++){
                mat[k][curr_col + j] = cst->val[k][cst_col + ncols_to_reverse -1 -j];
            }
        }

        curr_col += ncols_to_reverse ;
        cst_col += prog->max_array_dim;

        //copy npar and const cols
        for(k=0;k<prog->npar+1;k++)
            for(j=0;j<cst->nrows;j++)
                mat[j][curr_col + k] = cst->val[j][cst_col+k];

        curr_col += prog->npar + 1;
        cst_col += prog->npar + 1;
    }

    for(i=0;i<cst->nrows;i++){
        mat[i][curr_col] = cst->val[i][cst_col];
    }

    cst_col++;
    curr_col++;

    assert(cst_col == cst->ncols);

    int64** prev = cst->val;
    int prev_ncols = cst->ncols;
    cst->val = mat;
    cst->ncols = curr_col;

    print_polylib_visual_sets("temp", cst);
    sol = pluto_constraints_lexmin(cst, DO_NOT_ALLOW_NEGATIVE_COEFF);

    //TODO: remap the format of sol

    cst->val = prev;
    cst->ncols = prev_ncols;

    return sol;
}

PlutoConstraints* pluto_dist_convert_constraints(PlutoProg *prog, PlutoConstraints *basecst){
    int i, j, k;
    int stmt_cst_width = prog->nvar + prog->npar + 1;
    int arr_cst_width = prog->max_array_dim + prog->npar + 1;
    int cst_width = (prog->npar + 1) * (prog->narrays + 1) + prog->nstmts * (stmt_cst_width) + prog->narrays *(arr_cst_width) + 1;

    // PlutoConstraints *basecst = get_permutability_constraints(prog->deps,prog->ndeps, prog);

    assert(basecst->ncols == (prog->npar + 1) + prog->nstmts * (prog->nvar + 1) + 1);

    PlutoConstraints *cst = pluto_constraints_alloc(basecst->nrows, cst_width);
    cst->nrows = basecst->nrows;
    cst->ncols = cst_width;

    int cst_curr_col = 0, base_curr_col = 0;

    //copy npar + 1 col
    for(i=0; i<prog->npar+1; i++){
        for(j=0;j<basecst->nrows;j++)
            cst->val[j][i] = basecst->val[j][i];
    }
    cst_curr_col += (prog->npar+1) * (prog->narrays + 1);
    base_curr_col += prog->npar+1;

    //for each statement copy comp hyperplane coeffs
    for(i=0;i<prog->nstmts;i++){
        for(j=0;j<prog->nvar;j++){
            for(k=0;k<basecst->nrows;k++){
                cst->val[k][cst_curr_col+j] = basecst->val[k][base_curr_col+j];
            }
        }

        cst_curr_col += prog->nvar + prog->npar;
        base_curr_col += prog->nvar;

        for(j=0;j<basecst->nrows;j++)
            cst->val[j][cst_curr_col] = basecst->val[j][base_curr_col];

        cst_curr_col++;
        base_curr_col++;
    }

    cst_curr_col += prog->narrays * arr_cst_width ;

    for(i=0;i<basecst->nrows;i++)
        cst->val[i][cst_curr_col] = basecst->val[i][base_curr_col];

    cst_curr_col++;
    base_curr_col++;

    assert(base_curr_col == basecst->ncols);
    assert(cst_curr_col == cst->ncols);

    return cst;
}


int pluto_dist_is_valid_hyperplane(int64 *sol, int start, int len){

    assert(sol != NULL) ;

    int i;

    for (i = 0; i < len; ++i) {
        if(sol[start + i] != 0)
            return 1;
    }

    return 0;
}


void pluto_dist_add_hyperplanes(int64 *sol, int curr_hyperplane,  PlutoProg *prog)
{
    int i, j;
    int cur = 0;

    pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_LOOP);

    if(sol == NULL)
        return;

    //Upper bound on comm vol
    //	for(i=0;i<prog->npar + 1; i++)
    //		printf("%d\t", sol[i]);
    cur += (prog->npar + 1) * (prog->narrays + 1) ;
    Stmt *stmt;

    //Add the comp hyperplanes
    for(i=0;i<prog->nstmts;i++){

        stmt = prog->stmts[i];
        pluto_stmt_add_hyperplane(stmt, H_UNKNOWN, stmt->trans->nrows);

        //		if(pluto_stmt_get_num_ind_hyps(prog->stmts[i]) == prog->stmts[i]->dim_orig){
        //			cur += stmt->dim_orig + prog->npar + 1;
        //			continue;
        //		}

        for(j=0;j<stmt->dim_orig;j++){
            stmt->trans->val[stmt->trans->nrows-1][j] = sol[cur+ stmt->dim_orig  - 1 - j];
        }
        cur += j;

        for(j=0;j<prog->npar+1;j++){
            stmt->trans->val[stmt->trans->nrows-1][prog->nvar + j] = sol[cur + j];
        }
        cur += j;

        stmt->hyp_types[stmt->trans->nrows-1] =
            pluto_is_hyperplane_scalar(stmt, stmt->trans->nrows-1)?
            H_SCALAR: H_LOOP;
    }

    Array *arr;
    //	printf("\n Data Hyperplanes = \n");
    for(i=0;i<prog->narrays;i++){

        arr = prog->arrays[i];

        if(!pluto_dist_is_valid_hyperplane(sol, cur, arr->dim)){
            cur += arr->dim_orig + prog->npar + 1;
            continue;
        }

        if(arr->num_hyperplanes_found == arr->dim_orig){
            cur += arr->dim_orig + prog->npar + 1;
            continue;
        }

        /* Store the mapping of data and comp hyperplane.
        */
        arr->hyperplane_mapping[curr_hyperplane] = arr->num_hyperplanes_found;

        pluto_matrix_add_row(arr->trans, arr->trans->nrows);

        for(j=0;j<arr->dim_orig;j++){
            arr->trans->val[arr->trans->nrows-1][j] = sol[cur + arr->dim_orig -1 - j];
        }
        cur += j;

        for(j=0;j<prog->npar+1;j++){
            arr->trans->val[arr->trans->nrows-1][prog->max_array_dim + j] = sol[cur+ j];
        }
        cur += j;

        arr->num_hyperplanes_found++;
    }

    return;

}

int is_access_scalar(PlutoAccess *access){

    int i = 0;

    if(access->mat->nrows > 1)
        return 0;

    for(i=0; i<access->mat->ncols; i++)
        if(access->mat->val[0][i] != 0)
            return 0;

    return 1;

}



void pluto_dist_verify_data_hyperplanes(PlutoProg *prog){

    int i = 0;
    Array *arr;
    for(i=0;i<prog->narrays;i++){

        arr = prog->arrays[i];

        if(prog->arrays[i]->dim_orig > prog->arrays[i]->num_hyperplanes_found){

            assert(0);


            pluto_matrix_add_row(arr->trans, arr->trans->nrows);

            if(arr->dim_orig == 1){
                arr->hyperplane_mapping[0] = arr->num_hyperplanes_found;
                arr->trans->val[arr->trans->nrows-1][0] = 1;
            }
            else{
                arr->hyperplane_mapping[1] = arr->num_hyperplanes_found;
                arr->trans->val[arr->trans->nrows-1][1] = 1;
            }

            arr->num_hyperplanes_found++;

        }

    }
}
/* Top level automatic transformation algorithm for distributed memory
 * Find all the linearly independent permutable band of data and
 * computation hyperplanes at a level
 */
int pluto_dist_find_permutable_hyperplanes(PlutoProg *prog, bool lin_ind_mode,
        bool loop_search_mode, int max_sols, int num_ind_sols_found)
{

    Stmt *stmt;
    int i=0, j=0;
    int stmt_cst_width = prog->nvar + prog->npar + 1;
    int arr_cst_width = prog->max_array_dim + prog->npar + 1;
    int cst_width = (prog->npar + 1) * (prog->narrays + 1) + prog->nstmts * (stmt_cst_width) + prog->narrays *(arr_cst_width) + 1;
    int cst_nrows = 0;
    int64 *sol;
    int num_sol_found = 0;
    //	max_sol = (prog->nvar>prog->max_array_dim)?prog->nvar:prog->max_array_dim;

    for(i=0; i<prog->nstmts; i++){
        stmt = prog->stmts[i];
        cst_nrows += stmt->nreads;
        cst_nrows += stmt->nwrites;
    }


    PlutoConstraints *cst = pluto_constraints_alloc(cst_nrows * 2 * stmt_cst_width + cst_width, cst_width);
    PlutoConstraints *temp, *pluto;

    for(i=0; i<prog->nstmts; i++){
        stmt = prog->stmts[i];

        //add constraints for read array accesses
        for(j=0; j<stmt->nreads; j++){

            if(is_access_scalar(stmt->reads[j]) || (stmt->reads[j]->mat->nrows == 1))
                continue;

            temp = pluto_dist_add_constraints_array_access(stmt->reads[j],stmt, prog);
            print_polylib_visual_sets( stmt->reads[j]->name, temp);
            pluto_constraints_add(cst, temp);
            pluto_constraints_free(temp);
            //			 sol = pluto_constraints_lexmin(cst,DO_NOT_ALLOW_NEGATIVE_COEFF);

        }
        //add constraints for write array accesses
        for(j=0; j<stmt->nwrites; j++){

            if(is_access_scalar(stmt->writes[j]) || (stmt->writes[j]->mat->nrows == 1))
                continue;

            temp = pluto_dist_add_constraints_array_access(stmt->writes[j], stmt, prog);
            print_polylib_visual_sets( stmt->reads[j]->name, temp);
            pluto_constraints_add(cst, temp);
            pluto_constraints_free(temp);
        }
    }


    PlutoConstraints *pluto_basecst = get_permutability_constraints(prog);
    print_polylib_visual_sets("pluto", pluto_basecst);

    pluto_dist_add_upper_lower_bounds(cst, prog);

    pluto = pluto_dist_convert_constraints(prog,pluto_basecst);

    pluto_constraints_add(cst, pluto);
    pluto_constraints_free(pluto);

    PlutoConstraints* currcst = pluto_constraints_alloc(cst->nrows, cst->ncols);
    PlutoConstraints *curr_pluto = pluto_constraints_alloc(pluto_basecst->nrows, pluto_basecst->ncols);

    PlutoConstraints *nzcst, *non_zero_cst;
    do{

        pluto_constraints_copy(currcst, cst);
        pluto_constraints_copy(curr_pluto, pluto_basecst);

        //Add constraints to avoid trivial solun for data hyperplanes
        non_zero_cst = pluto_dist_get_nontrivial_solu_constraints(prog);
        pluto_constraints_add(currcst, non_zero_cst);
        pluto_constraints_free(non_zero_cst);

        //Add constraints to avoid trivial solun for comp hyperplanes
        nzcst = get_non_trivial_sol_constraints(prog, EAGER);
        print_polylib_visual_sets("nzcst", nzcst);
        pluto_constraints_add(curr_pluto, nzcst);

        //Add linear independent for data hyperplanes
        PlutoConstraints *ind_cst_data = pluto_dist_get_linear_ind_constraints(prog, currcst , EAGER);
        print_polylib_visual_sets("ind", ind_cst_data);
        pluto_constraints_add(currcst, ind_cst_data);
        pluto_constraints_free(ind_cst_data);

        PlutoConstraints *indcst = NULL;
        print_polylib_visual_sets("curr", curr_pluto);
        PlutoConstraints *ind_cst = get_linear_ind_constraints(prog, curr_pluto, EAGER);
        print_polylib_visual_sets("ind", ind_cst);
        pluto_constraints_add(curr_pluto, ind_cst);

        pluto_constraints_add(ind_cst, nzcst);
        indcst = pluto_dist_convert_constraints(prog, ind_cst);
        print_polylib_visual_sets("ind", indcst);
        pluto_constraints_free(ind_cst);
        pluto_constraints_free(nzcst);
        pluto_constraints_add(currcst, indcst);

        if (indcst != NULL && indcst->nrows == 0) {
            /* If you don't have any independence constraints, we would end up finding
             * the same solution that was found earlier; so we won't find anything
             * new */
            sol = NULL;
        }else{

            sol = pluto_dist_constraints_solve(currcst, prog);
            pluto_dist_pretty_print(sol, prog);

            if(sol != NULL){
                pluto_dist_add_hyperplanes(sol,num_ind_sols_found+num_sol_found,prog);

                free(sol);
                num_sol_found++;
            }

        }
        pluto_constraints_free(indcst);

    }while(num_sol_found < max_sols && sol != NULL);

    //add the dep cosntraints on comp hyperplanes

    pluto_constraints_free(pluto_basecst);
    pluto_constraints_free(currcst);
    return  num_sol_found;
}

void pluto_dist_add_upper_lower_bounds(PlutoConstraints *cst, PlutoProg *prog){
    int i=0;

    //int stmt_cst_width = prog->nvar + prog->npar + 1;
    //int arr_cst_width = prog->max_array_dim + prog->npar + 1;
    //int cst_width = (prog->npar + 1) * (prog->narrays + 1) + prog->nstmts * (stmt_cst_width) + prog->narrays *(arr_cst_width) + 1;
    int num_bounding_coeff = (prog->npar + 1) *( 1 + prog->narrays);

    /* Add upper bounds for transformation coefficients */
    //TODO: Add function to get upper bound for data hyperplanes.
    int ub = pluto_prog_get_largest_const_in_domains(prog);

    //Check the upper bound
    if(ub < 10) ub = 25;

    /* Putting too small an upper bound can prevent useful transformations;
     * also, note that an upper bound is added for all statements globally due
     * to the lack of an easy way to determine bounds for each coefficient to
     * prevent spurious transformations that involve shifts proportional to
     * loop bounds
     */
    if (ub >= 10)   {
        for (i=0; i<cst->ncols - num_bounding_coeff -1; i++)  {
            IF_DEBUG2(printf("Adding upper bound %d for transformation coefficients\n", ub););
            pluto_constraints_add_ub(cst, num_bounding_coeff+i, ub);
        }
    }
    for (i=0; i< cst->ncols- num_bounding_coeff -1 ; i++)  {
        pluto_constraints_add_lb(cst, num_bounding_coeff+i, -5);
    }

    return;
}


void pluto_data_dist_identity_trans(PlutoProg *prog){

    int i = 0, j = 0;
    int max_dim = 0;

    for(i=0;i<prog->nstmts;i++){

        max_dim = prog->stmts[i]->dim > max_dim ? prog->stmts[i]->dim : max_dim;
    }

    for(i=0;i<prog->narrays;i++){
        Array *arr = prog->arrays[i];

        arr->trans = pluto_matrix_identity(arr->dim_orig);

        for(j=0;j<prog->npar+1;j++)
            pluto_matrix_add_col(arr->trans, arr->trans->ncols);

        arr->num_hyperplanes_found = arr->dim_orig;

        if(arr->hyperplane_mapping == NULL)
            arr->hyperplane_mapping = malloc(max_dim* sizeof(int));

        for(j=0;j<max_dim;j++){
            arr->hyperplane_mapping[j] = -1;
        }

        for(j=0;j<max_dim;j++){
            arr->hyperplane_mapping[j] = j%arr->dim_orig;
        }
    }

}
/* 
 * Top-level automatic transformation algoritm 
 *
 * All dependences are reset to unsatisfied before starting
 *
 */ 
int pluto_auto_transform(PlutoProg *prog)
{
    int nsols, i, j, conc_start_found, depth;
    /* The maximum number of independent solutions needed across all stmts */
    int num_ind_sols;
    bool lin_ind_mode;
    bool loop_search_mode;

    int first_loop_hyp;
    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i]->satisfied = false;
    }

    /* Create the data dependence graph */
    prog->ddg = ddg_create(prog);
    ddg_compute_scc(prog);

    Graph *ddg = prog->ddg;
    int nvar = prog->nvar;
    int npar = prog->npar;

    if (nstmts == 0)  return 0;

    normalize_domains(prog);

    PlutoMatrix **orig_trans = malloc(nstmts*sizeof(PlutoMatrix *));
    int orig_num_hyperplanes = prog->num_hyperplanes;
    HyperplaneProperties *orig_hProps = prog->hProps;

    lin_ind_mode = EAGER;
    loop_search_mode = EAGER;

    /* Get rid of any existing transformation */
    for (i=0; i<nstmts; i++) {
        Stmt *stmt = prog->stmts[i];
        /* Save the original transformation */
        orig_trans[i] = stmt->trans;
        /* Pre-allocate a little more to prevent frequent realloc */
        stmt->trans = pluto_matrix_alloc(2*stmt->dim+1, stmt->dim+npar+1);
        stmt->trans->nrows = 0;
    }

    if (options->data_dist) {
        for (i=0; i<prog->narrays; i++) {
            Array *arr = prog->arrays[i];
            arr->trans = pluto_matrix_alloc(2*arr->dim, arr->dim+prog->npar+1);
            arr->trans->nrows = 0;
            arr->num_hyperplanes_found = 0;
        }
    }

    prog->num_hyperplanes = 0;
    prog->hProps = NULL;

    /* The number of independent solutions required for the deepest 
     * statement */
    nsols = 0;
    for (i=0; i<nstmts; i++)    {
        nsols = PLMAX(nsols, stmts[i]->dim);
    }

    depth = 0;

    if (precut(prog, ddg, depth))   {
        /* Distributed based on .fst or .precut file (customized user-supplied
         * fusion structure */
        num_ind_sols = pluto_get_max_ind_hyps(prog);
        printf("[Pluto] Forced custom fusion structure from .fst/.precut\n");
        IF_DEBUG(fprintf(stdout, "%d ind solns in .precut file\n", 
                    num_ind_sols));
    }else{
        num_ind_sols = 0;
        if (options->fuse == SMART_FUSE)    {
            cut_scc_dim_based(prog,ddg);
        }
    }

    if (options->data_dist){
        for (i=0; i<prog->narrays; i++) {
            Array *arr = prog->arrays[i];
            if(arr->hyperplane_mapping == NULL)
                arr->hyperplane_mapping = malloc(nsols* sizeof(int));

            for (j=0; j<nsols; j++) {
                arr->hyperplane_mapping[j] = -1;
            }
        }
    }

    /* Diamond tiling */
    conc_start_found = 0;
    
    do{
        int num_sols_found, num_sols_left, s;

        if (options->fuse == NO_FUSE)   {
            ddg_compute_scc(prog);
            cut_all_sccs(prog, ddg);
        }

        /*
         * nsols - num_ind_sols is not the number of remaining hyperplanes
         * to be found in the LAZY mode (it is for the EAGER
         * mode). In LAZY mode, there may be more to be found for *some*
         * statements
         */
        num_sols_left = 0;
        for (s=0; s<nstmts; s++) {
            num_sols_left = PLMAX(num_sols_left, stmts[s]->dim_orig
                    - pluto_stmt_get_num_ind_hyps(stmts[s]));
        }
        assert(lin_ind_mode == LAZY || num_sols_left == nsols - num_ind_sols);

        if (options->data_dist && !options->identity_data_dist) {
            num_sols_found = pluto_dist_find_permutable_hyperplanes(prog, 
                    lin_ind_mode, loop_search_mode,
                    num_sols_left, num_ind_sols);
        }else{
            num_sols_found = find_permutable_hyperplanes(prog, lin_ind_mode,
                    loop_search_mode, num_sols_left);
        }

        IF_DEBUG(fprintf(stdout, "Level: %d; \t%d hyperplanes found\n",
                    depth, num_sols_found));
        IF_DEBUG2(pluto_transformations_pretty_print(prog));

        num_ind_sols = pluto_get_max_ind_hyps(prog);
        first_loop_hyp = 0;

        /* Diamond tiling: done for the first band of permutable loops */
        if (options->lbtile && num_ind_sols >= 2 && !conc_start_found) {
            conc_start_found = pluto_diamond_tile(prog);
        }
        
        if (num_sols_found >= 1) {
            for (j=0; j<num_sols_found; j++)      {
                /* Mark dependences satisfied by this solution */
                dep_satisfaction_update(prog,
                        stmts[0]->trans->nrows-num_sols_found+j);
                ddg_update(ddg, prog);
            }
        }else{
            /* Satisfy inter-scc dependences via distribution since we have 
             * no more fusable loops */

            ddg_compute_scc(prog);

            if (ddg->num_sccs >= 2) {
                if (options->fuse == NO_FUSE)  {
                    /* No fuse */
                    cut_all_sccs(prog, ddg);
                }else if (options->fuse == SMART_FUSE)  {
                    /* Smart fuse (default) */
                    cut_smart(prog, ddg);
                }else{
                    /* Max fuse */
                    if (depth >= 2*nvar+1) cut_all_sccs(prog, ddg);
                    else cut_conservative(prog, ddg);
                }
            }else{
                /* Only one SCC */
                if (lin_ind_mode == EAGER)   {
                    IF_DEBUG(printf("Switching to LAZY mode\n"););
                    lin_ind_mode = LAZY;
                    /* loop_search_mode = LAZY; */
                }else{
                    /* LAZY mode */
                    assert(lin_ind_mode == LAZY);
                    /* There is a problem; solutions should have been found */
                    if (options->debug) {
                        printf("Number of unsatisfied deps: %d\n", 
                                get_num_unsatisfied_deps(prog->deps, prog->ndeps));
                        printf("Number of unsatisfied inter-stmt deps: %d\n", 
                                get_num_unsatisfied_inter_stmt_deps(prog->deps, prog->ndeps));
                        IF_DEBUG(pluto_stmts_print(stdout, prog->stmts, prog->nstmts););
                        fprintf(stderr, "[Pluto] Unfortunately, pluto cannot find any more hyperplanes.\n");
                        fprintf(stderr, "\tThis is usually a result of (1) a bug in the dependence tester,\n");
                        fprintf(stderr, "\tor (2) a bug in Pluto's auto transformation,\n");
                        fprintf(stderr, "\tor (3) an inconsistent .fst/.precut in your working directory.\n");
                        fprintf(stderr, "\tor (4) or a case where the PLUTO algorithm doesn't succeed\n");
                        pluto_transformations_pretty_print(prog);
                        pluto_compute_dep_directions(prog);
                        pluto_print_dep_directions(prog);
                    }
                    denormalize_domains(prog);
                    printf("[Pluto] WARNING: working with original (identity) transformation (if they exist)\n");
                    /* Restore original ones */
                    for (i=0; i<nstmts; i++) {
                        stmts[i]->trans = orig_trans[i];
                        prog->num_hyperplanes = orig_num_hyperplanes;
                        prog->hProps = orig_hProps;
                    }
                    return 1;
                }
            }
        }
        depth++;
    }while (!pluto_transformations_full_ranked(prog) || 
            !deps_satisfaction_check(prog->deps, prog->ndeps));


    if (options->lbtile && !conc_start_found) {
        PLUTO_MESSAGE(printf("[Pluto] Diamond tiling not possible/useful\n"););
    }

    if (options->data_dist) {
        copy_comp_to_data_hyperplanes(prog->stmts[0], prog->arrays[0], first_loop_hyp);
    }

    denormalize_domains(prog);

    if(options->data_dist && !options->identity_data_dist)
        pluto_dist_verify_data_hyperplanes(prog);

    for (i=0; i<nstmts; i++)    {
        pluto_matrix_free(orig_trans[i]);
    }
    free(orig_trans);
    if (orig_hProps) free(orig_hProps);

    return 0;
}


int get_num_unsatisfied_deps(Dep **deps, int ndeps)
{
    int i, count;

    count = 0;
    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i]->type))   continue;
        if (!deps[i]->satisfied)  {
            IF_DEBUG(printf("Unsatisfied dep %d\n", i+1));
            count++;
        }
    }

    return count;

}


int get_num_unsatisfied_inter_stmt_deps(Dep **deps, int ndeps)
{
    int i;

    int count = 0;
    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i]->type))   continue;
        if (deps[i]->src == deps[i]->dest)    continue;
        if (!deps[i]->satisfied)  {
            IF_DEBUG(printf("Unsatisfied dep %d\n", i+1));
            count++;
        }
    }

    return count;

}


void ddg_print(Graph *g)
{
    pluto_matrix_print(stdout, g->adj);
}


/* Update the DDG - should be called when some dependences 
 * are satisfied */
void ddg_update(Graph *g, PlutoProg *prog)
{
    int i, j;
    Dep *dep;

    for (i=0; i<g->nVertices; i++) 
        for (j=0; j<g->nVertices; j++)
            g->adj->val[i][j] = 0;

    for (i=0; i<prog->ndeps; i++)   {
        dep = prog->deps[i];
        if (IS_RAR(dep->type)) continue;
        /* Number of unsatisfied dependences b/w src and dest is stored in the
         * adjacency matrix */
        g->adj->val[dep->src][dep->dest] += !dep_is_satisfied(dep);
    }
}


/* 
 * Create the DDG (RAR deps not included) from the unsatisfied deps
 */
Graph *ddg_create(PlutoProg *prog)
{
    int i;

    Graph *g = graph_alloc(prog->nstmts);

    for (i=0; i<prog->ndeps; i++)   {
        Dep *dep = prog->deps[i];
        /* no input dep edges in the graph */
        if (IS_RAR(dep->type)) continue;
        /* remember it's a multi-graph */
        g->adj->val[dep->src][dep->dest] += !dep_is_satisfied(dep);
    }

    return g;
}


/* 
 * Get the dimensionality of the stmt with max dimensionality in the SCC
 */
static int get_max_orig_dim_in_scc(PlutoProg *prog, int scc_id)
{
    int i;

    int max = -1;
    for (i=0; i<prog->nstmts; i++)  {
        Stmt *stmt = prog->stmts[i];
        if (stmt->scc_id == scc_id) {
            max = PLMAX(max,stmt->dim_orig);
        }
    }

    return max;
}

/* Number of vertices in a given SCC */
static int get_scc_size(PlutoProg *prog, int scc_id)
{
    int i;
    Stmt *stmt;

    int num = 0;
    for (i=0; i<prog->nstmts; i++)  {
        stmt = prog->stmts[i];
        if (stmt->scc_id == scc_id) {
            num++;
        }
    }

    return num;
}


/* Compute the SCCs of a graph */
void ddg_compute_scc(PlutoProg *prog)
{
    int i;

    Graph *g = prog->ddg;

    dfs(g);

    Graph *gT = graph_transpose(g);

    dfs_for_scc(gT);

    g->num_sccs = gT->num_sccs;

    for (i=0; i<g->nVertices; i++)  {
        g->vertices[i].scc_id = gT->vertices[i].scc_id;
        int stmt_id = gT->vertices[i].id;
        assert(stmt_id == i);
        prog->stmts[i]->scc_id = g->vertices[i].scc_id;
    }

    for (i=0; i<g->num_sccs; i++)  {
        g->sccs[i].max_dim = get_max_orig_dim_in_scc(prog, i);
        g->sccs[i].size = get_scc_size (prog, i);
        g->sccs[i].id = gT->sccs[i].id;
    }

    graph_free(gT);

    graph_print_sccs(g);
}

/* Get this arrays's schedule
 * Schedule format
 * [num sched functions | orig dim iters | params | const ]
 * Number of rows == num sched functions (each row for one hyperplane)
 */
PlutoConstraints *pluto_arr_get_schedule(const Array *arr)
{
    int i;

    PlutoMatrix *sched, *trans;
    PlutoConstraints *schedcst;

    trans = arr->trans_orig;
    sched = pluto_matrix_dup(trans);

    for (i=0; i<sched->nrows; i++)  {
        pluto_matrix_negate_row(sched, sched->nrows-1-i);
        pluto_matrix_add_col(sched, 0);
        sched->val[trans->nrows-1-i][0] = 1;
    }

    schedcst = pluto_constraints_from_equalities(sched);

    pluto_matrix_free(sched);

    return schedcst;
}


/* Get this statement's schedule
 * Schedule format
 * [num sched functions | orig dim iters | params | const ]
 * Number of rows == num sched functions (each row for one hyperplane)
 */
PlutoConstraints *pluto_stmt_get_schedule(const Stmt *stmt)
{
    int i;

    PlutoMatrix *sched, *trans;
    PlutoConstraints *schedcst;

    trans = stmt->trans;
    sched = pluto_matrix_dup(trans);

    for (i=0; i<sched->nrows; i++)  {
        pluto_matrix_negate_row(sched, sched->nrows-1-i);
        pluto_matrix_add_col(sched, 0);
        sched->val[trans->nrows-1-i][0] = 1;
    }

    schedcst = pluto_constraints_from_equalities(sched);

    pluto_matrix_free(sched);

    return schedcst;
}


PlutoConstraints *pluto_get_transformed_dpoly(const Dep *dep, Stmt *src, Stmt *dest)
{
    int i, npar;
    PlutoConstraints *src_sched, *dest_sched;

    npar = src->domain->ncols - src->dim - 1;

    //pluto_constraints_print(stdout, dep->dpolytope);
    //printf("%d %d\n", src->dim, dest->dim);
    assert(dep->dpolytope->ncols == src->dim + dest->dim + npar +1);

    PlutoConstraints *dpoly = pluto_constraints_dup(dep->dpolytope);

    // IF_DEBUG(printf("Original dpoly is \n"););
    // IF_DEBUG(pluto_constraints_print(stdout, dpoly););

    for (i=0; i<src->trans->nrows; i++)  {
        pluto_constraints_add_dim(dpoly, 0, NULL);
    }
    for (i=0; i<dest->trans->nrows; i++)  {
        pluto_constraints_add_dim(dpoly, src->trans->nrows+src->dim, NULL);
    }

    src_sched = pluto_stmt_get_schedule(src);
    dest_sched = pluto_stmt_get_schedule(dest);

    for (i=0; i<dest->trans->nrows+dest->dim; i++) {
        pluto_constraints_add_dim(src_sched, src->trans->nrows+src->dim, NULL);
    }

    for (i=0; i<src->trans->nrows+src->dim; i++) {
        pluto_constraints_add_dim(dest_sched, 0, NULL);
    }

    pluto_constraints_add(dpoly, src_sched);
    pluto_constraints_add(dpoly, dest_sched);

    // IF_DEBUG(printf("New pre-domain is \n"););
    // IF_DEBUG(pluto_constraints_print(stdout, newdom););

    pluto_constraints_project_out(dpoly, 
            src->trans->nrows, src->dim);

    pluto_constraints_project_out(dpoly, 
            src->trans->nrows+dest->trans->nrows, dest->dim);

    // IF_DEBUG(printf("New domain is \n"););
    // IF_DEBUG(pluto_constraints_print(stdout, newdom););

    pluto_constraints_free(src_sched);
    pluto_constraints_free(dest_sched);

    return dpoly;
}


/* Compute region(s) of data accessed by 'acc' with 'copy_level' number of
 * outer loops as parameters
 * 'copy_level' outer dimensions will be treated as parameters in addition 
 * to global ones 
 * domain: set (iterations of stmt) accessing data - in transformed space
 * acc: original access function
 * Input format: [copy_level, stmt->trans->nrows, prog->npar, 1] 
 *                or [copy_level, stmt->trans->nrows-copy_level, prog->npar, 1]
 * 
 * Output format:  [copy_level, acc->nrows, prog->npar + 1]
 * */
PlutoConstraints *pluto_compute_region_data(const Stmt *stmt, 
        const PlutoConstraints *domain, const PlutoAccess *acc, 
        int copy_level, const PlutoProg *prog)
{
    int i, k, npar, *divs;

    assert(acc->mat != NULL);
    assert(copy_level >= 0 && copy_level <= stmt->trans->nrows);

    npar = prog->npar;

    assert((stmt->trans->nrows+npar+1 == domain->ncols) 
            || (copy_level+stmt->trans->nrows+npar+1 == domain->ncols));

    PlutoMatrix *newacc = pluto_get_new_access_func(stmt, acc->mat, &divs);

    PlutoConstraints *datadom = pluto_constraints_dup(domain);

    assert(newacc->ncols == stmt->trans->nrows + npar + 1);

    for (k=0; k<newacc->nrows; k++) {
        pluto_matrix_negate_row(newacc, newacc->nrows-1-k);
        pluto_matrix_add_col(newacc, stmt->trans->nrows);
        newacc->val[newacc->nrows-1-k][stmt->trans->nrows] = divs[k];

        pluto_constraints_add_dim(datadom, domain->ncols-prog->npar-1, NULL);
    }

    PlutoConstraints *acc_cst = pluto_constraints_from_equalities(newacc);

    for (i=0; i<domain->ncols-stmt->trans->nrows-npar-1; i++) {
        pluto_constraints_add_dim(acc_cst, 0, NULL);
    }

    pluto_constraints_add_to_each(datadom, acc_cst);

    pluto_constraints_project_out(datadom, copy_level,
            datadom->ncols-copy_level-npar-1-newacc->nrows);

    //IF_DEBUG(printf("compute_region_data: data set written to\n"););
    //IF_DEBUG(pluto_constraints_print(stdout, datadom););

    pluto_constraints_free(acc_cst);
    pluto_matrix_free(newacc);

    if (domain->next != NULL)   {
        datadom->next = 
            pluto_compute_region_data(stmt, domain->next, acc, copy_level, prog);
    }

    return datadom;
}


/* Update a dependence with a new constraint added to the statement domain */
void pluto_update_deps(Stmt *stmt, PlutoConstraints *cst, PlutoProg *prog)
{
    int i, c;

    Stmt **stmts = prog->stmts;

    assert(cst->ncols == stmt->domain->ncols);

    for (i=0; i<prog->ndeps; i++) {
        Dep *dep = prog->deps[i];
        if (stmts[dep->src] == stmt) {
            PlutoConstraints *cst_l = pluto_constraints_dup(cst);
            Stmt *tstmt = stmts[dep->dest];
            for (c=0; c<tstmt->dim; c++) {
                pluto_constraints_add_dim(cst_l, stmt->dim, NULL);
            }
            pluto_constraints_add(dep->dpolytope, cst_l);
            pluto_constraints_free(cst_l);
        }
        if (stmts[dep->dest] == stmt) {
            PlutoConstraints *cst_l = pluto_constraints_dup(cst);
            Stmt *sstmt = stmts[dep->src];
            for (c=0; c<sstmt->dim; c++) {
                pluto_constraints_add_dim(cst_l, 0, NULL);
            }
            pluto_constraints_add(dep->dpolytope, cst_l);
            pluto_constraints_free(cst_l);
        }
    }

    for (i=0; i<prog->ntransdeps; i++) {
        Dep *dep = prog->transdeps[i];
        if (stmts[dep->src] == stmt) {
            PlutoConstraints *cst_l = pluto_constraints_dup(cst);
            Stmt *tstmt = stmts[dep->dest];
            for (c=0; c<tstmt->dim; c++) {
                pluto_constraints_add_dim(cst_l, stmt->dim, NULL);
            }
            pluto_constraints_add(dep->dpolytope, cst_l);
            pluto_constraints_free(cst_l);
        }
        if (stmts[dep->dest] == stmt) {
            PlutoConstraints *cst_l = pluto_constraints_dup(cst);
            Stmt *sstmt = stmts[dep->src];
            for (c=0; c<sstmt->dim; c++) {
                pluto_constraints_add_dim(cst_l, 0, NULL);
            }
            pluto_constraints_add(dep->dpolytope, cst_l);
            pluto_constraints_free(cst_l);
        }
    }
}

/*
 * Detect scattering functions that map to a single value, and modify the
 * scattering function to set it to that value; if the domain is  empty, this 
 * function ends up setting all of the scattering functions to zero */
void pluto_detect_scalar_dimensions(PlutoProg *prog)
{
    int s, d;

    for (s=0; s<prog->nstmts; s++) {
        Stmt *stmt = prog->stmts[s];
        PlutoConstraints *newdom = pluto_get_new_domain(stmt);

        /* If it's empty, we'll just change it to zero */
        int is_empty = pluto_constraints_is_empty(newdom);
        if (is_empty) {
            for (d=0; d<stmt->trans->nrows; d++) {
                pluto_matrix_zero_row(stmt->trans, d);
            }
            continue;
        }

        for (d=0; d<stmt->trans->nrows; d++) {
            int is_const_ub, is_const_lb;
            int64 ub, lb;
            is_const_lb = pluto_constraints_get_const_lb(newdom, d, &lb);
            is_const_ub = pluto_constraints_get_const_ub(newdom, d, &ub);
            if (is_const_lb && is_const_ub) {
                if (ub == lb) {
                    pluto_matrix_zero_row(stmt->trans, d);
                    stmt->trans->val[d][stmt->trans->ncols-1] = -lb;
                }
            }
        }
        pluto_constraints_free(newdom);
    }
}


/* Are these statements completely fused until the innermost level */
int pluto_are_stmts_fused(Stmt **stmts, int nstmts, const PlutoProg *prog)
{
    int num;

    if (prog->num_hyperplanes <= 1) return 1;

    Ploop **loops = pluto_get_loops_under(stmts, nstmts, prog->num_hyperplanes-2, prog, &num);
    // pluto_loops_print(loops, num);
    pluto_loops_free(loops, num);

    return num==1;
}


/* 
 * For diamond tiling
 */
int pluto_diamond_tile(PlutoProg *prog)
{
    int b, nbands, conc_start_found;

    conc_start_found = 0;

    /* Get the permutability constraints since a call to
     * detect_transformation_properties with update dep satisfaction levels
     * and we won't get the constraints we want */

    /* Don't free basecst */
    PlutoConstraints *basecst = get_permutability_constraints(prog);

    /* 
     * For partial concurrent, if we are trying to find a hyperplane c 
     * such that the face f is in the cone of c and a previously found 
     * hyperplane h, then h is the cone_complement of c
     */
    pluto_detect_transformation_properties(prog);
    Band **bands = pluto_get_outermost_permutable_bands(prog, &nbands);

    for (b=0; b<nbands; b++) {
        PlutoMatrix **cone_complement_hyps;
        Band *band = bands[b];
        int conc_start_found_band, evict_pos;
        int i, first_loop_hyp, cone_complement_pos;

        /* Band should not have outer parallelism */
        if (pluto_loop_is_parallel(prog, band->loop)) continue;

        /* Band should have inner parallelism */
        int ni, s;
        Ploop **iloops = pluto_get_loops_immediately_inner(band->loop, prog, &ni);
        for (i=0; i<ni; i++) {
            for (s=0; s<band->loop->nstmts; s++) {
                if (!pluto_loop_is_parallel_for_stmt(prog, iloops[i], 
                            band->loop->stmts[s])) break;
            }
            if (s<band->loop->nstmts) break;
        }
        if (i<ni) {
            pluto_loops_free(iloops, ni);
            continue;
        }

        /* Pure Inner parallelism should be lost via tiling */
        for (i=0; i<ni; i++) {
            if (!pluto_loop_has_satisfied_dep_with_component(
                        prog, iloops[i])) break;
        }
        pluto_loops_free(iloops, ni);
        if (i<ni) continue;

        cone_complement_hyps = malloc(
                prog->nstmts*sizeof(PlutoMatrix *));
        for (i=0; i<prog->nstmts; i++) {
            cone_complement_hyps[i] = NULL;
        }
        for (i=0; i<band->loop->nstmts; i++) {
            cone_complement_hyps[i] = NULL;
        }

        first_loop_hyp = band->loop->depth;
        /* 
         * Find hyperplane that will be replaced by the newly found
         * hyperplane
         * Concurrent start pertains to the first band alone
         */
        evict_pos = find_hyperplane_to_be_evicted(prog, 
                first_loop_hyp, band->width);
        /* If we haven't yet found the cone_complement_pos, just 
         * choose the first one as the cone_complement_pos */
        cone_complement_pos = first_loop_hyp;

        /* If first_loop_hyp hyperplane itself is to be replaced, 
         * choose the next one as cone_complement_pos */
        if (evict_pos == first_loop_hyp) cone_complement_pos++ ;
        conc_start_found_band = find_cone_complement_hyperplane(evict_pos,
                cone_complement_hyps, cone_complement_pos, basecst, prog);

        /* Re-arrange the transformation matrix if concurrent start
         * was found, store the replaced hyperplane so that it can be 
         * put back for the right intra-tile order */
        if (conc_start_found_band) {
            IF_DEBUG(printf("[Pluto] Transformations before concurrent start enable\n")); 
            IF_DEBUG(pluto_transformations_pretty_print(prog););
            for (i=0; i<band->loop->nstmts; i++){
                Stmt *stmt = band->loop->stmts[i];
                /* Since we do concurrent start only once */
                assert(stmt->evicted_hyp == NULL);
                stmt->evicted_hyp = pluto_matrix_alloc(1, stmt->trans->ncols);
                copy_hyperplane(stmt->evicted_hyp->val[0], 
                        stmt->trans->val[evict_pos], stmt->trans->ncols);
                copy_hyperplane(stmt->trans->val[evict_pos], 
                        cone_complement_hyps[stmt->id]->val[0], stmt->trans->ncols);
            }
            prog->evicted_hyp_pos = evict_pos;
            PLUTO_MESSAGE(printf("[Pluto] Concurrent start hyperplanes found\n"););
            IF_DEBUG(printf("[Pluto] Transformations after concurrent start enable\n")); 
            IF_DEBUG(pluto_transformations_pretty_print(prog););
        }

        conc_start_found |= conc_start_found_band;

        for (i=0; i<prog->nstmts; i++) {
            pluto_matrix_free(cone_complement_hyps[i]);
        }
        free(cone_complement_hyps);
    }

    pluto_bands_free(bands, nbands);

    return conc_start_found;
}
