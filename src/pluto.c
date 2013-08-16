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

/*
 * Returns the number of (new) satisfied dependences at this level
 */
int dep_satisfaction_update(PlutoProg *prog, int level)
{
    int i;
    int num_new_carried;

    int ndeps = prog->ndeps;
    Dep **deps = prog->deps;

    num_new_carried=0;

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
        prog->deps[i]->satisfaction_level =  prog->num_hyperplanes;
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
                dep->satvec[level] =  1;
            }else{
                dep->satvec[level] =  0;
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


/*
 * This calls pluto_constraints_solve, but before doing that does some preprocessing
 * - removes variables that we know will be assigned 0 - also do some
 *   permutation of the variables to get row-wise access
 */
int *pluto_prog_constraints_solve(PlutoConstraints *cst, PlutoProg *prog)
{
    Stmt **stmts;
    int nstmts, nvar, npar;

    stmts  = prog->stmts;
    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;

    /* Remove redundant variables - that don't appear in your outer loops */
    int redun[npar+1+nstmts*(nvar+1)+1];
    int i, j, k, q;
    int *sol, *fsol;
    PlutoConstraints *newcst;

    assert(cst->ncols-1 == npar+1+nstmts*(nvar+1));


    newcst = pluto_constraints_alloc(cst->nrows, CST_WIDTH);

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

    q=0;
    for (j=0; j<cst->ncols; j++) {
        if (!redun[j])  {
            for (i=0; i<cst->nrows; i++) {
                newcst->val[i][q] = cst->val[i][j];
            }
            q++;
        }
    }
    newcst->nrows = cst->nrows;
    newcst->ncols = q;

    /* Add upper bounds for transformation coefficients */
    int ub = get_coeff_upper_bound(prog);

    /* Putting too small an upper bound can prevent useful transformations;
     * also, note that an upper bound is added for all statements globally due
     * to the lack of an easy way to determine bounds for each coefficient to
     * prevent spurious transformations that involve shifts proportional to
     * loop bounds
     */
    if (ub >= 10)   {
        for (i=0; i<newcst->ncols-npar-1-1; i++)  {
            IF_DEBUG2(printf("Adding upper bound %d for transformation coefficients\n", ub););
            pluto_constraints_add_ub(newcst, npar+1+i, ub);
        }
    }

    /* Reverse the variable order for stmts */
    PlutoMatrix *perm_mat = pluto_matrix_alloc(newcst->ncols, newcst->ncols);
    PlutoMatrix *newcstmat = pluto_matrix_alloc(newcst->nrows, newcst->ncols);

    for (i=0; i<newcst->ncols; i++) {
        bzero(perm_mat->val[i], sizeof(int)*newcst->ncols);
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
    /* matrix_print(stdout, newcst->val, newcst->nrows, newcst->ncols); */
    /* matrix_print(stdout, newcstmat, newcst->nrows, newcst->ncols); */

    /* Save it so that it can be put back and freed correctly */
    int **save = newcst->val;
    newcst->val = newcstmat->val;

    // IF_DEBUG(dump_poly(newcst));
    sol = pluto_constraints_solve(newcst, DO_NOT_ALLOW_NEGATIVE_COEFF);
    /* Put it back so that it can be freed correctly */
    newcst->val = save;

    pluto_matrix_free(newcstmat);

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

        fsol = (int *)malloc(cst->ncols*sizeof(int));
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

    IF_DEBUG(printf("Cutting between SCC id %d and id %d\n", scc1, scc2));

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

    IF_DEBUG(printf("Cutting between all SCCs\n"));

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

    int i, j;

    int num_new_carried = 0;

    /* First time, cut between SCCs of different dimensionalities */
    if (cut_scc_dim_based(prog,ddg))   {
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
    int *bestsol;
    PlutoConstraints *basecst, *nzcst;
    PlutoConstraints *currcst;

    int ndeps = prog->ndeps;
    int nstmts = prog->nstmts;
    Stmt **stmts = prog->stmts;
    Dep **deps = prog->deps;
    int nvar = prog->nvar;
    int npar = prog->npar;

    IF_DEBUG(fprintf(stdout, "Finding hyperplanes: max %d\n", max_sols));

    assert(max_sols >= 0);

    if (max_sols == 0) return 0;

    basecst = get_permutability_constraints(deps, ndeps, prog);

    num_sols_found = 0;
    /* We don't expect to add a lot to basecst - just ortho constraints
     * and trivial soln avoidance constraints */
    currcst = pluto_constraints_alloc(basecst->nrows+nstmts+nvar*nstmts, CST_WIDTH);

    do{
        pluto_constraints_copy(currcst, basecst);
        nzcst = get_non_trivial_sol_constraints(prog, loop_search_mode);
        pluto_constraints_add(currcst, nzcst);
        pluto_constraints_free(nzcst);

        PlutoConstraints *indcst = get_linear_ind_constraints(prog, currcst, lin_ind_mode);

        if (indcst->nrows == 0) {
            /* If you don't have any independence constraints, we would end up finding
             * the same solution that was found earlier; so we won't find anything
             * new */
            bestsol = NULL;
        }else{
            pluto_constraints_add(currcst, indcst);
            IF_DEBUG2(printf("Solving for %d solution\n", num_sols_found+1));
            IF_DEBUG2(pluto_constraints_pretty_print(stdout, currcst));
            bestsol = pluto_prog_constraints_solve(currcst, prog);
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


    pluto_constraints_free(basecst);
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
                            fscanf(precut, "%d", &stmts[i]->trans->val[stmts[i]->trans->nrows-1][j]);
                        }else{
                            stmts[i]->trans->val[stmts[i]->trans->nrows-1][j] = 0;
                        }
                    }
                    for (j=0; j<npar; j++)    {
                        fscanf(precut, "%d", &ignore);
                        stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar] = 0;
                    }
                    /* Constant part */
                    fscanf(precut, "%d", &stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar]);

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


void pluto_print_depsat_vectors(Dep **deps, int ndeps, int levels)
{
    int i, j;

    printf("\nSatisfaction vectors for transformed program\n");

    for (i=0; i<ndeps; i++) {
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
        printf(") Satisfied: %d, Sat level: %d\n", deps[i]->satisfied, deps[i]->satisfaction_level);

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
    }
}


/* Pad statement domains to maximum domain depth to make it easier to construct
 * scheduling constraints. These will be removed before autopoly returns.
 * Also, corresponding dimensions from ILP space will be removed before ILP
 * calls
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
        for (i=0; i<prog->nstmts; i++)    {
            PlutoConstraints *copy = 
                pluto_constraints_alloc(2*prog->stmts[i]->domain->nrows, prog->stmts[i]->domain->ncols);
            pluto_constraints_copy(copy, prog->stmts[i]->domain);
            for (j=0; j<prog->stmts[i]->dim_orig; j++)    {
                fourier_motzkin_eliminate(copy, 0);
            }
            assert(copy->ncols == npar+1);
            count += copy->nrows;

            if (count <= prog->nstmts*npar)    {
                pluto_constraints_add(context, copy);
            }else{
                pluto_constraints_free(copy);
                break;
            }
            pluto_constraints_free(copy);
        }
        pluto_constraints_simplify(context);
        if (options->debug) {
            printf("Global constraint context\n");
            pluto_constraints_print(stdout, context );
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


/* Top-level automatic transformation algoritm */
int pluto_auto_transform(PlutoProg *prog)
{
    int nsols, i, j;
    int num_ind_sols, depth;
    bool lin_ind_mode;
    bool loop_search_mode;

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

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
    prog->num_hyperplanes = 0;
    prog->hProps = NULL;

    /* The number of independent solutions required for the deepest 
     * statement */
    nsols = 0;
    for (i=0; i<nstmts; i++)    {
        nsols = PLMAX(nsols, stmts[i]->dim);
    }

    num_ind_sols = 0;
    depth=0;

    if (precut(prog, ddg, depth))   {
        /* Distributed based on .fst or .precut file (customized user-supplied
         * fusion structure */
        num_ind_sols = pluto_get_max_ind_hyps(prog);
        printf("[Pluto] Forced custom fusion structure from .fst/.precut\n");
        IF_DEBUG(fprintf(stdout, "%d ind solns in .precut file\n", 
                    num_ind_sols));
    }else{
        if (options->fuse == SMART_FUSE)    {
            cut_scc_dim_based(prog,ddg);
        }
    }

    do{
        int sols_found;

        if (options->fuse == NO_FUSE)   {
            ddg_compute_scc(prog);
            cut_all_sccs(prog, ddg);
        }

        sols_found = find_permutable_hyperplanes(prog, lin_ind_mode, 
                loop_search_mode, PLMAX(1,nsols-num_ind_sols));

        IF_DEBUG(fprintf(stdout, "Level: %d; \t%d hyperplanes found\n", 
                    depth, sols_found));
        IF_DEBUG2(pluto_transformations_pretty_print(prog));
        num_ind_sols = pluto_get_max_ind_hyps(prog);

        if (sols_found > 0) {
            for (j=0; j<sols_found; j++)      {
                /* Mark dependences satisfied by this solution */
                dep_satisfaction_update(prog,
                        stmts[0]->trans->nrows-sols_found+j);
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
                    IF_DEBUG2(printf("Switching to incremental ortho constr mode\n"););
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
                    fprintf(stdout, "[Pluto] working with original (identity) transformation (if they exist)\n");
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

    denormalize_domains(prog);

    //pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);

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

void pluto_transformations_pretty_print(const PlutoProg *prog)
{
    int nstmts, i, j;

    nstmts = prog->nstmts;

    for (i=0; i<nstmts; i++) {
        Stmt *stmt = prog->stmts[i];
        fprintf(stdout, "T(S%d): ", i+1);
        int level;
        printf("(");
        for (level=0; level<stmt->trans->nrows; level++) {
            char **vars = malloc((stmt->dim+prog->npar)*sizeof(char *));
            for (j=0; j<stmt->dim; j++) {
                vars[j] = stmt->iterators[j];
            }
            for (j=0; j<prog->npar; j++) {
                vars[stmt->dim+j] = prog->params[j];
            }
            if (level > 0) printf(", ");
            pretty_print_affine_function(stdout, stmt->trans->val[level], 
                    stmt->dim+prog->npar, vars);
            free(vars);
        }
        printf(")\n");

        pluto_matrix_print(stdout, stmt->trans);

        printf("loop types (");
        for (level=0; level<stmt->trans->nrows; level++) {
            if (level > 0) printf(", ");
            if (stmt->hyp_types[level] == H_SCALAR) printf("scalar");
            else if (stmt->hyp_types[level] == H_LOOP) printf("loop");
            else if (stmt->hyp_types[level] == H_TILE_SPACE_LOOP) printf("tloop");
            else printf("unknown");
        }
        printf(")\n\n");
    }
}


/* List properties of newly found hyperplanes */
void pluto_print_hyperplane_properties(const PlutoProg *prog)
{
    int j, numH;
    HyperplaneProperties *hProps;

    hProps = prog->hProps;
    numH = prog->num_hyperplanes;

    if (numH == 0)  {
        fprintf(stdout, "No hyperplanes\n");
    }

    /* Note that loop properties are calculated for each dimension in the
     * transformed space (common for all statements) */
    for (j=0; j<numH; j++)  {
        fprintf(stdout, "t%d --> ", j+1);
        switch (hProps[j].dep_prop)    {
            case PARALLEL:
                fprintf(stdout, "parallel ");
                break;
            case SEQ:
                fprintf(stdout, "serial   ");
                break;
            case PIPE_PARALLEL:
                fprintf(stdout, "fwd_dep  ");
                break;
            default:
                fprintf(stdout, "unknown  ");
                break;
        }
        switch (hProps[j].type) {
            case H_LOOP:
                fprintf(stdout, "loop  ");
                break;
            case H_SCALAR:
                fprintf(stdout, "scalar");
                break;
            case H_TILE_SPACE_LOOP:
                fprintf(stdout, "tLoop ");
                break;
            default:
                fprintf(stdout, "unknown  ");
                // assert(0);
                break;
        }
        fprintf(stdout, " (band %d)", hProps[j].band_num);
        fprintf(stdout, hProps[j].unroll? "ujam":"no-ujam"); 
        fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}


/*
 * Pretty prints a one-dimensional affine function
 * ndims: number of variables
 * func should have ndims+1 elements (affine function)
 * vars: names of the variables; if either of them are NULL, x0, x1, ... are used
 */
void pretty_print_affine_function(FILE *fp, int *func, int ndims, char **vars)
{
    char *var[ndims];
    int j;

    for (j=0; j<ndims; j++)  {
        if (vars[j]) {
            var[j] = strdup(vars[j]);
        }else{
            var[j] = malloc(5);
            sprintf(var[j], "x%d", j+1);
        }
    }

    int sign_flag = 0;
    for (j=0; j<ndims; j++)   {
        if (func[j] == 1)  {
            if (sign_flag) fprintf(fp, "+");
            fprintf(fp, "%s", var[j]);
        }else if (func[j] == -1)  {
            if (sign_flag) fprintf(fp, "-");
            fprintf(fp, "%s", var[j]);
        }else if (func[j] != 0)  {
            if (sign_flag) fprintf(fp, "+");
            if (func[j] > 0) {
                fprintf(fp, "%d%s", func[j], var[j]);
            }else{
                fprintf(fp, "(%d%s)", func[j], var[j]);
            }
        }
        if (func[j] != 0) sign_flag = 1;
    }
    if (func[ndims] >= 1)  {
        if (sign_flag) fprintf(fp, "+");
        fprintf(fp, "%d", func[ndims]);
    }else{
        if (!sign_flag) fprintf(fp, "%d", func[ndims]);
    }

    for (j=0; j<ndims; j++)  {
        free(var[j]);
    }
}


void pluto_transformations_print(const PlutoProg *prog)
{
    int i;

    for (i=0; i<prog->nstmts; i++)    {
        printf("T_(S%d) \n", prog->stmts[i]->id+1);
        pluto_matrix_print(stdout, prog->stmts[i]->trans);
    }
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
                pluto_constraints_add_dim(cst_l, stmt->dim);
            }
            pluto_constraints_add(dep->dpolytope, cst_l);
            pluto_constraints_free(cst_l);
        }
        if (stmts[dep->dest] == stmt) {
            PlutoConstraints *cst_l = pluto_constraints_dup(cst);
            Stmt *sstmt = stmts[dep->src];
            for (c=0; c<sstmt->dim; c++) {
                pluto_constraints_add_dim(cst_l, 0);
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
                pluto_constraints_add_dim(cst_l, stmt->dim);
            }
            pluto_constraints_add(dep->dpolytope, cst_l);
            pluto_constraints_free(cst_l);
        }
        if (stmts[dep->dest] == stmt) {
            PlutoConstraints *cst_l = pluto_constraints_dup(cst);
            Stmt *sstmt = stmts[dep->src];
            for (c=0; c<sstmt->dim; c++) {
                pluto_constraints_add_dim(cst_l, 0);
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
            int is_const_ub, is_const_lb, ub, lb;
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
