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

#include <sys/time.h>

#include "pluto.h"
#include "math_support.h"
#include "constraints.h"
#include "post_transform.h"
#include "program.h"
#include "transforms.h"
#include "ddg.h"
#include "version.h"

void pluto_print_colours(int *colour,PlutoProg *prog);
bool* innermost_dep_satisfaction_dims(PlutoProg *prog, bool *tile_preventing_deps);
bool colour_scc(int scc_id, int *colour, int c, int stmt_pos, int pv, PlutoProg *prog);


bool dep_satisfaction_test(Dep *dep, PlutoProg *prog, int level);

int get_num_unsatisfied_deps(Dep **deps, int ndeps);
int get_num_unsatisfied_inter_stmt_deps(Dep **deps, int ndeps);
int get_num_unsatisfied_inter_scc_deps(PlutoProg *prog);

int pluto_diamond_tile(PlutoProg *prog);

static double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

/*
 * Returns the number of (new) satisfied dependences at this level
 *
 * NOTE: for every unsatisfied dependence, this function tests if the entire
 * dependence has been satisfied at 'level'
 *
 */
int dep_satisfaction_update(PlutoProg *prog, int level)
{
    int i;
    int num_new_carried;

    int ndeps = prog->ndeps;
    Dep **deps = prog->deps;

    num_new_carried = 0;

    IF_DEBUG(printf("[pluto] dep_satisfaction_update (level %d)\n", level););

    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];
        if (!dep_is_satisfied(dep)) {
            dep->satisfied = dep_satisfaction_test(dep, prog, level);
            if (dep->satisfied) {
                IF_MORE_DEBUG(printf("[pluto] dep_satisfaction_update: dep %d satisfied\n", i+1););
                if (!IS_RAR(dep->type)) num_new_carried++;
                dep->satisfaction_level = level;
            }
        }
    }
    IF_DEBUG(printf("\t %d dep(s) satisfied\n", num_new_carried););

    return num_new_carried;
}

/* Check whether all deps are satisfied */
int deps_satisfaction_check(PlutoProg *prog)
{
    int i;

    for (i=0; i<prog->ndeps; i++) {
        if (IS_RAR(prog->deps[i]->type)) continue;
        if (!dep_is_satisfied(prog->deps[i])) {
            return false;
        }
    }
    return true;
}

void pluto_dep_satisfaction_reset(PlutoProg *prog)
{
    int i;

    IF_DEBUG(printf("[pluto] pluto_dep_satisfaction_reset\n"););

    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i]->satisfied = false;
        prog->deps[i]->satisfaction_level = -1;
    }
}

/*
 * Conservative but powerful enough: until a dependence has been completely
 * satisfied (a level at which it is completely satisifed), a non-zero
 * dependence component would set satvec for that level to one
 */
void pluto_compute_dep_satisfaction(PlutoProg *prog)
{
    int i, level;

    IF_DEBUG(printf("[pluto] pluto_compute_dep_satisfaction\n"););

    pluto_dep_satisfaction_reset(prog);
    
    for (i=0; i<prog->num_hyperplanes; i++) {
        dep_satisfaction_update(prog, i);
    }

    /* Create and set satisfaction vectors */
    for (i=0; i<prog->ndeps; i++) {
        Dep *dep = prog->deps[i];
        if (IS_RAR(dep->type)) continue;


        if (dep->satvec) free(dep->satvec);
        dep->satvec = (int *) malloc(prog->num_hyperplanes*sizeof(int));

        /* Direction vectors should be available */
        assert(dep->dirvec != NULL);

        for (level=0; level<prog->num_hyperplanes; level++) {
            if (dep->dirvec[level] != DEP_ZERO &&
                    (dep->satisfaction_level >= level || dep->satisfaction_level == -1)) {
                dep->satvec[level] = 1;
            }else{
                dep->satvec[level] = 0;
            }
        }
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
        if (deps[i].src != deps[i].dest) {
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
            IF_DEBUG2(printf("[pluto_get_coeff_bound_for_cone_complement] Adding lower bound %d for stmt dim coefficients\n", -4););
            pluto_constraints_add_lb(cst, npar+1+s*(nvar+1)+i, -4);
        }
        /* Translation coefficients need not be negative */
        IF_DEBUG2(printf("[pluto_get_coeff_bound_for_cone_complement] Adding lower bound %d for stmt translation coefficient\n", 0););
        pluto_constraints_add_lb(cst, npar+1+s*(nvar+1)+nvar, 0);
    }
    return cst;
}


PlutoMatrix* construct_cplex_objective(const PlutoConstraints *cst, const PlutoProg *prog)
{
    int npar = prog->npar;
    int nvar = prog->nvar;
    int i, j, k;

    PlutoMatrix *obj = pluto_matrix_alloc(1, cst->ncols-1);
    pluto_matrix_set(obj, 0);

    /* u */
    for (j=0; j<npar; j++) {
        obj->val[0][j] = 5*5*nvar*prog->nstmts;
    }
    /* w */
    obj->val[0][npar] = 5*nvar*prog->nstmts;

    for (i=0, j=npar+1; i<prog->nstmts; i++) {
        for (k=j; k<j+prog->stmts[i]->dim_orig; k++) {
            obj->val[0][k] = (nvar+2)*(prog->stmts[i]->dim_orig-(k-j));
        }
        /* constant shift */
        obj->val[0][k] = 1;
        j += prog->stmts[i]->dim_orig + 1;
    }
    return obj;
}


/*
 * This calls pluto_constraints_lexmin, but before doing that does some preprocessing
 * - removes variables that we know will be assigned 0 - also do some
 *   permutation/substitution of variables
 */
int64 *pluto_prog_constraints_lexmin(PlutoConstraints *cst, PlutoProg *prog)
{
    Stmt **stmts;
    int i, j, k;
    int nstmts, nvar, npar, del_count;
    int64 *sol, *fsol;
    PlutoConstraints *newcst;
    double t_start;

    stmts = prog->stmts;
    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;

    assert(cst->ncols - 1 == CST_WIDTH - 1);

    /* Remove redundant variables - that don't appear in your outer loops */
    int redun[npar+1+nstmts*(nvar+1)+1];
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

    del_count = 0;
    newcst = pluto_constraints_dup(cst);
    for (j = 0; j < cst->ncols-1; j++) {
        if (redun[j]) {
            pluto_constraints_remove_dim(newcst, j-del_count);
            del_count++;
        }
    }

    /* Permute the constraints so that if all else is the same, the original
     * hyperplane order is preserved (no strong reason to do this) */
    /* We do not need to permute in case of pluto-lp-dfp */
    if (!options->dfp){
    j = npar + 1;
    for (i=0; i<nstmts; i++)    {
        for (k=j; k<j+(stmts[i]->dim_orig)/2; k++) {
            pluto_constraints_interchange_cols(newcst, k, j + (stmts[i]->dim_orig - 1 - (k-j)));

        }
        j += stmts[i]->dim_orig+1;
    }

    }
    IF_DEBUG(printf("[pluto] pluto_prog_constraints_lexmin (%d variables, %d constraints)\n",
                cst->ncols-1, cst->nrows););

    /* Solve the constraints using chosen solvers*/
    if (options->islsolve) {
        t_start = rtclock(); 
        sol = pluto_constraints_lexmin_isl(newcst, DO_NOT_ALLOW_NEGATIVE_COEFF);
        prog->mipTime += rtclock()-t_start;
    }
#ifdef GLPK
    else if (options->glpk || options->lp || options->dfp || options->gurobi) {
         
        double **val = NULL;
        int **index = NULL;
        int num_ccs, nrows;
        num_ccs = 0;

        PlutoMatrix *obj = construct_cplex_objective(newcst, prog);

        if (options->lp) {
            nrows = newcst->ncols-1-npar-1;
            populate_scaling_csr_matrices_for_pluto_program(&index, &val, nrows, prog);
            num_ccs = prog->ddg->num_ccs;
        }

        t_start = rtclock(); 
        if (options->glpk) {
        sol = pluto_prog_constraints_lexmin_glpk(newcst, obj, val, index, npar, num_ccs);
        } 
        else if (options->gurobi) {
            sol = pluto_prog_constraints_lexmin_gurobi(newcst, obj, val, index, npar, num_ccs);
        }
        prog->mipTime += rtclock()-t_start;

        pluto_matrix_free(obj);
        if (options->lp) {
            for (i=0; i<nrows; i++) {
                free(val[i]);
                free(index[i]);
            }
            free(val);
            free(index);
        }
    }
#endif
    else {
        t_start = rtclock(); 
        sol = pluto_constraints_lexmin_pip(newcst, DO_NOT_ALLOW_NEGATIVE_COEFF);
        prog->mipTime += rtclock()-t_start;
    }
    /* sol = pluto_constraints_lexmin(newcst, DO_NOT_ALLOW_NEGATIVE_COEFF); */
    /* print_polylib_visual_sets("csts", newcst); */


    fsol = NULL;
    if (sol) {
        int k1, k2, q;
        int64 tmp;
        /* Permute the solution in line with the permuted cst */
        if(!options->dfp){
            j = npar + 1;
            for (i=0; i<nstmts; i++)    {
                for (k=j; k<j+(stmts[i]->dim_orig)/2; k++) {
                    k1 = k;
                    k2 = j + (stmts[i]->dim_orig - 1 - (k-j));
                    tmp = sol[k1];
                    sol[k1] = sol[k2];
                    sol[k2] = tmp;
                }
                j += stmts[i]->dim_orig+1;
            }
        }

        fsol = (int64 *) malloc((cst->ncols-1)*sizeof(int64));

        /* Fill the soln with zeros for the redundant variables */
        q = 0;
        for (j=0; j<cst->ncols-1; j++) {
            fsol[j] = redun[j]? 0: sol[q++];
        }
        free(sol);
    }

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
        IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes-1););
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
        IF_DEBUG(printf("\tonly one SCC\n"));
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
    IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes-1););
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
 *
 * (Due to the algorithm used for computing SCCs, there are no edges 
 * between SCC<i> and SCC<j> where i > j)
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
        IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes-1););
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

    IF_DEBUG(printf("[pluto] get_linear_ind_constraints\n"););

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
 * See sub-functions for hyp_search_mode and lin_ind_mode
 *
 * If all statements already have enough linearly independent solutions, no
 * independence constraints will be generated, and since no non-trivial
 * solution constraints are added in such a case, the trivial zero solution will 
 * end up being returned.
 * */
int find_permutable_hyperplanes(PlutoProg *prog, bool hyp_search_mode, 
        int max_sols, int band_depth)
{
    int num_sols_found, j, k;
    int64 *bestsol;
    PlutoConstraints *basecst, *nzcst, *boundcst;
    PlutoConstraints *currcst;

    int nstmts = prog->nstmts;
    Stmt **stmts = prog->stmts;
    int nvar = prog->nvar;
    int npar = prog->npar;

    IF_DEBUG(fprintf(stdout, "[pluto] find_permutable_hyperplanes: max solution(s): %d; band depth: %d\n", max_sols, band_depth));

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
        nzcst = get_non_trivial_sol_constraints(prog, hyp_search_mode);
        pluto_constraints_add(currcst, nzcst);
        pluto_constraints_free(nzcst);

        PlutoConstraints *indcst = get_linear_ind_constraints(prog, currcst, hyp_search_mode);
        // print_polylib_visual_sets("ind", indcst);
        IF_DEBUG2(printf("linear independence constraints\n"));
        IF_DEBUG2(pluto_constraints_pretty_print(stdout, indcst););

        if (indcst->nrows == 0) {
            /* If you don't have any independence constraints, we would end 
             * up finding the same solution that was found earlier; so we 
             * won't find anything new */
            IF_DEBUG(printf("No linearly independent rows\n"););
            bestsol = NULL;
        }else{
            pluto_constraints_add(currcst, indcst);
            IF_DEBUG(printf("[pluto] (Band %d) Solving for hyperplane #%d\n", band_depth+1, num_sols_found+1));
            // IF_DEBUG2(pluto_constraints_pretty_print(stdout, currcst));
            bestsol = pluto_prog_constraints_lexmin(currcst, prog);
        }
        pluto_constraints_free(indcst);

        if (bestsol != NULL)    {
            IF_DEBUG(fprintf(stdout, "[pluto] find_permutable_hyperplanes: found a hyperplane\n"));
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
            IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes-1););
        }else{
            IF_DEBUG(fprintf(stdout, "[pluto] find_permutable_hyperplanes: No hyperplane found\n"));
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

            assert(ignore == prog->nstmts);

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

                    pluto_stmt_add_hyperplane(stmts[i], H_UNKNOWN, stmts[i]->trans->nrows);

                    for (j=0; j<nvar; j++)    {
                        if (stmts[i]->is_orig_loop[j])  {
                            fscanf(precut, "%lld", &stmts[i]->trans->val[stmts[i]->trans->nrows-1][j]);
                        }else{
                            stmts[i]->trans->val[stmts[i]->trans->nrows-1][j] = 0;
                        }
                    }
                    for (j=0; j<npar; j++)    {
                        fscanf(precut, "%d", &ignore);
                        stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar+j] = 0;
                    }
                    /* Constant part */
                    fscanf(precut, "%lld", &stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar+npar]);
                    if (get_loop_type(stmts[i], stmts[i]->trans->nrows-1)
                            == H_LOOP)  {
                        stmts[i]->hyp_types[stmts[i]->trans->nrows-1] = H_LOOP;
                    }else{
                        stmts[i]->hyp_types[stmts[i]->trans->nrows-1] = H_SCALAR;
                    }

                }

                /* Number of tiling levels */
                fscanf(precut, "%d", &ignore);

                /* FIXME: to tile or not is specified depth-wise, why? Just
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
                        prog->hProps[prog->num_hyperplanes-1].type = H_LOOP;
                    }else{
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

    /* OUTDATED */
    assert(0);

    IF_DEBUG(printf("[pluto] pluto_detect_transformation_properties\n"););

    if (prog->nstmts == 0) return;

    HyperplaneProperties *hProps = prog->hProps;

    assert(prog->num_hyperplanes == stmts[0]->trans->nrows);

    // pluto_deps_print(stdout, prog);

    /* First compute satisfaction levels */
    pluto_compute_dep_directions(prog);
    pluto_compute_dep_satisfaction(prog);

    band = 0;
    num_loops_in_band = 0;
    int bandStart = 0;

    for (level=0; level < prog->num_hyperplanes; ) {
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
                /* Dependences satisfied at scalar dimensions earlier are fine (even 
                 * if at dimensions in the same band */
                if (deps[i]->satisfaction_level < level && 
                        hProps[deps[i]->satisfaction_level].type == H_SCALAR) continue;
                /* Dependences satisfied within the band or not satisfied yet
                 * should have non-negative components */
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
                 *
                 * CORRECTION: under partial satisfaction, this can happen.
                 */
                if (num_loops_in_band == 0) {
                    fprintf(stdout, "[pluto] Unfortunately, the transformation computed has violated a dependence.\n");
                    fprintf(stdout, "\tPlease make sure there is no inconsistent/illegal .fst file in your working directory.\n");
                    fprintf(stdout, "\tIf not, this usually is a result of a bug in the dependence tester,\n");
                    fprintf(stdout, "\tor a bug in Pluto's auto transformation.\n");
                    fprintf(stdout, "\tPlease send this input file to the author if possible.\n");
                    pluto_transformations_pretty_print(prog);
                    pluto_print_hyperplane_properties(prog);
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
    }

    if (num_loops_in_band == 1) {
        if (hProps[level-1].dep_prop == PIPE_PARALLEL)
            hProps[level-1].dep_prop = SEQ;
    }

    /* 
     * This functionality is obsolete since Ploop based functions are now used
     * Permutable bands of loops could have inner parallel loops; they 
     * all have been detected as fwd_dep (except the outer parallel one of a band); 
     * we just modify those to parallel */
    for (i=0; i<prog->num_hyperplanes; i++)  {
        for (j=0; j<prog->ndeps; j++) {
            if (IS_RAR(deps[j]->type)) continue;
            if (deps[j]->satisfaction_level >= i && deps[j]->dirvec[i] != DEP_ZERO) break;
        }

        if (j==prog->ndeps)   {
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
    int i, j, ndeps, nlevels;
    Dep **deps;
   
    deps = prog->deps;
    ndeps = prog->ndeps;
    nlevels = prog->num_hyperplanes;

    printf("\nDirection vectors for transformed program\n");

    for (i=0; i<ndeps; i++) {
        printf("Dep %d: S%d to S%d: ", i+1, deps[i]->src+1, deps[i]->dest+1);
        printf("(");
        for (j=0; j<nlevels; j++) {
            printf("%c, ", deps[i]->dirvec[j]);
        }
        assert(deps[i]->satvec != NULL);
        printf(") satisfied: %s, satvec: (", 
                deps[i]->satisfied? "yes":" no");
        for (j=0; j<nlevels; j++) {
            printf("%d, ", deps[i]->satvec[j]);
        }
        printf("), sat level: %d", deps[i]->satisfaction_level);

        for (j=0; j<nlevels; j++) {
            if (deps[i]->dirvec[j] > 0)  {
                break;
            }
            /* Just because this is not printed does not mean a dependence has
             * not been violated; it could still be */
            if (deps[i]->dirvec[j] < 0) {
                printf("Dep %d violated: S%d to S%d\n", i, deps[i]->src+1, deps[i]->dest+1);
                printf("%d %d\n", deps[i]->satisfaction_level, deps[i]->satisfied);
            }else if (deps[i]->dirvec[j] < 0) {
                printf("Dep %d violated: S%d to S%d\n", i, deps[i]->src+1, deps[i]->dest+1);
                printf("%d %d\n", deps[i]->satisfaction_level, deps[i]->satisfied);
            }
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
                pluto_constraints_add(context, copy);
                pluto_constraints_free(copy);
            }else{
                pluto_constraints_free(copy);
                break;
            }
        }
        pluto_constraints_simplify(context);
        IF_DEBUG(printf("[pluto] parameter context from domains\n"););
        IF_DEBUG(pluto_constraints_compact_print(stdout, context););

        /* Add context to every dep polyhedron */
        for (i=0; i<prog->ndeps; i++) {
            PlutoConstraints *bounding_poly = 
                pluto_constraints_dup(prog->deps[i]->dpolytope);

            for (k=0; k<context->nrows; k++) {
                pluto_constraints_add_inequality(bounding_poly);

                /* Already initialized to zero */

                for (j=0; j<npar+1; j++){
                    bounding_poly->val[bounding_poly->nrows-1][j+bounding_poly->ncols-(npar+1)] = 
                        context->val[k][j];
                }
            }
            /* Update reference, add_row can resize */
            pluto_constraints_free(prog->deps[i]->bounding_poly);
            prog->deps[i]->bounding_poly = bounding_poly;
        }
        pluto_constraints_free(context);
    }else{
        IF_DEBUG(printf("\tNo global context\n"));
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
     * values */
    bool *neg = malloc(sizeof(bool)*npar);
    for (k=0; k<prog->ndeps; k++) {
        Dep *dep = prog->deps[k];
        PlutoConstraints *dpoly = dep->dpolytope;

        bzero(neg, npar*sizeof(bool));

        if (dpoly->nrows == 0) continue;

        for (j=2*nvar; j<2*nvar+npar; j++)  {
            int min = dpoly->val[0][j];
            int max = dpoly->val[0][j];
            for (i=1; i<dpoly->nrows; i++)  {
                min = PLMIN(dpoly->val[i][j], min);
                max = PLMAX(dpoly->val[i][j], max);
            }

            if (min < 0 && max <= 0)    {
                neg[j-2*nvar] = true;
                IF_DEBUG(printf("Dep %d has negative coeff's for parameter %s\n",
                            dep->id, prog->params[j-2*nvar]));
            }
        }

        /* For parameters appearing with negative coefficients in upper bounds */
        for (j=0; j<npar; j++)  {
            if (neg[j])   {
                pluto_constraints_add_inequality(dep->bounding_poly);
                dep->bounding_poly->val[dep->bounding_poly->nrows-1][2*nvar+j] = 1;
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

/*
 * Finds domain face that allows concurrent start (for diamond tiling)
 * FIXME: iteration space boundaries are assumed to be corresponding to
 * rectangular ones
 *
 * Returns: matrix with row i being the concurrent start face for Stmt i
 */
PlutoMatrix *get_face_with_concurrent_start(PlutoProg *prog, Band *band)
{
    PlutoConstraints *bcst;
    int s, _s, j, nz, nvar, npar;
    PlutoMatrix *conc_start_faces;

    IF_DEBUG(printf("[pluto] get_face_with_concurrent_start for band\n\t"););
    IF_DEBUG(pluto_band_print(band););

    npar = prog->npar;
    nvar = prog->nvar;

    PlutoConstraints *fcst = get_feautrier_schedule_constraints(prog, 
            band->loop->stmts, band->loop->nstmts);

    bcst = get_coeff_bounding_constraints(prog);
    pluto_constraints_add(fcst, bcst);
    pluto_constraints_free(bcst);

    int64 *sol = pluto_prog_constraints_lexmin(fcst, prog);
    pluto_constraints_free(fcst);

    if (!sol) {
        IF_DEBUG(printf("[pluto] get_face_with_concurrent_start: no valid 1-d schedules \n"););
        return NULL;
    }

    conc_start_faces = pluto_matrix_alloc(band->loop->nstmts, nvar+npar+1);
    pluto_matrix_set(conc_start_faces, 0);

    for (s=0, _s=0; s<prog->nstmts; s++) {
        if (!pluto_stmt_is_member_of(s, band->loop->stmts, 
                    band->loop->nstmts)) continue;
        assert(_s <= band->loop->nstmts-1);
        for (j=0; j<nvar; j++) {
            conc_start_faces->val[_s][j] = sol[npar+1+s*(nvar+1)+j];
        }
        conc_start_faces->val[_s][nvar+npar] = sol[npar+1+s*(nvar+1)+nvar];
        _s++;
    }
    free(sol);

    IF_DEBUG(printf("[pluto] get_face_with_concurrent_start: 1-d schedules\n"););
    for (s=0; s<band->loop->nstmts; s++) {
        IF_DEBUG(printf("\tf(S%d) = ", band->loop->stmts[s]->id+1););
        IF_DEBUG(pluto_affine_function_print(stdout, conc_start_faces->val[s], nvar+npar,
                    band->loop->stmts[s]->domain->names););
        IF_DEBUG(printf("\n"););
    }

    /* 1-d schedule should be parallel to an iteration space boundary 
     * FIXME: assuming canonical boundaries */
    for (s=0; s<band->loop->nstmts; s++) {
        nz = 0;
        for (j=0; j<nvar; j++) {
            if (conc_start_faces->val[s][j]) nz++;
        }
        if (nz != 1) break;
    }

    if (s < band->loop->nstmts) {
        pluto_matrix_free(conc_start_faces);
        IF_DEBUG(printf("[pluto] No iteration space faces with concurrent start for all statements\n"););
        return NULL;
    }

    IF_DEBUG(printf("[pluto] faces with concurrent start found for all statements\n"););

    return conc_start_faces;
}


/*
 * Find hyperplane inside the cone  of previously found hyperplanes 
 * and the face allowing concurrent start
 *
 * conc_start_faces[i]: concurrent start face for statement $i$
 *
 * evict_pos: position of the hyperplane to be evicted by the one that will
 * enable concurrent start
 *
 * cone_complement_pos: in case of partial concurrent start, the 
 * hyperplane that will form the cone with the conc start hyperplane
 *
 * cone_complement_hyps will set to the cone complement hyperplanes found
 * for statements in the band
 */
int find_cone_complement_hyperplane(Band *band, PlutoMatrix *conc_start_faces, int evict_pos, 
        int cone_complement_pos, PlutoConstraints *basecst, PlutoProg *prog, 
        PlutoMatrix **cone_complement_hyps)
{
    int i, s, j, k, lambda_k, nstmts, nvar, npar;
    int64 *bestsol;
    PlutoConstraints *con_start_cst, *lastcst;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = band->loop->nstmts;

    IF_DEBUG(printf("[pluto] find_cone_complement_hyperplane for band\n\t"););
    IF_DEBUG(pluto_band_print(band););

    /* lastcst is the set of additional constraints */
    lastcst = pluto_constraints_alloc(2*nvar*nstmts, 
            (npar+1+prog->nstmts*(nvar+1)+1) + nvar*nstmts);

    /* all lambdas >=1 */
    for (i=0; i<nstmts; i++) {
        int stmt_offset = npar+1+prog->nstmts*(nvar+1)+i*nvar;
        for (j=0; j<nvar; j++)  {
            pluto_constraints_add_inequality(lastcst);
            lastcst->val[lastcst->nrows-1][stmt_offset+j] =1;
            lastcst->val[lastcst->nrows-1][lastcst->ncols-1] = -1;
        }
    }

    /* Now, add the constraints for the new hyperplane to be in the cone
     * of the face and the negatives of the hyperplanes already found
     * (excluding the one being evicted: at `evict_pos') */
    for (s=0; s<nstmts; s++) {
        Stmt *stmt = band->loop->stmts[s];
        int stmt_offset1= npar+1+stmt->id*(nvar+1);
        int stmt_offset2= npar+1+prog->nstmts*(nvar+1)+s*nvar;
        for (j=0; j<nvar; j++)  {
            pluto_constraints_add_equality(lastcst);
            lastcst->val[lastcst->nrows-1][stmt_offset1+j] =1;

            lastcst->val[lastcst->nrows-1][stmt_offset2] = -(conc_start_faces->val[s][j]);

            if (options->partlbtile) {
                lastcst->val[lastcst->nrows-1][stmt_offset2+1] = 
                    stmt->trans->val[cone_complement_pos][j];
            }else{
                lambda_k = 0;
                /* Just for the band depth hyperplanes */
                for (k=band->loop->depth; k < band->loop->depth + band->width; k++){
                    if (k != evict_pos && stmt->hyp_types[k] != H_SCALAR) {
                        lastcst->val[lastcst->nrows-1][stmt_offset2+lambda_k+1] = stmt->trans->val[k][j];
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
        printf("[pluto] No concurrent start possible\n");
    }else{
        IF_DEBUG(printf("[pluto] Concurrent start possible\n"););
        for (j=0; j<nstmts; j++) {
            Stmt *stmt = band->loop->stmts[j];
            cone_complement_hyps[j] =
                pluto_matrix_alloc(1, stmt->dim+npar+1);
            for (k=0; k<nvar; k++)    {
                cone_complement_hyps[j]->val[0][k] =
                    bestsol[npar+1+ stmt->id*(nvar+1) + k];
            }
            /* No parametric shifts */
            for (k=nvar; k<nvar+npar; k++)    {
                cone_complement_hyps[j]->val[0][k] = 0;
            }
            cone_complement_hyps[j]->val[0][nvar+npar] =
                bestsol[npar+1+ stmt->id*(nvar+1) + nvar];

            IF_DEBUG(printf("\tcone_complement(S%d) = ", stmt->id+1););
            IF_DEBUG(pluto_affine_function_print(stdout, cone_complement_hyps[j]->val[0], nvar, stmt->iterators););
            IF_DEBUG(printf("\n"););
        }
        free(bestsol);
    }

    return (cone_complement_hyps[0] == NULL)? 0:1; 
}

/* 
 * Check if the d'th row for any statement in the band is parallel to the
 * face allowing concurrent start
 */
int is_concurrent_start_face(Band *band, PlutoMatrix *conc_start_faces, int d)
{
    int s;

    for (s=0; s<band->loop->nstmts; s++){
        if (pluto_vector_is_parallel(band->loop->stmts[s]->trans, d,
                conc_start_faces, s)) {
            return 1;
        }
    }
    return 0;
}

/* 
 * Find the hyperplace parallel to the concurrent start face 
 * that will be evicted; if there is none, return
 * the last hyperplane
 */
int find_hyperplane_to_be_evicted(Band *band, PlutoMatrix *conc_start_faces, PlutoProg *prog)
{
    int d;
    for (d=band->loop->depth; d < band->loop->depth + band->width-1; d++){
        if (is_concurrent_start_face(band, conc_start_faces, d)) return d;
    }
    /* Return the last one */
    return band->loop->depth + band->width - 1; 
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

/* 
 * Top-level automatic transformation algoritm 
 *
 * All dependences are reset to unsatisfied before starting
 *
 */ 
int pluto_auto_transform(PlutoProg *prog)
{
    int i, j, s, nsols, conc_start_found, depth;
    /* The maximum number of linearly independent solutions needed across all
     * statements */
    int num_ind_sols_req;

    /* The number of linearly independent solutions found (max across all
     * statements) */
    int num_ind_sols_found;
    /* Pluto algo mode -- LAZY or EAGER */
    bool hyp_search_mode;

#ifdef GLPK
    Graph* fcg;
    int *colour, nVertices;
#endif

    Stmt **stmts = prog->stmts;
    int nstmts = prog->nstmts;

    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i]->satisfied = false;
    }

    /* Create the data dependence graph */
    prog->ddg = ddg_create(prog);
    ddg_compute_scc(prog);
    for (i=0; i<prog->ddg->num_sccs; i++){
        prog->ddg->sccs[i].vertices = NULL;
    }

    Graph *ddg = prog->ddg;
    int nvar = prog->nvar;
    int npar = prog->npar;

    prog->cst_solve_time = 0.0;
    prog->cst_const_time = 0.0;
    prog->scaling_cst_sol_time = 0.0;
    prog->mipTime = 0.0;
    prog->ilpTime = 0.0;
    prog->skew_time = 0.0;
    prog->cst_write_time = 0.0;
    prog->fcg_const_time = 0.0;
    prog->fcg_update_time = 0.0;
    prog->fcg_colour_time = 0.0;
    prog->fcg_dims_scale_time = 0.0;
    prog->fcg_cst_alloc_time = 0.0;

    prog->num_lp_calls = 0;

    if (nstmts == 0)  return 0;

    PlutoMatrix **orig_trans = malloc(nstmts*sizeof(PlutoMatrix *));
    PlutoHypType **orig_hyp_types = malloc(nstmts*sizeof(PlutoHypType *));
    int orig_num_hyperplanes = prog->num_hyperplanes;
    HyperplaneProperties *orig_hProps = prog->hProps;

    /* Get rid of any existing transformation */
    for (i=0; i<nstmts; i++) {
        Stmt *stmt = prog->stmts[i];
        /* Save the original transformation */
        orig_trans[i] = stmt->trans;
        orig_hyp_types[i] = stmt->hyp_types;
        /* Pre-allocate a little more to prevent frequent realloc */
        stmt->trans = pluto_matrix_alloc(2*stmt->dim+1, stmt->dim+npar+1);
        stmt->trans->nrows = 0;
        stmt->hyp_types = NULL;
        stmt->intra_stmt_dep_cst = NULL;
    }

    normalize_domains(prog);

    hyp_search_mode = EAGER;

    prog->num_hyperplanes = 0;
    prog->hProps = NULL;

    /* The number of independent solutions required for the deepest 
     * statement */
    num_ind_sols_req = 0;
    for (i=0; i<nstmts; i++)    {
        num_ind_sols_req = PLMAX(num_ind_sols_req, stmts[i]->dim);
    }

    depth = 0;

    if (precut(prog, ddg, depth))   {
        /* Distributed based on .fst or .precut file (customized user-supplied
         * fusion structure */
        num_ind_sols_found = pluto_get_max_ind_hyps(prog);
        printf("[pluto] Forced custom fusion structure from .fst/.precut\n");
        IF_DEBUG(fprintf(stdout, "%d ind solns in .precut file\n", 
                    num_ind_sols_found));
    }else{
        num_ind_sols_found = 0;
        if (options->fuse == SMART_FUSE)    {
            cut_scc_dim_based(prog,ddg);
        }
    }
    


    /* For diamond tiling */
    conc_start_found = 0;

    if (options->dfp) {
#ifdef GLPK
        if(options->fuse == NO_FUSE) {
            ddg_compute_scc(prog);
            cut_all_sccs(prog, ddg);
        }
        compute_scc_vertices(prog->ddg);
        IF_DEBUG(printf("[Pluto] Initial DDG\n"););
        IF_DEBUG(pluto_matrix_print(stdout, prog->ddg->adj););
        /* ddg_compute_scc(prog); */
        if (!options->silent) {
            printf("[Pluto] Building fusion conflict graph\n");
        }

        nVertices = 0;
        for (i=0; i<nstmts; i++) {
            ddg->vertices[i].fcg_stmt_offset = nVertices;
            nVertices += stmts[i]->dim_orig;
        }

        colour = (int*) malloc(nVertices*sizeof(int));
        for(i=0; i<nVertices;i++){
            colour[i] = 0;
        }

        PlutoConstraints *permutecst = get_permutability_constraints(prog);
        IF_DEBUG(pluto_constraints_cplex_print(stdout,permutecst););

        /* Yet to start colouring hence the current_colour can be either 0 or 1 */
        prog->fcg = build_fusion_conflict_graph(prog,colour, nVertices, 0);


        fcg = prog->fcg;
        fcg->num_coloured_vertices = 0;
        fcg->to_be_rebuilt = false;

        IF_DEBUG(printf("[pluto] Fusion Conflict graph\n"););
        IF_DEBUG(pluto_matrix_print(stdout, fcg->adj););

        prog->total_coloured_stmts = (int*) malloc(nvar*sizeof(int));
        prog->scaled_dims = (int*) malloc(nvar*sizeof(int));
        prog->coloured_dims = 0;
        for (i=0;i<nvar;i++) {
            prog->total_coloured_stmts[i] = 0;
            prog->scaled_dims[i] = 0;
        }

        /* printf("[Pluto]: Num hyperplanes found so far %d\n", prog->num_hyperplanes); */
        find_permutable_dimensions_scc_based(colour, prog);

        IF_DEBUG(printf("[Pluto] Colouring Successful\n"););
        IF_DEBUG(pluto_print_colours(colour,prog););

        if(!options->silent && options->debug) {
            printf("[Pluto]: Transformations before skewing \n");
            pluto_transformations_pretty_print(prog);
        }

        introduce_skew(prog);

        free(colour);
        free(prog->total_coloured_stmts);
        free(prog->scaled_dims);
#endif
    } else {

        do{
            /* Number of linearly independent solutions remaining to be found
             * (maximum across all statements) */
            int num_sols_left;

            if (options->fuse == NO_FUSE)   {
                ddg_compute_scc(prog);
                cut_all_sccs(prog, ddg);
            }

            num_sols_left = 0;
            for (s=0; s<nstmts; s++) {
                /* Num linearly independent hyperplanes remaining to be
                 * found for a statement; take max across all */
                num_sols_left = PLMAX(num_sols_left, stmts[s]->dim_orig
                        - pluto_stmt_get_num_ind_hyps(stmts[s]));
            }
            /* Progress in the EAGER mode is made every time a solution is found;
             * thus, the maximum number of linearly independent solutions
             * remaining to be found is the difference between the number required
             * for the deepest statement and the number found so far for the
             * deepest statement (since in EAGER mode, if there was a statement
             * that had fewer than num_ind_sols_found linearly independent hyperplanes,
             * it means it didn't need that many hyperplanes and all of its
             * linearly independent solutions had been found */
            assert(hyp_search_mode == LAZY || num_sols_left == num_ind_sols_req - num_ind_sols_found);

            nsols = find_permutable_hyperplanes(prog, hyp_search_mode,
                    num_sols_left, depth);

            IF_DEBUG(fprintf(stdout, "[pluto] pluto_auto_transform: band level %d; %d hyperplane(s) found\n",
                        depth, nsols));
            IF_DEBUG2(pluto_transformations_pretty_print(prog));

            num_ind_sols_found = pluto_get_max_ind_hyps(prog);

            if (nsols >= 1) {
                /* Diamond tiling: done for the first band of permutable loops */
                if (options->lbtile && nsols >= 2 && !conc_start_found) {
                    conc_start_found = pluto_diamond_tile(prog);
                }

                for (j=0; j<nsols; j++)      {
                    /* Mark dependences satisfied by this solution */
                    dep_satisfaction_update(prog, stmts[0]->trans->nrows - nsols + j);
                    ddg_update(ddg, prog);
                }
            }else{
                /* Satisfy inter-scc dependences via distribution since we have 
                 * no more fusable loops */

                ddg_compute_scc(prog);

                if (get_num_unsatisfied_inter_scc_deps(prog) >= 1) {
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
                    /* Only one SCC or multiple SCCs with no unsatisfied inter-SCC
                     * deps, and no solutions found  */
                    if (hyp_search_mode == EAGER)   {
                        IF_DEBUG(printf("[pluto] Switching to LAZY mode\n"););
                        hyp_search_mode = LAZY;
                    }else if (!deps_satisfaction_check(prog)) {
                        assert(hyp_search_mode == LAZY);
                        /* There is a problem; solutions should have been found if
                         * there were no inter-scc deps, and some unsatisfied deps
                         * existed */
                        if (options->debug || options->moredebug) {
                            printf("\tNumber of unsatisfied deps: %d\n",
                                    get_num_unsatisfied_deps(prog->deps, prog->ndeps));
                            printf("\tNumber of unsatisfied inter-scc deps: %d\n",
                                    get_num_unsatisfied_inter_scc_deps(prog));
                            fprintf(stdout, "[pluto] WARNING: Unfortunately, pluto cannot find any more hyperplanes.\n");
                            fprintf(stdout, "\tThis is usually a result of (1) a bug in the dependence tester,\n");
                            fprintf(stdout, "\tor (2) a bug in Pluto's auto transformation,\n");
                            fprintf(stdout, "\tor (3) an inconsistent .fst/.precut in your working directory.\n");
                            fprintf(stdout, "\tTransformation found so far:\n");
                            pluto_transformations_pretty_print(prog);
                            pluto_compute_dep_directions(prog);
                            pluto_compute_dep_satisfaction(prog);
                            pluto_print_dep_directions(prog);
                        }
                        denormalize_domains(prog);
                        printf("[pluto] WARNING: working with original (identity) transformation (if they exist)\n");
                        /* Restore original ones */
                        for (i=0; i<nstmts; i++) {
                            stmts[i]->trans = orig_trans[i];
                            stmts[i]->hyp_types = orig_hyp_types[i];
                            prog->num_hyperplanes = orig_num_hyperplanes;
                            prog->hProps = orig_hProps;
                        }
                        return 1;
                    }
                }
            }
            /* Under LAZY mode, do a precise dep satisfaction check to take 
             * care of partial satisfaction (rarely needed) */
            if (hyp_search_mode == LAZY) pluto_compute_dep_satisfaction_precise(prog);
            depth++;
        }while (!pluto_transformations_full_ranked(prog) || 
                !deps_satisfaction_check(prog));
    }

 /* Deallocate the fusion conflict graph */
    if (options->dfp){
#ifdef GLPK
        ddg = prog->ddg;
        for(i=0; i<ddg->num_sccs; i++){
            free(ddg->sccs[i].vertices);
        }
        graph_free(prog->fcg);
#endif
    }
    if (options->lbtile && !conc_start_found) {
        PLUTO_MESSAGE(printf("[pluto] Diamond tiling not possible/useful\n"););
    }

    denormalize_domains(prog);

    for (i=0; i<nstmts; i++)    {
        pluto_matrix_free(orig_trans[i]);
        free(orig_hyp_types[i]);
    }
    free(orig_trans);
    free(orig_hyp_types);
    free(orig_hProps);

    IF_DEBUG(printf("[pluto] pluto_auto_transform: successful, done\n"););

    return 0;
}


int get_num_unsatisfied_deps(Dep **deps, int ndeps)
{
    int i, count;

    count = 0;
    for (i=0; i<ndeps; i++) {
        if (IS_RAR(deps[i]->type))   continue;
        if (!deps[i]->satisfied)  {
            IF_DEBUG(printf("\tUnsatisfied dep %d\n", i+1));
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

int get_num_unsatisfied_inter_scc_deps(PlutoProg *prog)
{
    int i;

    int count = 0;
    for (i=0; i<prog->ndeps; i++) {
        Dep *dep = prog->deps[i];
        if (IS_RAR(dep->type))   continue;
        Stmt *src_stmt = prog->stmts[dep->src];
        Stmt *dest_stmt = prog->stmts[dep->dest];
        if (src_stmt->scc_id != dest_stmt->scc_id && !dep->satisfied)  {
            count++;
        }
    }

    return count;
}


void ddg_print(Graph *g)
{
    pluto_matrix_print(stdout, g->adj);
}


/*
 * Update the DDG - should be called when some dependences
 * are satisfied
 **/
void ddg_update(Graph *g, PlutoProg *prog)
{
    int i, j;
    Dep *dep;

    IF_DEBUG(printf("[pluto] updating DDG\n"););

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

/* Compute the connected components of the graph */
void ddg_compute_cc(PlutoProg *prog)
{
    int i;
    int cc_id = -1;
    int num_cc = 0;
    int stmt_id;
    int time = 0;
    IF_DEBUG(printf("[pluto] ddg_compute_cc\n"););
    Graph *g = prog->ddg;
    /* Make the graph undirected. */
    Graph *gU = get_undirected_graph(g);
    for (i=0; i<gU->nVertices; i++){
        gU->vertices[i].vn = 0;
    }
    for (i=0; i<gU->nVertices; i++){
        if (gU->vertices[i].vn == 0){
            cc_id++;
            num_cc++;
            gU->vertices[i].cc_id = cc_id;
            dfs_vertex(gU,&gU->vertices[i],&time);
            gU->vertices[i].cc_id = cc_id;
        }
        g->vertices[i].cc_id = gU->vertices[i].cc_id;
        stmt_id = g->vertices[i].id;
        assert(stmt_id == i);
        prog->stmts[i]->cc_id = g->vertices[i].cc_id;
    }
    g->num_ccs = num_cc;
    graph_free(gU);
}

/* Compute the SCCs of a graph (usig Kosaraju's algorithm) */
void ddg_compute_scc(PlutoProg *prog)
{
    int i;

    IF_DEBUG(printf("[pluto] ddg_compute_scc\n"););

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
        g->sccs[i].sol = NULL;
        g->sccs[i].is_parallel = 0;
    }

    graph_free(gT);

    graph_print_sccs(g);
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

/* Are these statements completely fused until the innermost level */
int pluto_are_stmts_fused(Stmt **stmts, int nstmts, const PlutoProg *prog)
{
    int num;

    if (prog->num_hyperplanes <= 1) return 1;

    Ploop **loops = pluto_get_loops_under(stmts, nstmts, prog->num_hyperplanes-2, prog, &num);
    pluto_loops_free(loops, num);

    return num==1;
}


/* 
 * Diamond Tiling
 */
int pluto_diamond_tile(PlutoProg *prog)
{
    int b, d, nbands, conc_start_enabled, conc_start_enabled_band;

    IF_DEBUG(printf("[pluto] pluto_diamond_tile\n")); 

    conc_start_enabled = 0;

    /* Get the permutability constraints since a call to
     * detect_transformation_properties with update dep satisfaction levels
     * and we won't get the constraints we want */

    /* Don't free basecst */
    PlutoConstraints *basecst = get_permutability_constraints(prog);

    pluto_compute_dep_directions(prog);
    pluto_compute_dep_satisfaction(prog);

    Band **bands = pluto_get_outermost_permutable_bands(prog, &nbands);

    for (b=0; b<nbands; b++) {
        PlutoMatrix **cone_complement_hyps;
        Band *band = bands[b];
        int evict_pos;
        int i, first_loop_hyp, cone_complement_pos, ni, s;

        /* Band should not have outer parallelism */
        if (pluto_loop_is_parallel(prog, band->loop)) continue;

        /* Band should have inner parallelism */
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

        /* Domains should allows point-wise concurrent start */
        PlutoMatrix *conc_start_faces = get_face_with_concurrent_start(prog, band);
        if (!conc_start_faces) continue;

        /* face with concurrent start shouldn't be normal to all hyperplanes
         * of all statements in this band */
        for (s=0; s<band->loop->nstmts; s++) {
            for (d=band->loop->depth; d<band->loop->depth + band->width; d++) {
                if (!pluto_vector_is_normal(band->loop->stmts[s]->trans, d,
                        conc_start_faces, s)) break;
            }
            if (d < band->loop->depth + band->width) break;
        }
        if (s == band->loop->nstmts) {
            printf("row normal\n");
            continue;
        }

        cone_complement_hyps = malloc(
                band->loop->nstmts*sizeof(PlutoMatrix *));
        for (i=0; i<band->loop->nstmts; i++) {
            cone_complement_hyps[i] = NULL;
        }

        first_loop_hyp = band->loop->depth;
        /* 
         * Find hyperplane that will be replaced by the newly found
         * hyperplane
         * Concurrent start pertains to the first band alone
         */
        evict_pos = find_hyperplane_to_be_evicted(band, conc_start_faces, prog);

        /* If we haven't yet found the cone_complement_pos, just 
         * choose the first one as the cone_complement_pos */
        cone_complement_pos = first_loop_hyp;

        /* If first_loop_hyp hyperplane itself is to be replaced, 
         * choose the next one as cone_complement_pos */
        if (evict_pos == first_loop_hyp) cone_complement_pos++ ;
        conc_start_enabled_band = find_cone_complement_hyperplane(band, conc_start_faces, evict_pos,
                cone_complement_pos, basecst, prog, cone_complement_hyps);
        pluto_matrix_free(conc_start_faces);

        /* Re-arrange the transformation matrix if concurrent start
         * was found, store the replaced hyperplane so that it can be 
         * put back for the right intra-tile order */
        if (conc_start_enabled_band) {
            IF_DEBUG(printf("[pluto] Transformations before concurrent start enable\n")); 
            IF_DEBUG(pluto_transformations_pretty_print(prog););
            IF_DEBUG(pluto_print_hyperplane_properties(prog););
            for (i=0; i<band->loop->nstmts; i++){
                Stmt *stmt = band->loop->stmts[i];
                /* Since we do concurrent start only once */
                assert(stmt->evicted_hyp == NULL);
                stmt->evicted_hyp = pluto_matrix_alloc(1, stmt->trans->ncols);
                copy_hyperplane(stmt->evicted_hyp->val[0], 
                        stmt->trans->val[evict_pos], stmt->trans->ncols);
                copy_hyperplane(stmt->trans->val[evict_pos], 
                        cone_complement_hyps[i]->val[0], stmt->trans->ncols);
                stmt->evicted_hyp_pos = evict_pos;
            }
        }

        conc_start_enabled |= conc_start_enabled_band;

        for (i=0; i<band->loop->nstmts; i++) {
            pluto_matrix_free(cone_complement_hyps[i]);
        }
        free(cone_complement_hyps);
    }

    pluto_bands_free(bands, nbands);

    if (conc_start_enabled) {
        pluto_dep_satisfaction_reset(prog);
        PLUTO_MESSAGE(printf("[pluto] Concurrent start hyperplanes found\n"););
        IF_DEBUG(printf("[pluto] Transformations after concurrent start enable\n")); 
        IF_DEBUG(pluto_transformations_pretty_print(prog););
    }

    return conc_start_enabled;
}
