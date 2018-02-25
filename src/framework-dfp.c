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

#include <sys/time.h>

#include "math_support.h"
#include "constraints.h"
#include "pluto.h"
#include "program.h"

#include <isl/constraint.h>
#include <isl/mat.h>
#include <isl/set.h>
#include "candl/candl.h"

static double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

/* Checks for feasibility of constraints.
 * If feasible then return the solution else returns NULL */
/* double* pluto_fusion_constraints_feasibility_solve_glpk(PlutoConstraints *cst, PlutoProg *prog, int src_dim, int target_dim, int fcg_src, int fcg_dest) */
double* pluto_fusion_constraints_feasibility_solve_glpk(PlutoConstraints *cst, PlutoMatrix *obj)
{
    /* int i, j, nstmts, nvar, npar, stmt_offset; */
    /* PlutoMatrix *obj; */
    /* Stmt **stmts; */
    /* glp_prob *lp; */
    double* sol;
    /* Create the data dependence graph */
    /* prog->ddg = ddg_create(prog); */
    /* ddg_compute_scc(prog); */
    /* ddg_compute_cc(prog); */

    /* nstmts = prog->nstmts; */
    /* nvar = prog->nvar; */
    /* npar = prog->npar; */
    /* stmts = prog->stmts; */


    /* obj = construct_cplex_objective(cst,prog); */
    sol = pluto_fcg_constraints_lexmin_glpk(cst, obj);
    return sol;

    /* lp = create_glpk_problem_from_pluto_constraints(cst); */

    /* for (j=0; j<obj->ncols; j++) { */
    /*     glp_set_obj_coef(lp, j+1, obj->val[0][j]); */
    /* } */
    /* pluto_matrix_free(obj); */

    /* Set the lower bounds of all variables to zero in the LP problem. 
     * The eqality constraints for loop and constant part are set in the constraint matrix. */
    /* for (i=0; i<glp_get_num_cols(lp); i++) { */
    /*     glp_set_col_bnds(lp, i+1, GLP_LO, 0.0,0.0); */
    /*     if (options->ilp) { */
    /*         glp_set_col_kind(lp,i+1, GLP_IV); */
    /*     } */
    /* } */

/* Debug code. To be cleaned up after testing */
    /* for (i=0; i<nstmts; i++) { */
    /*     stmt_offset = npar+1+(nvar+1)*i; */
    /*     for (j=0; j<stmts[i]->dim_orig; j++){ */
    /*         if(((stmt_offset+j)!=src_dim) && ((stmt_offset+j)!=target_dim)){ */
    /*             glp_set_col_bnds(lp, stmt_offset+j+1, GLP_FX, 0.0,0.0); */
    /*         } */
    /*         else */
    /*             glp_set_col_bnds(lp, stmt_offset+j+1, GLP_FR, 0.0, 0.0); */
    /*         if (options->ilp) { */
    /*             glp_set_col_kind(lp,stmt_offset+j+1, GLP_IV); */
    /*         } */
    /*     } */
    /*     #<{(| Set bound for the constant |)}># */
    /*     glp_set_col_bnds(lp,stmt_offset+nvar+1, GLP_FR, 0.0,0.0); */
    /*     if (options->ilp) { */
    /*         glp_set_col_kind(lp,stmt_offset+nvar+1, GLP_IV); */
    /*     } */
    /* } */

    /* if (options->debug) { */
    /*     glp_write_lp(lp,NULL,"pairwise-debug.lp"); */
    /* } */
    /*  */
    /* if (!options->debug && !options->moredebug) { */
    /*     glp_term_out(GLP_OFF); */
    /* } */
    /*  */
    /*  */
    /* glp_smcp parm; */
    /* glp_init_smcp(&parm); */
    /* parm.presolve = GLP_ON; */
    /*  */
    /* parm.msg_lev = GLP_MSG_OFF; */
    /* IF_DEBUG(parm.msg_lev = GLP_MSG_ON;); */
    /* IF_MORE_DEBUG(parm.msg_lev = GLP_MSG_ALL;); */
    /*  */
    /* glp_scale_prob(lp, GLP_SF_AUTO); */
    /* glp_adv_basis(lp, 0); */
    /* glp_simplex(lp, &parm); */

    /* int lp_status = glp_get_status(lp); */

    /* The cost matrix isnt currently used anywhere. It was used for experimental purposes */
    /* if (lp_status == GLP_INFEAS || lp_status == GLP_UNDEF) { */
    /*     glp_delete_prob(lp); */
    /*     return NULL; */
    /* } else { */
    /*     glp_iocp iocp; */
    /*     glp_init_iocp(&iocp); */
        /* The default is 1e-5; one may need to reduce it even further
         * depending on how large a coefficient we might see */
    /*     iocp.tol_int = 1e-7; */
    /*     IF_DEBUG(printf("Setting GLPK integer tolerance to %e\n", iocp.tol_int)); */
    /*  */
    /*  */
    /*     iocp.msg_lev = GLP_MSG_OFF; */
    /*     IF_DEBUG(iocp.msg_lev = GLP_MSG_ON;); */
    /*     IF_MORE_DEBUG(iocp.msg_lev = GLP_MSG_ALL;); */
    /*  */
    /*     glp_intopt(lp, &iocp); */
    /*  */
    /*     double* sol; */
    /*     sol = (double*) malloc((cst->ncols-1)*sizeof(double)); */
    /*     for(i=0; i<glp_get_num_cols(lp); i++){ */
    /*         sol[i] = glp_mip_col_val(lp, i+1); */
    /*     } */
    /*  */
    /*     glp_delete_prob(lp); */
    /*     return sol; */
    /* } */
}


/* Adds edges in FCG corresponding to the satements represented by the nodes v1 and v2 in DDG*/
void fcg_add_pairwise_edges(Graph *fcg, int v1, int v2, PlutoProg *prog, int *colour, PlutoConstraints *boundcst, int current_colour, PlutoConstraints **conflicts, PlutoMatrix *obj)
{
    Graph *ddg;
    int i,j,ndeps,nstmts,nvar,npar, src_offset,dest_offset,fcg_offset1,fcg_offset2;
    double *sol;
    int row_offset;
    Stmt **stmts;
    PlutoConstraints *conflictcst;

    ddg = prog->ddg;
    Dep **deps, *dep;
    ndeps = prog->ndeps;
    deps = prog->deps;

    stmts = prog->stmts;
    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;



    double tstart = rtclock();
    assert (*conflicts != NULL);
    conflictcst = *conflicts;

    prog->fcg_cst_alloc_time += rtclock() - tstart;
    row_offset = conflictcst->nrows-CST_WIDTH-1;
    /* conflictcst->ncols = CST_WIDTH; */

    /* conflictcst->val[row_offset][CST_WIDTH-1] = -1; */
    /* conflictcst->val[row_offset +1][CST_WIDTH-1] = -1; */

    tstart = rtclock();
    for (i=0; i<ndeps; i++){
        dep = deps[i];
        /*if (options->varliberalize && dep->skipdep) {
          continue;
          }*/
        if(dep_is_satisfied(dep)){
            continue;
        }
        if ((dep->src == v1 && dep->dest == v2)||(dep->src==v2 && dep->dest ==v1)){
            if(dep->cst == NULL){
                compute_pairwise_permutability(dep,prog);
            }
            IF_DEBUG(printf("Adding Constraints for dependence %d\n",i););
            pluto_constraints_add(conflictcst, dep->cst);
        }
    }
    if(stmts[v1]->intra_stmt_dep_cst != NULL){
        pluto_constraints_add(conflictcst,stmts[v1]->intra_stmt_dep_cst);
    }
    if(stmts[v2]->intra_stmt_dep_cst != NULL){
        pluto_constraints_add(conflictcst,stmts[v2]->intra_stmt_dep_cst);
    }


    src_offset = npar+1+(nvar+1)*v1;
    dest_offset = npar+1+(nvar+1)*v2;

    fcg_offset1 = ddg->vertices[v1].fcg_stmt_offset;
    fcg_offset2 = ddg->vertices[v2].fcg_stmt_offset;

    /* Solve Pluto LP by setting corresponding coeffs to 0 without any objective.
     * This is the check for fusability of two dimensions */
    for(i=0; i<stmts[v1]->dim_orig; i++){
        /* note that the vertex should not be coloured. Even if the vertex has a 
         * self edge, it must be considered during construction of the FCG. This
         * is because,even after satisfying the permute preventing dep, it might
         * still prevent fusion. */
        if (colour[fcg_offset1 + i] == 0 || colour[fcg_offset1 + i] == current_colour) {

            /* Set the lower bound of i^th dimension of v1 to 1 */
            conflictcst->val[row_offset + src_offset+i][CST_WIDTH-1] = -1;
            conflictcst->is_eq[row_offset + src_offset+i] = 0;

            for (j=0; j<stmts[v2]->dim_orig; j++) {
                if (colour[fcg_offset2 + j] == 0 || colour[fcg_offset1 + i] == current_colour) {

                    /* Set the lower bound of i^th dimension of v1 to 1 */
                    conflictcst->val[row_offset + dest_offset+j][CST_WIDTH-1] = -1;
                    conflictcst->is_eq[row_offset + dest_offset+j] = 0;

                    /* Check if fusing ith dimesion of the source with ith dimension
                     * of the target is valid */

                    /* conflictcst->val[row_offset + 1][dest_offset+j] = 1; */

                    /* printf("[pluto]: Fusion conflict constraints:\n"); */
                    /* pluto_constraints_cplex_print(stdout,conflictcst); */
                    /* Solve these constraints using GLPK */
                    prog->num_lp_calls ++;
                    tstart = rtclock();
                    sol = pluto_fusion_constraints_feasibility_solve_glpk(conflictcst, obj);
                    prog->cst_solve_time += rtclock() - tstart;

                    /* If no solutions, then dimensions are not fusable. Add an edge in the conflict graph. */
                    if(sol == NULL)
                    {
                        IF_DEBUG(printf("Unable to fuse Dimesnion %d of statement %d with dimension %d of statement %d \n",i,v1,j,v2););
                        IF_DEBUG(printf(" Adding edge %d to %d in fcg\n",fcg_offset1+i,fcg_offset2+j););
                        fcg->adj->val[fcg_offset1+i][fcg_offset2+j] = 1;
                    } else {
                        free(sol);
                    }
                    /* Unset the lowerbound for the coefficient of c_i. The same constraint matrix is reused for all coeffs. */
                    conflictcst->val[row_offset+dest_offset+j][CST_WIDTH-1] = 0;
                    conflictcst->is_eq[row_offset+dest_offset+j] = 1;
                    /* conflictcst->val[row_offset+1][dest_offset+j] = 0; */
                }
            }

            /* Unset the lowerbound for the coefficient of c_i. The same constraint matrix is reused for all coeffs. */
            conflictcst->val[row_offset + src_offset+i][CST_WIDTH-1] = 0;
            conflictcst->is_eq[row_offset + src_offset+i] = 1;
            /* conflictcst->val[row_offset+0][src_offset+i] = 0; */
        }
    }
    /* conflictcst->nrows = row_offset+2; */
    /* conflictcst->ncols = CST_WIDTH; */
    /* conflictcst->val[row_offset][CST_WIDTH-1] = 0; */
    /* conflictcst->val[row_offset+1][CST_WIDTH-1] = 0; */
    return;
}


/* Computes intra statement dependence constraints for every unstisfied dependence */
void compute_intra_stmt_deps(PlutoProg *prog)
{
    int ndeps, src_stmt,dest_stmt,i;
    Dep **deps;
    Dep *dep;
    Stmt **stmts;
    Stmt *stmt;

    deps = prog->deps;
    ndeps = prog->ndeps;
    stmts = prog->stmts;
    for (i=0; i<ndeps; i++) {
        dep = deps[i];
        if (options->rar ==0 && IS_RAR(dep->type)) {
            continue;
        }

        /* if (options->varliberalize && dep->skipdep) {
            continue;
        } */

        if (dep_is_satisfied(dep)) {
            continue;
        }

        src_stmt = dep->src;
        dest_stmt = dep->dest;
        if (src_stmt == dest_stmt) {
            stmt = stmts[src_stmt];
            IF_DEBUG(printf("Computing Intra statement deps for statement: %d Dep: %d\n",src_stmt,dep->id););
            if (dep->cst==NULL) {
                compute_pairwise_permutability(dep,prog);
            }
            if (stmt->intra_stmt_dep_cst == NULL) {
                stmt->intra_stmt_dep_cst = pluto_constraints_alloc(dep->cst->nrows,dep->cst->ncols);
                stmt->intra_stmt_dep_cst->nrows = dep->cst->nrows;
                stmt->intra_stmt_dep_cst -> ncols = dep->cst->ncols;
                pluto_constraints_copy(stmt->intra_stmt_dep_cst, dep->cst);
            } else {
                pluto_constraints_add(stmt->intra_stmt_dep_cst, dep->cst);
            }
        }
    }

}

/* Adds permute preventing edges for intra statement dependences.  These edges
 * are added as self loops on FCG vertices.  These vertices can not be coloured.
 * Inter statement permute preventing deps do not cause a problem as they will be
 * represented by the inter statement edges */
#if 0
void add_permute_preventing_edges(Graph* fcg, int *colour, PlutoProg *prog, PlutoConstraints* boundcst, int current_colour)
{
    int nstmts,nvar,npar,i,j,stmt_offset,fcg_stmt_offset;
    double *sol;
    Stmt **stmts;
    PlutoConstraints *intra_stmt_dep_cst, *coeff_bounds;

    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;

    stmts = prog->stmts;

    /* Compute the intra statment dependence constraints */
    compute_intra_stmt_deps(prog);

    fcg_stmt_offset = 0;
    for (i=0; i<nstmts; i++) {
        if (stmts[i]->intra_stmt_dep_cst!=NULL) {
            /* Constraints to check permutability are added in the first row */
            double tstart = rtclock();
            coeff_bounds = pluto_constraints_alloc(1,CST_WIDTH);
            prog->fcg_cst_alloc_time += rtclock() - tstart;
            coeff_bounds->nrows = 1;
            coeff_bounds->ncols = CST_WIDTH;
            /* Add the intra statement dependence constraints and bounding constraints */
            intra_stmt_dep_cst = stmts[i]->intra_stmt_dep_cst;

            pluto_constraints_add(coeff_bounds,intra_stmt_dep_cst);
            pluto_constraints_add(coeff_bounds,boundcst);

            stmt_offset = (npar+1)+ i*(nvar+1);
            coeff_bounds->val[0][CST_WIDTH-1] = -1;

            for (j=0; j<stmts[i]->dim_orig; j++) {
                if (colour[fcg_stmt_offset + j]==0 || colour[fcg_stmt_offset+j] == current_colour) {
                    IF_DEBUG(printf("[Permute_preventing_edges]: Checking permutability of dimension %d of statement %d \n",j,i););
                    coeff_bounds->val[0][stmt_offset+j] = 1;
                    prog->num_lp_calls++;

                    sol = pluto_fusion_constraints_feasibility_solve_glpk(coeff_bounds,prog,stmt_offset+j, stmt_offset+j, fcg_stmt_offset+j, fcg_stmt_offset+j);
                    /* If the constraints are infeasible then add a self edge in the FCG */

                    if (sol == NULL) {
                        IF_DEBUG(printf("Dimension %d of Statment %d is not permutable\n",j,i););
                        fcg->adj->val[fcg_stmt_offset+j][fcg_stmt_offset+j] = 1;
                    } else {
                        free(sol);
                    }
                    /* reset the coeff bound of this dimension */
                    coeff_bounds->val[0][stmt_offset+j] = 0;
                }
            }
            tstart = rtclock();
            pluto_constraints_free(coeff_bounds);
            prog->fcg_cst_alloc_time += rtclock() - tstart;
        }
        fcg_stmt_offset += stmts[i]->dim_orig;
    }
}
#endif

/* Removes all the edges in the FCG from a dimension of a statement that is in an 
 * SCC whose ID is less than or equal to scc1 to the dimension of a statement 
 * present in a SCC greater than or equal to scc2. */
void update_fcg_between_sccs(Graph *fcg, int scc1, int scc2, PlutoProg *prog)
{
    int nstmts, i, j, k, l, stmt_offset1, stmt_offset2;
    int nvar, npar;
    Graph *ddg;
    Stmt **stmts;
    double tstart;

    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;
    ddg = prog->ddg;
    stmts = prog->stmts;

    assert (fcg->to_be_rebuilt == false);
    tstart = rtclock();

    if (nstmts == 1) {
        return;
    }
    /* Assumes that the DDG has already been cut. */
    if (options->fuse == NO_FUSE) {
        for (i=1; i<nstmts; i++) {
            for (j=0; j<i ; j++) {
                if (stmts[i]->trans->val[stmts[i]->trans->nrows-1][nvar + npar] !=
                        stmts[j]->trans->val[stmts[j]->trans->nrows-1][nvar + npar]) {
                    stmt_offset1 = ddg->vertices[i].fcg_stmt_offset;
                    stmt_offset2 = ddg->vertices[i].fcg_stmt_offset;
                    for(k=0; k<stmts[i]->dim_orig; k++) {
                        for(l=0; l<stmts[j]->dim_orig; l++) {
                            fcg->adj->val[stmt_offset1+k][stmt_offset2+l] = 0;
                            fcg->adj->val[stmt_offset2+l][stmt_offset1+k] = 0;
                        }
                    }
                }
            }
        }
    } else {
        IF_DEBUG(printf("Updating FCG between SCCs%d and %d\n",scc1,scc2););
        for (i=0; i<nstmts; i++) {
            for(j=0; j<nstmts; j++) {
                if ((stmts[i]->scc_id >= scc2 && stmts[j]->scc_id<scc2) ||
                        (stmts[j]->scc_id>=scc2 && stmts[i]->scc_id < scc2)) {
                    stmt_offset1 = ddg->vertices[i].fcg_stmt_offset;
                    stmt_offset2 = ddg->vertices[j].fcg_stmt_offset;
                    for(k=0; k<stmts[i]->dim_orig; k++) {
                        for(l=0; l<stmts[j]->dim_orig; l++) {
                            fcg->adj->val[stmt_offset1+k][stmt_offset2+l] = 0;
                            fcg->adj->val[stmt_offset2+l][stmt_offset1+k] = 0;
                        }
                    }
                }
            }
        }
    }

    prog->fcg_update_time += rtclock() - tstart;
}


/* Build the fusion conflict graph for a given program.  The current colour is 
 * used to rebuild FCG for the current level.  This is need in case we are 
 * separating out construction of FCG for permute preventing dependence and
 * fusion preventing dependences */
Graph* build_fusion_conflict_graph(PlutoProg *prog, int *colour, int num_nodes, int current_colour)
{
    int i,j,k,stmt_offset,nstmts, nvar,npar, nrows;
    Stmt **stmts;
    Graph *ddg;
    Graph *fcg;
    double t_start, t_start2;
    PlutoConstraints *boundcst, **conflicts;
    PlutoMatrix *obj;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;

    ddg = prog->ddg;

    t_start = rtclock();

    fcg = graph_alloc(num_nodes);

    boundcst = get_coeff_bounding_constraints(prog);
    /* Add premutation preventing intra statement dependence edges in the FCG.
     * These are self loops on vertices of the FCG. */ 

    t_start2 = rtclock();
    /* add_permute_preventing_edges(fcg, colour, prog, boundcst, current_colour); */
    IF_DEBUG(printf("[Pluto] Build Fusion Conflict graph: FCG add permute preventing edges: %0.6lfs\n",rtclock()-t_start2););
    IF_DEBUG(printf("[Pluto] Build Fusion Conflict graph: Number of LP calls to check dimension permutability: %ld\n",prog->num_lp_calls););

    /* Add inter statement fusion and permute preventing edges.  */
    t_start2 = rtclock();
    conflicts = (PlutoConstraints**)malloc(sizeof(PlutoConstraints*));

    /* The last CST_WIDTH-1 number of rows represent the bounds on the coeffcients  */
    *conflicts = pluto_constraints_alloc(CST_WIDTH-1 + boundcst->nrows,CST_WIDTH);
    (*conflicts)->ncols = CST_WIDTH;

    
    obj = construct_cplex_objective(*conflicts, prog);

    pluto_constraints_add(*conflicts, boundcst);
    assert((*conflicts)->nrows == boundcst->nrows);

    nrows = boundcst->nrows;
    (*conflicts)->nrows = boundcst->nrows + CST_WIDTH-1;

    /* u and w are lower bounded by 0 */
    for (i=0; i<npar+1; i++) {
        (*conflicts)-> val[nrows+i][i] = 1;
    }

    /* In the last CST_WIDTH-npar+1 number of rows, correspond to equality constraints.
     * These are changed during dimension wise computation of edges of the FCG*/
    for (i=npar+1; i<CST_WIDTH-1; i++) {
        (*conflicts)->is_eq[nrows+i] = 0;
        (*conflicts)->val[nrows+i][i] = 1;
    }
    
    pluto_constraints_cplex_print(stdout, *conflicts);

    for (i=0; i<nstmts-1; i++) {
        /* The lower bound for  constant shift of i^th statement is 0 */
        (*conflicts)->is_eq[npar+1+i*(nvar+1)+nvar] = 0;
        for (j=i+1; j<nstmts; j++) {
            if (is_adjecent(ddg,i,j)) {
                /* Set the lower bound of the constant shift to be 1. */
                (*conflicts)->is_eq[npar+1+j*(nvar+1)+nvar] = 0;
                fcg_add_pairwise_edges(fcg,i,j,prog, colour, boundcst, current_colour, conflicts, obj);
                (*conflicts)->is_eq[npar+1+j*(nvar+1)+nvar] = 1;
            }
        }
        (*conflicts)->is_eq[npar+1+i*(nvar+1)+nvar] = 1;
    }
    IF_DEBUG(printf("[Pluto] Build Fusion Conflict graph: FCG add parwise edges: %0.6lfs\n", rtclock()-t_start2););

    /* Add egdes between different dimensions of the same statement */
    stmt_offset=0;
    for (i=0; i<nstmts;i++) {
        for (j=stmt_offset; j<stmt_offset+stmts[i]->dim_orig; j++) {
            fcg->vertices[j].fcg_stmt_offset = i;
            for (k=j+1; k<stmt_offset+stmts[i]->dim_orig;k++) {
                fcg->adj->val[j][k] = 1;
                fcg->adj->val[k][j] = 1;
            }
        }
        stmt_offset += stmts[i]->dim_orig;


        /* Remove the intra statement dependence constraints. Else the permutability constraints 
         * might be incorrect for rebuilding the fusion conflict graph.  */

        pluto_constraints_free(stmts[i]->intra_stmt_dep_cst);
        stmts[i]->intra_stmt_dep_cst = NULL;
    }

    pluto_constraints_free(boundcst);
    pluto_constraints_free(*conflicts);
    free(conflicts);
    prog->fcg_const_time += rtclock() - t_start;

    IF_DEBUG(printf("[Pluto] Build FCG: Total number of LP calls in building the FCG: %ld\n",prog->num_lp_calls););
    return fcg;
}
