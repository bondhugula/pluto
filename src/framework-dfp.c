/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * Author: Aravind Acharya
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

int scale_shift_permutations(PlutoProg *prog, int *colour, int c);

static double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

/*********************************** FCG construction routines *****************************************/

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
    row_offset = conflictcst->nrows-CST_WIDTH+1;
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
void add_permute_preventing_edges(Graph* fcg, int *colour, PlutoProg *prog, PlutoConstraints* boundcst, int current_colour, PlutoMatrix *obj)
{
    int nstmts,nvar,npar,i,j,stmt_offset,fcg_stmt_offset;
    int nrows;
    double *sol;
    Stmt **stmts;
    PlutoConstraints *intra_stmt_dep_cst, *coeff_bounds;

    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;

    stmts = prog->stmts;

    nrows = boundcst->nrows-CST_WIDTH+1;


    /* Compute the intra statment dependence constraints */
    compute_intra_stmt_deps(prog);

    fcg_stmt_offset = 0;
    for (i=0; i<nstmts; i++) {
        if (stmts[i]->intra_stmt_dep_cst!=NULL) {
            /* Constraints to check permutability are added in the first row */
            double tstart = rtclock();
            coeff_bounds = pluto_constraints_alloc(1,CST_WIDTH);
            prog->fcg_cst_alloc_time += rtclock() - tstart;
            coeff_bounds->nrows = 0;
            coeff_bounds->ncols = CST_WIDTH;
            /* Add the intra statement dependence constraints and bounding constraints */
            intra_stmt_dep_cst = stmts[i]->intra_stmt_dep_cst;

            pluto_constraints_add(coeff_bounds,boundcst);

            pluto_constraints_add(coeff_bounds,intra_stmt_dep_cst);

            stmt_offset = (npar+1)+ i*(nvar+1);
            /* coeff_bounds->val[0][CST_WIDTH-1] = -1; */

            for (j=0; j<stmts[i]->dim_orig; j++) {
                if (colour[fcg_stmt_offset + j]==0 || colour[fcg_stmt_offset+j] == current_colour) {
                    IF_DEBUG(printf("[Permute_preventing_edges]: Checking permutability of dimension %d of statement %d \n",j,i););
                    coeff_bounds->is_eq[nrows+stmt_offset+j] = 0;
                    coeff_bounds->val[nrows+stmt_offset+j][CST_WIDTH-1] = -1;
                    /* coeff_bounds->val[0][stmt_offset+j] = 1; */
                    prog->num_lp_calls++;

                    sol = pluto_fusion_constraints_feasibility_solve_glpk(coeff_bounds,obj);
                    /* If the constraints are infeasible then add a self edge in the FCG */

                    if (sol == NULL) {
                        IF_DEBUG(printf("Dimension %d of Statment %d is not permutable\n",j,i););
                        fcg->adj->val[fcg_stmt_offset+j][fcg_stmt_offset+j] = 1;
                    } else {
                        free(sol);
                    }
                    /* reset the coeff bound of this dimension */
                    coeff_bounds->is_eq[nrows+stmt_offset+j] = 1;
                    coeff_bounds->val[nrows+stmt_offset+j][CST_WIDTH-1] = 0;
                    /* coeff_bounds->val[0][stmt_offset+j] = 0; */
                }
            }
            tstart = rtclock();
            pluto_constraints_free(coeff_bounds);
            prog->fcg_cst_alloc_time += rtclock() - tstart;
        }
        fcg_stmt_offset += stmts[i]->dim_orig;
    }
}

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
        (*conflicts)->is_eq[nrows+i] = 1;
        (*conflicts)->val[nrows+i][i] = 1;
    }
    
    pluto_constraints_cplex_print(stdout, *conflicts);

    /* Add premutation preventing intra statement dependence edges in the FCG.
     * These are self loops on vertices of the FCG. */ 
    add_permute_preventing_edges(fcg, colour, prog, *conflicts, current_colour, obj);

    IF_DEBUG(printf("[Pluto] Build Fusion Conflict graph: FCG add permute preventing edges: %0.6lfs\n",rtclock()-t_start2););
    IF_DEBUG(printf("[Pluto] Build Fusion Conflict graph: Number of LP calls to check dimension permutability: %ld\n",prog->num_lp_calls););

    /* Add inter statement fusion and permute preventing edges.  */
    t_start2 = rtclock();

    for (i=0; i<nstmts-1; i++) {
        /* The lower bound for  constant shift of i^th statement is 0 */
        (*conflicts)->is_eq[nrows + npar+1+i*(nvar+1)+nvar] = 0;
        for (j=i+1; j<nstmts; j++) {
            if (is_adjecent(ddg,i,j)) {
                /* Set the lower bound of the constant shift to be 1. */
                (*conflicts)->is_eq[nrows + npar+1+j*(nvar+1)+nvar] = 0;
                fcg_add_pairwise_edges(fcg,i,j,prog, colour, boundcst, current_colour, conflicts, obj);
                (*conflicts)->is_eq[nrows + npar+1+j*(nvar+1)+nvar] = 1;
            }
        }
        (*conflicts)->is_eq[nrows + npar+1+i*(nvar+1)+nvar] = 1;
    }
    IF_DEBUG(printf("[Pluto] Build Fusion Conflict graph: FCG add parwise edges: %0.6lfs\n", rtclock()-t_start2););

    pluto_matrix_free(obj);

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




/******************  FCG Colouring Routines **********************************/

/* Prints colour of each vertex of the FCG */
void pluto_print_colours(int *colour,PlutoProg *prog)
{
    int nstmts,i,j,stmt_offset;
    Stmt **stmts;

    nstmts = prog->nstmts;
    stmts = prog->stmts;

    stmt_offset = 0;


    for (i=0; i<nstmts;i++){
        for (j=0;j<stmts[i]->dim_orig;j++){
            printf("Colour of Dimension %d of Stmt %d: %d\n",j,i,colour[stmt_offset+j]);
        }
        stmt_offset+=j;
    }
}

/* Check if it is valid to give colour c to a vertex v in the fcg.
 * Colour is the array containing the colours assigned to each vertex */
bool is_valid_colour(int v, int c, Graph *fcg, int * colour)
{
    int i, fcg_nVertices;
    fcg_nVertices = fcg->nVertices;
    for (i=0;i< fcg_nVertices;i++){
        if((fcg->adj->val[i][v]==1||fcg->adj->val[v][i]==1) && colour[i]==c){
            return false;
        }
    }
    return true;
}

bool is_discarded(int v, int list[], int num)
{
    int i;
    for (i=0; i<num; i++){
        if(list[i]==v)
            return true;
    }
    return false;
}


/* Routine that returns the next vertex to be coloured. Currently returns the next vertex in the ordered list. */
int get_next_min_vertex(int fcg_stmt_offset, int stmt_id, int *list, int num, int pv, PlutoProg *prog){
    int i, min, npar, nvar, stmt_offset;
    Stmt **stmts;
    int scc_id;
    double *sol;

    nvar = prog->nvar;
    npar = prog->npar;
    stmts = prog->stmts;
    min = 0;

    /* if(pv == -1){ */
    /*     return min; */
    /* } */
    for (i=0; i<stmts[stmt_id]->dim_orig; i++) {
        if(!is_discarded(fcg_stmt_offset+i, list, num)) {
            if (options->lpcolour) {
                scc_id = stmts[stmt_id]->scc_id;
                sol = prog->ddg->sccs[scc_id].sol;
                assert (sol != NULL);
                stmt_offset = npar+1+(nvar+1)*stmt_id+i;
                if (sol[stmt_offset] == 0.0f) {
                    continue;
                }
            }
            min = i;
            break;
        }
        /* if(!is_discarded(fcg_stmt_offset + i, list, num)){ */
        /*     if(pv == -1){ */
        /*         return i; */
        /*     } */
        /*    #<{(| Set the first value of min_cost |)}># */
        /*     if(min_cost < 0){ */
        /*         min_cost = prog->costMatrix->val[fcg_stmt_offset+i][pv]; */
        /*         min = i; */
        /*     } */
        /*  */
        /*     if(prog->costMatrix->val[fcg_stmt_offset+i][pv]< min_cost){ */
        /*         min = i; */
        /*         min_cost = prog->costMatrix->val[fcg_stmt_offset+i][pv]; */
        /*     } */
        /* }  */
    }
    return min;
}



/* Colours the input SCC recursively.  The statement pos refers to the position 
 * of the statement in the list of vertices in the scc and pv refers to the 
 * previous vertex.  Returns true if the colouring is successful; 
 * else returns false.  */

bool colour_scc(int scc_id, int *colour, int c, int stmt_pos, int pv, PlutoProg *prog)
{
    int j, v, fcg_offset, stmt_id, nvar;
    Graph *ddg,*fcg;
    Scc *sccs;

    nvar = prog->nvar;
    ddg = prog->ddg;
    fcg = prog->fcg;
    sccs = ddg->sccs;

    int list[nvar];
    int num_discarded = 0;
    /* memset(list, -1, nvar); */

    /* ToDo: Check if this condition can really happen.  */
    if(stmt_pos >= sccs[scc_id].size){
        return true;
    }

    if(prog->coloured_dims >= sccs[scc_id].max_dim) {
        if(prog->coloured_dims > sccs[scc_id].max_dim){
            return true;
        }
        IF_DEBUG(printf("[colour SCC]: All Dimensions of statment %d in SCC %d have been coloured\n", sccs[scc_id].vertices[stmt_pos], scc_id););
        /* cut if the scc's are not already distributed and you can not colour further. */
        /* for each SCC which is greater than the current scc, if there is a dep edge 
         * between these scc's then cut between these scc's. The cut has to respect 
         * the existing dependence */

        /* This is just for experimental purposes. The following if code can be removed if the assert never fails */
        if(sccs[scc_id].size !=1 ){
            printf("SCC %d has size %d\n", scc_id,sccs[scc_id].size);
            int i;
            for(i=0;i<sccs[scc_id].size; i++){
                printf("S%d,,",sccs[scc_id].vertices[i]);
            }
            printf("\n");
        }
        assert (sccs[scc_id].size ==1);

        if(sccs[scc_id].size == 1){
            for(j=0;j<ddg->num_sccs; j++){
                if(scc_id!=j){
                    if((j < scc_id) && ddg_sccs_direct_connected(ddg,prog,j,scc_id)) {
                        IF_DEBUG(printf("[colour SCC]: Cutting between scc %d and %d\n",j,scc_id););
                        if(options->fuse == NO_FUSE) { 
                            cut_all_sccs(prog,ddg);
                        } else {
                            cut_between_sccs(prog,ddg,j,scc_id);
                            /* cut_smart(prog,ddg); */
                            /* You also need to cut between a successor node as well */
                            for (j=scc_id+1; j<ddg->num_sccs; j++) {
                                if (ddg_sccs_direct_connected(ddg, prog, scc_id, j)) {
                                    IF_DEBUG(printf("[colour SCC]: Cutting between scc %d and %d\n",scc_id,j););
                                    /* cut_between_sccs(prog, ddg, scc_id, j); */
                                    cut_all_sccs(prog, ddg);
                                    /* cut_smart(prog,ddg); */
                                    break;
                                }
                            }
                            break;
                        }
                    } else if (ddg_sccs_direct_connected(ddg, prog, scc_id, j)) {
                        IF_DEBUG(printf("[colour SCC]: Cutting between scc %d and %d\n",scc_id,j););

                        if (options->fuse == NO_FUSE) {
                            cut_all_sccs(prog,ddg);
                        }
                        else {
                            cut_between_sccs(prog, ddg, scc_id, j);
                        }
                        /* cut_smart(prog,ddg); */
                        break;
                    }
                }
            }
        }

        return true;
    }

    stmt_id = sccs[scc_id].vertices[stmt_pos];
    fcg_offset = ddg->vertices[stmt_id].fcg_stmt_offset;

    while(num_discarded!=nvar){
        j = get_next_min_vertex(fcg_offset, stmt_id, list, num_discarded, pv, prog);
        IF_DEBUG(printf("[Colour SCC] Trying Colouring dimension %d of statement %d with colour %d\n",j,stmt_id,c););

        v = fcg_offset+j;

        /* If the dimension is already coloured with a different colour. 
         * Else it tries to check if the existing colour is fine. This is done 
         * as opposed to undoing the existing colour and then redoing it in 
         * the next step once FCG is rebuilt */
        if(colour[v]>0 && colour[v]!=c){
            IF_DEBUG(printf("[Colour SCC]Dimension %d of statement %d already coloured with colour %d\n",j,stmt_id,colour[v]););
            list[num_discarded] = v;
            num_discarded++;
            continue;
        }

        /* Can not colour a vertex with a self edge. 
         * This dimension is not permutable */
        if (fcg->adj->val[v][v] != 0) {
            list[num_discarded] = v;
            num_discarded++;
            continue;
        }

        if (pv>=0 && is_adjecent(fcg,v,pv)) {
            list[num_discarded] = v;
            num_discarded++;
            continue;
        }

        /* Check if this is a valid colour */
        if (is_valid_colour(v,c,fcg,colour)) {
            colour[v] = c;
            /* If this is a valid colour, then try colouring the next vertex in the SCC */
            if(colour_scc(scc_id,colour,c,stmt_pos+1, v, prog)){
                IF_DEBUG(printf("[Colour SCC]Colouring dimension %d of statement %d with colour %d\n",j,stmt_id,c););
                return true;
            } else { 
                list[num_discarded] = v;
                num_discarded++;
                IF_DEBUG(printf("[Colour SCC] Unable to Colour dimension %d of statement %d with colour %d\n",j,stmt_id,c););
                /* Undo the colouring. Try the next vertex. */
                colour[v] = 0;
            }
        } else {
            colour[v] = 0;
            list[num_discarded] = v;
            num_discarded++;
        }
    }
    return false;
}



/* Colours all scc's with a colour c.  Returns true if the SCC's had to be distributed.  Else returns false */
bool colour_fcg_scc_based(int c, int *colour, PlutoProg *prog)
{
    int i,j,nsccs,prev_scc;
    bool is_distributed;
    Graph *ddg,*fcg;
    double t_start;

    ddg = prog->ddg;
    fcg = prog->fcg;
    nsccs = ddg->num_sccs;

    is_distributed = false;
    prev_scc = -1;


    for (i=0; i<nsccs; i++) {
        IF_DEBUG(printf("[colour_fcg_scc_based]: Colouring Scc %d of Size %d with colour %d\n",i,ddg->sccs[i].size, c););
        /* If colouring fails in the fist SCC */
        if(!colour_scc(i,colour,c,0, -1, prog)) {
            IF_DEBUG(printf("Unable to colour SCC %d\n",i););

            fcg = prog->fcg;
            /* In case of first scc, no inter scc deps can be satisfied. A permute 
             * preventing dependence has prevented colouring. 
             * Update the DDG whenever an inter SCC is satisfied dependence is
             * satisfied.  Note that dependencies that are satisfied by previous dimensions 
             * are updated in the DDG.  However, updating the FCG is delayed in order to 
             * account for permute preventing dependences.  Whenever the colouring fails, 
             * one has to update FCG with respect to the dependences that have already been 
             * satisfied along with the dependences those satisfied by the cut*/
            if (fcg->to_be_rebuilt == true || i == 0) {
                IF_DEBUG(printf("FCG Before Reconstruction\n"););
                IF_DEBUG(pluto_matrix_print(stdout, fcg->adj););

                if (options->fuse == NO_FUSE) {
                    cut_all_sccs(prog, ddg);
                }
                prog->fcg_colour_time += rtclock() - t_start;
                /* Current colour that is being used to colour the fcg is c */
                printf("FCG to be rebuilt due to a permute preventing dep: Colouring with colour %d\n",c);
                prog->fcg = build_fusion_conflict_graph(prog, colour, fcg->nVertices, c);

                t_start = rtclock();
                prog->fcg->num_coloured_vertices = fcg->num_coloured_vertices;
                /* need not update the FCG till the next hyperplane is found */
                prog->fcg->to_be_rebuilt = false;
                graph_free(fcg);
                fcg = prog->fcg;
                IF_DEBUG(printf("[Pluto]: Fcg After reconstruction\n"););
                IF_DEBUG( pluto_matrix_print(stdout, fcg->adj););
                /* Needed only if it is not the first SCC */
                if (i!=0) {
                    is_distributed = colour_scc(i, colour, c, 0, -1, prog);
                    if (!is_distributed){
                        /* Colouring was prevented by a fusion preventing dependence. 
                         * Therefore cut DDG then update FCG and then colour */
                        IF_DEBUG(printf("FCG Before Updating\n"););
                        IF_DEBUG(pluto_matrix_print(stdout, fcg->adj););
                        IF_DEBUG(printf("[colour_fcg_scc_based]:Total Number of SCCs %d\n",nsccs););

                        if (options->fuse == NO_FUSE) {
                            cut_all_sccs(prog,ddg);
                            /* TODO: Update this call testing is done */
                            update_fcg_between_sccs(fcg, 0, 0, prog);
                        } else {
                            for(j=prev_scc; j>=0; j--) {
                                if (ddg_sccs_direct_connected(ddg, prog, j, i)) {

                                    IF_DEBUG(printf("[colour_fcg_scc_based]:Cutting between SCC %d and %d\n",i,j););
                                    cut_between_sccs(prog,ddg,j,i);
                                    break;
                                }
                            }

                            update_fcg_between_sccs(fcg,prev_scc,i,prog);
                        }
                        IF_DEBUG(printf("DDG after Cut\n"););
                        IF_DEBUG(pluto_matrix_print(stdout, ddg->adj););
                        IF_DEBUG(printf("[Pluto] Colour_fcg_dim_based: Updating FCG\n"););

                        IF_DEBUG( printf("FCG after Updating \n"););
                        IF_DEBUG( pluto_matrix_print(stdout, fcg->adj););
                        is_distributed = colour_scc(i,colour,c,0, -1, prog);
                    }
                } else {
                    /* If the colouring of first SCC had failed previously */ 
                    is_distributed = colour_scc(i, colour, c, 0,  -1, prog);
                }
            } else {
                IF_DEBUG(printf("FCG Before Updating\n"););
                IF_DEBUG(pluto_matrix_print(stdout, fcg->adj););
                IF_DEBUG(printf("[Pluto] Colour_fcg_dim_based: Updating FCG\n"););
                if (options->fuse == NO_FUSE) {
                    cut_all_sccs(prog,ddg);

                    update_fcg_between_sccs(fcg, 0, 0, prog);
                } else {
                    for(j=prev_scc; j>=0; j--) {
                        if (ddg_sccs_direct_connected(ddg, prog, j, i)) {

                            IF_DEBUG(printf("[colour_fcg_scc_based]:Cutting between SCC %d and %d\n",i,j););
                            cut_between_sccs(prog,ddg,j,i);
                            break;
                        }
                    }
                    update_fcg_between_sccs(fcg,prev_scc,i,prog);
                }
                IF_DEBUG(printf("DDG after Cut\n"););
                IF_DEBUG(pluto_matrix_print(stdout, ddg->adj););
                IF_DEBUG( printf("FCG after Updating \n"););
                IF_DEBUG( pluto_matrix_print(stdout, fcg->adj););
                is_distributed = colour_scc(i,colour,c,0, -1, prog);
            }

            /* Needed in case of partial satisfaction */
            if (is_distributed == false) {
                printf("Num Deps satisfied with precise check %d\n",pluto_compute_dep_satisfaction_precise(prog));

                pluto_transformations_pretty_print(prog);
                pluto_compute_dep_directions(prog);
                /* pluto_compute_dep_satisfaction(prog); */
                pluto_print_dep_directions(prog);

                prog->fcg_colour_time += rtclock() - t_start;
                prog->fcg = build_fusion_conflict_graph(prog, colour, fcg->nVertices,c);
                t_start = rtclock();
                prog->fcg->num_coloured_vertices = fcg->num_coloured_vertices;

                /* need not update the FCG till the next hyperplane is found */
                prog->fcg->to_be_rebuilt = false;
                graph_free(fcg);
                fcg = prog->fcg;
                IF_DEBUG(printf("[Pluto]: Fcg After reconstruction\n"););
                IF_DEBUG( pluto_matrix_print(stdout, fcg->adj););
                is_distributed = colour_scc(i,colour,c,0, -1, prog);

            }
            assert (is_distributed == true);
        }
        fcg->num_coloured_vertices += ddg->sccs[i].size;
        prog->total_coloured_stmts[c-1] += ddg->sccs[i].size;
        prev_scc = i;
        prog->fcg_colour_time += rtclock() - t_start;
    }
    for (i=0; i<nsccs; i++) {
        if (prog->ddg->sccs[i].sol != NULL) {
            free(prog->ddg->sccs[i].sol);
            prog->ddg->sccs[i].sol = NULL;
        }
    }

    return is_distributed;
}

void find_permutable_dimensions_scc_based(int *colour, PlutoProg *prog)
{
    int i,j,num_coloured_dims,max_colours;
    Stmt **stmts;
    Graph *ddg;
    double t_start;

    max_colours = prog->nvar;
    stmts = prog->stmts;

    for (i=1; i<=max_colours; i++) {
        colour_fcg_scc_based(i, colour, prog);

        t_start = rtclock();

        num_coloured_dims = scale_shift_permutations(prog, colour, i-1);

        prog->fcg_dims_scale_time += rtclock() - t_start;
        if (num_coloured_dims == 0){
            printf ("[Pluto]: Num hyperplanes found: %d\n", prog->num_hyperplanes);
            printf("[Pluto]: This appears to be a bug in Pluto FCG based auto-transformation.\n");
            printf("[Pluto]: Transformation found so far\n");
            pluto_transformations_pretty_print(prog);
            pluto_print_colours(colour,prog);
            pluto_compute_dep_directions(prog);
            pluto_compute_dep_satisfaction(prog);
            pluto_print_dep_directions(prog);
            assert(0);

        }
        IF_DEBUG(printf ("[Pluto]: Num hyperplanes found: %d\n", prog->num_hyperplanes););
        prog->scaled_dims[i-1] = 1;


        prog->coloured_dims += num_coloured_dims;
        t_start = rtclock();
        for (j=0; j<num_coloured_dims; j++) {
            dep_satisfaction_update(prog,stmts[0]->trans->nrows-num_coloured_dims+j);
        }

        prog->fcg->to_be_rebuilt = 1;

        /* Recompute the SCC's in the updated DDG */
        ddg = prog->ddg;
        IF_DEBUG(printf("[Find_permutable_dims_scc_based]: Updating SCCs \n"););
        free_scc_vertices(ddg);

        /* You can update the DDG but do not update the FCG.  Doing otherwise will remove 
         * edges wich prevents permutation which is unsound */
        ddg_update(ddg, prog);
        IF_DEBUG(printf("DDG after colouring with colour %d\n",i););
        IF_DEBUG(pluto_matrix_print(stdout, ddg->adj););
        ddg_compute_scc(prog);
        compute_scc_vertices(ddg);
        IF_DEBUG2(pluto_transformations_pretty_print(prog););
        IF_DEBUG2(pluto_compute_dep_directions(prog););
        IF_DEBUG2(pluto_compute_dep_satisfaction(prog););
        IF_DEBUG2(pluto_print_dep_directions(prog););
    }
    /* All dimensions have been coloured but still there are some deps that 
     * need to be satisfied at the innermost level by distribution.  */
    if (i == max_colours+1 && !deps_satisfaction_check(prog)) {
        cut_all_sccs(prog, prog->ddg);
    }

    IF_DEBUG(printf("[Pluto] Colouring Successful\n"););
    IF_DEBUG(pluto_print_colours(colour,prog););


    return;
}


/*************************** Scaling Routines ******************/

/* Once the permutation is found, it finds the scling and shifting factors for the permtation
 * Scales the dimensions in the with colour c+1. Returns 1 if scaling
 * was successful. Else returns 0. */
int scale_shift_permutations(PlutoProg *prog, int *colour, int c)
{
    int j, k, stmt_offset, select;
    int nvar, npar;
    int nstmts;
    double t_start;
    PlutoConstraints *basecst,*coeffcst,*boundcst;
    int64 *sol;

    Stmt **stmts;

    nvar = prog->nvar;
    npar = prog->npar;
    stmts = prog->stmts;
    nstmts = prog->nstmts;


    t_start = rtclock();
    basecst = get_permutability_constraints(prog);
    assert (basecst->ncols == CST_WIDTH);

    boundcst = get_coeff_bounding_constraints(prog);
    pluto_constraints_add(basecst,boundcst);
    pluto_constraints_free(boundcst);


    coeffcst = pluto_constraints_alloc(basecst->nrows + (nstmts*nvar), basecst->ncols);
    coeffcst->nrows = basecst->nrows;
    coeffcst->ncols = basecst->ncols;
    assert (coeffcst->ncols == CST_WIDTH);


    IF_DEBUG(printf("Num stmts coloured with colour %d: %d\n", c+1, prog->total_coloured_stmts[c]););

    if (prog->total_coloured_stmts[c] == nstmts) {
        coeffcst = pluto_constraints_copy(coeffcst,basecst);

        /* Pick a colour that you would start with. This is buggy. You need to pick a colour*/
        select = c+1;
        IF_DEBUG(printf("[pluto] Finding Scaling factors for colour %d\n",select););
        /* Stmt offset in the colour array. This corresponds to FCG */
        stmt_offset = 0;
        /* Add CST_WIDTH number of cols and set appropriate constraints to 1 and set the rest to 0
         * These redundant cols are then removed. */

        for (j=0; j<nstmts; j++){
            for (k=0; k<nvar; k++){
                if (stmts[j]->is_orig_loop[k] && colour[stmt_offset+k]==select){
                    pluto_constraints_add_lb(coeffcst,npar+1+j*(nvar+1)+k,1);
                }
                else {
                    pluto_constraints_add_equality(coeffcst);
                    coeffcst->val[coeffcst->nrows-1][npar+1+j*(nvar+1)+k] = 1;
                }
            }
            stmt_offset += stmts[j]->dim_orig;
        }

        /* Solve the constraints to find the hyperplane at this level */
        t_start = rtclock();

        if (options->glpk) {
            double **val = NULL;
            int **index = NULL;
            int num_ccs, nrows;
            num_ccs = 0;

            PlutoMatrix *obj = construct_cplex_objective(coeffcst, prog);

            if (!options->ilp) {
                nrows = coeffcst->ncols-1-npar-1;
                populate_scaling_csr_matrices_for_pluto_program(&index, &val, nrows, prog);
                num_ccs = prog->ddg->num_ccs;
            }

            sol = pluto_prog_constraints_lexmin_glpk(coeffcst, obj, val, index, npar, num_ccs);
            if (!options->ilp) {
                for (j=0; j<nrows; j++) {
                    free(val[j]);
                    free(index[j]);
                }
                free(val);
                free(index);
            }
        } else {
            sol = pluto_constraints_lexmin(coeffcst, DO_NOT_ALLOW_NEGATIVE_COEFF);
        }
        /* sol = solve_scaling_constraints_glpk(coeffcst,colour,select,prog);  */

        if (sol != NULL) {
            IF_DEBUG(fprintf(stdout, "[pluto] find_permutable_hyperplanes: found a hyperplane\n"));
            /* num_sols_found++; */

            /* if (options->varliberalize) { */
            /*     for (j=0; j<ndeps; j++) { */
            /*         #<{(| Check if it has to be c or c+1 |)}># */
            /*         if(deps[j]->temp_across && c < deps[j]->fuse_depth  */
            /*                 && pluto_domain_equality(stmts[deps[j]->src],stmts[deps[j]->dest])) { */
            /*             for(k=0;k<nvar;k++) { */
            /*                 if(sol[npar+1+(deps[j]->src)*(nvar+1)+k] != sol[npar+1+(deps[j]->dest)*(nvar+1)+k]) */
            /*                     break; */
            /*             } */
            /*             if(k!=nvar || (sol[npar+1+(deps[j]->src)*(nvar+1)+nvar]!=sol[npar+1+(deps[j]->dest)*(nvar+1)+nvar])) { */
            /*                 printf("Cutting between SCCs to prevent illegal transformation with var-lib"); */
            /*                 cut_between_sccs(prog,ddg,stmts[deps[j]->src]->scc_id, stmts[deps[j]->dest]->scc_id); */
            /*             } */
            /*         } */
            /*     } */
            /* } */
            pluto_prog_add_hyperplane(prog, prog->num_hyperplanes, H_LOOP);

            for (j=0; j<nstmts; j++)    {
                Stmt *stmt = stmts[j];
                pluto_stmt_add_hyperplane(stmt, H_UNKNOWN, stmt->trans->nrows);
                for (k=0; k<nvar; k++)    {
                    stmt->trans->val[stmt->trans->nrows-1][k] = sol[npar+1+j*(nvar+1)+k];
                }
                /* No parameteric shifts */
                for (k=nvar; k<nvar+npar; k++)    {
                    stmt->trans->val[stmt->trans->nrows-1][k] = 0;
                }
                /* Constant loop shift */
                stmt->trans->val[stmt->trans->nrows-1][nvar+npar] = sol[npar+1+j*(nvar+1)+nvar];

                stmt->hyp_types[stmt->trans->nrows-1] =  
                    pluto_is_hyperplane_scalar(stmt, stmt->trans->nrows-1)?
                    H_SCALAR: H_LOOP;

            }
            free(sol);
            IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes-1););
            pluto_constraints_free(coeffcst);
            prog->scaling_cst_sol_time += rtclock()-t_start;
            return 1;
        } else {
            printf("[pluto] No Hyperplane found\n");
            pluto_constraints_free(coeffcst);
            prog->scaling_cst_sol_time += rtclock()-t_start;
            return 0;
        }

    } else {
        IF_DEBUG(printf("Not All statements have been coloured\n"););
        pluto_constraints_free(coeffcst);
        return 0;
    }
    /* } */
}

/* Routines that introduce loop skewing after loop permutations, loop skewing
 * and loop shifting transfomations have been found. */
bool get_negative_components(Dep *dep, bool *dims_with_neg_components, PlutoProg *prog, int level)
{
    int i;
    bool has_negative_comp;
    HyperplaneProperties *hProps;
    int loop_dims;

    hProps = prog->hProps;
    has_negative_comp = false;
    loop_dims = 0;
    /* printf("[get_negative_comp]: Level :%d \n",level); */
    for (i=0; i<prog->num_hyperplanes; i++){
        if(hProps[i].type == H_SCALAR && i < level){
            continue;
        }
        if(hProps[i].type == H_LOOP && i < level){
            loop_dims++;
            continue;
        }
        if(hProps[i].type == H_SCALAR && i >= level){
            continue;
            /* return has_negative_comp; */
        }
        if (dep->dirvec[i] == DEP_MINUS || dep->dirvec[i] == DEP_STAR ) {
            /* printf("dep %d has a negative component in level %d \n", dep->id,i); */
            dims_with_neg_components[loop_dims] = 1;
            has_negative_comp = true;
            break;
        }
        loop_dims++;
    }
    /* printf("Negative component for Dep %d : %d\n", dep->id, has_negative_comp); */
    return has_negative_comp;
}


bool* dims_to_be_skewed(PlutoProg *prog, int scc_id, bool *tile_preventing_deps, int level)
{
    int i, ndeps, nvar;
    Stmt **stmts;
    Dep *dep;
    bool* dims_with_neg_components;
    
    nvar = prog->nvar;
    ndeps = prog->ndeps;
    stmts = prog->stmts;

    dims_with_neg_components = (bool*) malloc (nvar*sizeof(bool*));
    bzero(dims_with_neg_components, nvar*sizeof(bool));


    /* For each dep find whether it is satisfied by a cut or loop; */
    for (i=0; i<ndeps; i++) {
        dep = prog->deps[i];
        if (!options->rar && IS_RAR(dep->type))
            continue;
        if (!(stmts[dep->src]->scc_id == scc_id) || !(stmts[dep->dest]->scc_id ==scc_id))
            continue;
        /* if (options->varliberalize && dep->skipdep) { */
        /*     continue; */
        /* } */
      
        if (dep_is_satisfied(dep)) {
            if(get_negative_components(dep,dims_with_neg_components, prog, level)){
                tile_preventing_deps[i] = 1;
            }
        }
    }
    return dims_with_neg_components;
}


bool* innermost_dep_satisfaction_dims(PlutoProg *prog, bool *tile_preventing_deps)
{
    int i, j, ndeps, loop_dims;
    Dep *dep;
    bool* sat_dim;
    HyperplaneProperties *hProps;

    sat_dim = (bool*)malloc(prog->nvar*sizeof(bool));
    bzero(sat_dim, (prog->nvar)*sizeof(bool));
    ndeps = prog->ndeps;
    hProps = prog->hProps;

    for (i=0; i<ndeps; i++) {
        dep = prog->deps[i];
        loop_dims = 0;
        if (tile_preventing_deps[i]) {
            /* Update this. Make src_dims to be the levels instead of dimensions. */
            for(j=0; j<prog->num_hyperplanes; j++) {
                if (j == dep->satisfaction_level) {
                    break;
                }
                else if(hProps[j].type == H_LOOP) {
                    loop_dims++;
                }
                /* if(hProps[j].type == H_SCALAR){ */
                /*     continue; */
                /* } */
                /* if(dep->dirvec[j] == DEP_MINUS || dep->dirvec[j] == DEP_STAR) */
                /*     break; */
                /* if(dep->dirvec[j] == DEP_PLUS ) */
                /*     innermost_sat_dim = loop_dims; */
            }
            sat_dim[loop_dims] = 1;
        }
    }
    return sat_dim;
}

PlutoConstraints *get_skewing_constraints(bool *src_dims, bool* skew_dims, int scc_id, PlutoProg* prog, int level, PlutoConstraints *skewCst)
{
    int i, j, nvar, npar, nstmts;
    Stmt **stmts;
    /* PlutoConstraints *skewCst, *boundcst; */

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;

    assert (skewCst->ncols == CST_WIDTH);
    
    for (i=0; i<nstmts; i++) {
        for(j=0; j<stmts[i]->dim_orig; j++){
            /* printf("Scc_ids: %d %d\n",stmts[i]->scc_id, scc_id); */
            if (src_dims[j] && stmts[i]->scc_id == scc_id) {
                /* printf("Coeff %d of Statemtnt %d >=1\n", j, i); */
                pluto_constraints_add_lb(skewCst, npar+1+ i*(nvar+1)+j, 1);
            }
            else {
                /* printf("Coeff %d of Statemtnt %d  = %lld\n", j, i, stmts[i]->trans->val[level][j]); */
                pluto_constraints_add_equality(skewCst);
                skewCst->val[skewCst->nrows-1][npar+1+i*(nvar+1)+j]= 1;
                /* Set the value of the current coeff to the one that you have already found */
                skewCst->val[skewCst->nrows-1][CST_WIDTH-1] = -stmts[i]->trans->val[level][j];
            }
        }
        /* printf("Adding Coeff bound for the constant of statement %d \n", i); */
        pluto_constraints_add_lb(skewCst, npar+1+i*(nvar+1)+nvar, 0); 
    }
    /* prog->cst_const_time += rtclock() - tstart; */
    return skewCst;
}


/* Introduce Skewing Transformations if necessary: Called only when using the FCG based appraoch */
void introduce_skew(PlutoProg *prog)
{
    int i, j, k, num_sccs, nvar, npar, nstmts, level,ndeps;
    int initial_cuts;
    Graph *orig_ddg;
    int64* sol;
    PlutoConstraints *skewingCst, *basecst, *const_dep_check_cst;
    HyperplaneProperties *hProps;
    Stmt **stmts;
    bool *src_dims, *skew_dims, tile_preventing_deps[prog->ndeps];
    double tstart;


    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;
    ndeps = prog->ndeps;

    /* If there are zero or one hyperpane then you dont need to skew */
    if (prog->num_hyperplanes <=1) {
        return;
    }
    assert (prog->hProps != NULL);
    hProps = prog->hProps;

    if (!options->silent) {
        printf("[Pluto]: Tileabilty with skew\n");
    }
    /* Reset dependence satisfaction */
    /* pluto_dep_satisfaction_reset(prog); */
    tstart = rtclock();
    pluto_compute_dep_directions(prog);
    printf("Dep direction computation time %0.6lf\n", rtclock()-tstart);

    tstart = rtclock();
    /* pluto_transformations_pretty_print(prog); */
    /* pluto_compute_dep_satisfaction(prog); */
    /* printf("Dep satisfaction computation time %0.6lf\n", rtclock()-tstart); */
    /* pluto_print_dep_directions(prog); */
    pluto_dep_satisfaction_reset(prog);

    /* tstart = rtclock(); */
    Graph *newDDG = ddg_create (prog);
    orig_ddg = prog->ddg;
    prog->ddg = newDDG;

    for(i=0; i<prog->ndeps; i++) {
        tile_preventing_deps[i] = 0;
    }

    initial_cuts = 0;
    for (level = 0; level< prog->num_hyperplanes; level++) {
        if(hProps[level].type == H_LOOP){
            break;
        }
        initial_cuts ++; 
        dep_satisfaction_update(prog, level);
    }

    /* printf("Initial cuts: %d \n", initial_cuts); */
    /* Needed to handle the case when there are no loops  */
    if(initial_cuts == prog->num_hyperplanes) {
        return;
    }
    basecst = get_permutability_constraints(prog);
    ddg_update(newDDG, prog);

    assert(level == initial_cuts);
    ddg_compute_scc(prog);
    num_sccs = newDDG->num_sccs;
    /* prog->cst_const_time += rtclock() -tstart; */
    /* printf("Pre processing before skewing %0.6lf\n", rtclock()-tstart ); */

    /* printf("Num SCCs: %d\n", num_sccs); */
    const_dep_check_cst = pluto_constraints_alloc(ndeps*nvar+1, CST_WIDTH);
    skewingCst = pluto_constraints_alloc(basecst->nrows+nstmts*(nvar+1), basecst->ncols);
    skewingCst->nrows = 0;
    skewingCst->ncols = CST_WIDTH;
    pluto_constraints_add(skewingCst, basecst);
    dep_satisfaction_update(prog,level);
    for (i =0; i<num_sccs; i++) {
        IF_DEBUG(printf("-------Looking for skews in SCC %d -----------------\n", i););
        tstart = rtclock();

        /* dep_satisfaction_update(prog,level); */
        /* if (!constant_deps_in_scc(i, level, const_dep_check_cst, prog)) { */
        /*  */
        /*     IF_DEBUG(printf("Scc %d has atleast one non constant dep \n",i);); */
        /*     continue; */
        /* } */
        IF_DEBUG(printf("-------Analyzing for skews in SCC %d ----------------\n", i););
        skew_dims = dims_to_be_skewed(prog, i, tile_preventing_deps, level);
        src_dims = innermost_dep_satisfaction_dims(prog, tile_preventing_deps);
        /* prog->cst_const_time += rtclock() - tstart; */
        /* printf("Src Dims\n"); */
        /* for(j = 0; j<nvar; j++) { */
        /*     printf("%d %d\n", j, src_dims[j]); */
        /* } */
        level++;
        for (; level <prog->num_hyperplanes; level++) {
            if (hProps[level].type != H_LOOP) {
                continue;
            }

            /* printf("Skew dims\n"); */
            /* for( j=0; j<nvar; j++) { */
            /*     printf("%d %d \n", j, skew_dims[j]); */
            /* } */
            int skew_dim = 0;
            for(j=initial_cuts; j<prog->num_hyperplanes; j++) {
                if(prog->hProps[j].type == H_LOOP && skew_dims[skew_dim] == 1) {
                    level = j;
                    break;
                }
                else if (prog->hProps[j].type == H_LOOP){
                    skew_dim++;
                }
            }

            /* Skewing has to be done at level j+1 */
            if (j==prog->num_hyperplanes) {
                /* printf("[Pluto]: No skews required for SCC %d \n", i); */
                /* prog->ddg = orig_ddg; */
                /* graph_free(newDDG); */
                break;
            }
            /* printf("Skewing at level :%d\n", level); */

            skewingCst->nrows = basecst->nrows;
            get_skewing_constraints(src_dims, skew_dims, i, prog, level, skewingCst);
            /* pluto_constraints_cplex_print(stdout,skewingCst); */

            /* Solve Skewing Constraints */
            //tstart = rtclock();

            if (options->glpk) {
                double **val = NULL;
                int **index = NULL;
                int num_ccs, nrows;
                num_ccs = 0;

                PlutoMatrix *obj = construct_cplex_objective(skewingCst, prog);

                if (!options->ilp) {
                    nrows = skewingCst->ncols-1-npar-1;
                    populate_scaling_csr_matrices_for_pluto_program(&index, &val, nrows, prog);
                    num_ccs = prog->ddg->num_ccs;
                }

                sol = pluto_prog_constraints_lexmin_glpk(skewingCst, obj, val, index, npar, num_ccs);
                if (!options->ilp) {
                    for (j=0; j<nrows; j++) {
                        free(val[j]);
                        free(index[j]);
                    }
                    free(val);
                    free(index);
                }
            } else {
                sol = pluto_constraints_lexmin(skewingCst, DO_NOT_ALLOW_NEGATIVE_COEFF);
            }
            /* sol = solve_scaling_constraints_glpk(coeffcst,colour,select,prog);  */

            /* sol = pluto_prog_constraints_lexmin(skewingCst, prog); */

            if(sol){
                /* Set the Appropriate coeffs in the transformation matrix */
                for (j=0; j<nstmts; j++) {
                    for (k = 0; k<nvar; k++){
                        stmts[j]->trans->val[level][k] = sol[npar+1+j*(nvar+1)+k];
                    }
                    /* No parametric Shifts */
                    for (k=nvar; k<nvar+npar; k++)    {
                        stmts[j]->trans->val[level][k] = 0;
                    }
                    /* The constant Shift */
                    stmts[j]->trans->val[level][nvar+npar] = sol[npar+1+j*(nvar+1)+nvar];

                }
                //printf("Trransformation Matrix After skewing \n");
                //pluto_transformations_pretty_print(prog);
                dep_satisfaction_update(prog, level);
                pluto_compute_dep_directions(prog);
                /* pluto_print_dep_directions(prog); */
                free(skew_dims);
                if(level < prog->num_hyperplanes-1) {
                    skew_dims = dims_to_be_skewed(prog, i, tile_preventing_deps, level+1);
                    free(src_dims);
                    src_dims = innermost_dep_satisfaction_dims(prog, tile_preventing_deps);
                }
            }
            else {
                /* printf("No Solution found \n"); */
                /* The loop nest is not tilable */
                free(skew_dims);
                break;
            }
        }
        free(src_dims);
        level = initial_cuts;
    }

    pluto_constraints_free(const_dep_check_cst);
    pluto_constraints_free(skewingCst);
    /* printf("All SCC's tested for skew\n"); */
    prog->ddg = orig_ddg;
    graph_free(newDDG);
    if (!options->silent) {
        printf("[Pluto]: Post processing skewing complete\n");
    }
    return;
}

