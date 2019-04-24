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

#if defined GLPK || defined GUROBI
int scale_shift_permutations(PlutoProg *prog, int *colour, int c);
double * pluto_fusion_constraints_feasibility_solve(PlutoConstraints *cst, PlutoMatrix *obj);
bool colour_scc(int scc_id, int *colour, int c, int stmt_pos, int pv, PlutoProg *prog);
void pluto_print_colours(int *colour,PlutoProg *prog);
PlutoMatrix *par_preventing_adj_mat;
PlutoMatrix *dep_distance_mat;

static double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

/* Constructs linear independence constraints for each statement in SCC scc_id. */
PlutoConstraints* dfp_get_scc_ortho_constraints (int *colour, int scc_id,
       PlutoProg *prog)
{
    int nvar, npar, nstmts,i, j, q, stmt_offset;
    Stmt **stmts;
    PlutoConstraints *indcst;
    bool has_dim_to_be_coloured;
    has_dim_to_be_coloured = false;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    indcst = NULL;

    stmts = prog->stmts;
    q = 0;

    stmt_offset = npar+1;

    for (i=0; i<nstmts; i++) {
        if (stmts[i]->scc_id == scc_id) {
            for (j=0; j<stmts[i]->dim_orig; j++) {
                if (colour[q] == 0) {
                    if (indcst == NULL) {
                        indcst = pluto_constraints_alloc(nstmts, CST_WIDTH);
                        indcst->nrows = 0;
                        indcst->ncols = CST_WIDTH;
                    }
                    indcst->val[indcst->nrows][stmt_offset+i*(nvar+1)+j] = 1;
                    has_dim_to_be_coloured = true;
                }
                q++;
            }
            if (has_dim_to_be_coloured == true) {
                indcst->val[indcst->nrows][CST_WIDTH-1] = -1;
                indcst->nrows ++;
            }
        } else {
            q += stmts[i]->dim_orig;
        }
        has_dim_to_be_coloured = false;
    }
    return indcst;
}

static int64 get_dep_dist_from_pluto_sol(double *sol, int npar)
{
    unsigned int i;
    int64 sum = 0;
    for (i=0; i<npar+1; i++) {
        sum = sum + ceil(sol[i]);
    }
    return sum;
}


/* A hyperplane is parallel if u+w is zero. */
/* Returns true if the solution represents a parallel hyperplane */
static bool is_lp_solution_parallel(double *sol, int npar)
{
    int i;
    double tmp;
    tmp = 0.0f;
    for (i = 0; i<npar+1; i++) {
        tmp += sol[i];
    }
    if (tmp == 0.0f)
        return true;
    else 
        return false;
}

/* Routine to mark parallel SCCs. This is called in dfp approach
 * when SCC clustering is disabled */
void mark_parallel_sccs(int *colour, PlutoProg* prog) 
{
    int i, num_sccs;
    PlutoConstraints* indcst, *boundcst, *permutecst;
    PlutoMatrix *obj;
    double *sol;

    num_sccs = prog->ddg->num_sccs;


    boundcst = get_coeff_bounding_constraints(prog);
    obj = construct_cplex_objective(boundcst, prog); 

    for (i=0; i<num_sccs; i++) {
        IF_DEBUG(printf("[pluto] Checking parallelism for SCC %d\n", i););
        sol = NULL;
        permutecst = get_scc_permutability_constraints(i, prog);
        indcst = dfp_get_scc_ortho_constraints(colour, i, prog);

        /* If there are no deps or if there are no linear independence
         * constraints then the scc is parallel*/
        if (indcst != NULL) {
            pluto_constraints_add(indcst, boundcst);
            if (permutecst != NULL) {
                pluto_constraints_add(indcst, permutecst);
            }
            sol = pluto_fusion_constraints_feasibility_solve(indcst, obj); 

            /* If sol is null, test again with a precise satisfaction check */
            if (sol == NULL) {
                pluto_compute_dep_satisfaction_precise(prog);
                pluto_transformations_pretty_print(prog);
                ddg_update(prog->ddg, prog);
                ddg_compute_scc(prog);
                assert(num_sccs == prog->ddg->num_sccs);
                free(permutecst);
                free(indcst);
                permutecst = get_scc_permutability_constraints(i, prog);
                indcst = dfp_get_scc_ortho_constraints(colour, i, prog);
                pluto_constraints_add(indcst, boundcst);
                if (permutecst != NULL) {
                    pluto_constraints_add(indcst, permutecst);
                }
                sol = pluto_fusion_constraints_feasibility_solve(indcst, obj);
            }
            /* There must exist a hyperplane for a scc that weakly satisfies
             * all depenences in that SCC*/
            assert (sol != NULL);
            if (is_lp_solution_parallel(sol, prog->npar)) {
                prog->ddg->sccs[i].is_parallel = 1;
                IF_DEBUG(printf("SCC %d is parallel \n",i););
            } else {
                prog->ddg->sccs[i].is_parallel = 0;
            }
            pluto_constraints_free(indcst);
        } else {
            /* The case where there are no more dimensions
             * to be found for the SCC */
            prog->ddg->sccs[i].is_parallel = 1;
        }
        prog->ddg->sccs[i].sol = sol;
        if (permutecst != NULL) {
            pluto_constraints_free(permutecst);
        }
    }
    pluto_matrix_free(obj);
    pluto_constraints_free(boundcst);
}

/* Print ids and number of parallel SCCs */
void print_parallel_sccs(Graph *ddg)
{
    int i, num_par_sccs;
    num_par_sccs = 0;
    printf("Ids Parallel SCCs:");
    for (i=0; i<ddg->num_sccs; i++) {
        if (ddg->sccs[i].is_parallel) {
            printf (" %d",i );
            num_par_sccs++;
        }
    }
    printf("\n");
    printf("Number of SCCs:%d\n", ddg->num_sccs);
    printf("Number of Parallel SCCs:%d\n", num_par_sccs);
}


/*************** FCG construction routines *************************/

/* Checks for feasibility of constraints.
 * If feasible then return the solution else returns NULL */
double* pluto_fusion_constraints_feasibility_solve(PlutoConstraints *cst,
        PlutoMatrix *obj)
{
    double* sol;
    sol = NULL;
    if (options->gurobi) {
#ifdef GUROBI
        sol = pluto_fcg_constraints_lexmin_gurobi(cst, obj);
#endif
    }else{
#ifdef GLPK
        sol = pluto_fcg_constraints_lexmin_glpk(cst, obj);
#endif
    }
    return sol;
}

/* Adds edges in FCG corresponding to the statements
 * represented by the nodes v1 and v2 in DDG*/
void fcg_add_pairwise_edges(Graph *fcg, int v1, int v2, PlutoProg *prog,
       int *colour, PlutoConstraints *boundcst, int current_colour,
       PlutoConstraints **conflicts, PlutoMatrix *obj)
{
    Graph *ddg;
    int i,j,ndeps,nstmts,nvar,npar, src_offset, dest_offset;
    double *sol;
    int row_offset, src_scc_id, dest_scc_id, fcg_offset1, fcg_offset2;
    bool check_parallel;
    Stmt **stmts;
    PlutoConstraints *conflictcst;
    Scc *sccs;

    ddg = prog->ddg;
    Dep **deps, *dep;
    ndeps = prog->ndeps;
    deps = prog->deps;

    stmts = prog->stmts;
    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;
    sccs = ddg->sccs;

    check_parallel = false;

    double tstart = rtclock();
    assert (*conflicts != NULL);
    conflictcst = *conflicts;

    prog->fcg_cst_alloc_time += rtclock() - tstart;
    row_offset = conflictcst->nrows-CST_WIDTH+1;
    /* conflictcst->ncols = CST_WIDTH; */

    /* conflictcst->val[row_offset][CST_WIDTH-1] = -1; */
    /* conflictcst->val[row_offset +1][CST_WIDTH-1] = -1; */

    for (i=0; i<ndeps; i++) {
        dep = deps[i];
        /*if (options->varliberalize && dep->skipdep) {
          continue;
          }*/
        if (dep_is_satisfied(dep)) {
            continue;
        }
        if ((dep->src == v1 && dep->dest == v2) ||
                (dep->src==v2 && dep->dest ==v1)) {
            if (dep->cst == NULL){
                compute_pairwise_permutability(dep,prog);
            }
            IF_DEBUG(printf("Adding Constraints for dependence %d\n",i););
            pluto_constraints_add(conflictcst, dep->cst);
        }
    }

    if (stmts[v1]->intra_stmt_dep_cst != NULL) {
        pluto_constraints_add(conflictcst, stmts[v1]->intra_stmt_dep_cst);
    }

    if (stmts[v2]->intra_stmt_dep_cst != NULL) {
        pluto_constraints_add(conflictcst, stmts[v2]->intra_stmt_dep_cst);
    }


    src_offset = npar+1+(nvar+1)*v1;
    dest_offset = npar+1+(nvar+1)*v2;

    fcg_offset1 = ddg->vertices[v1].fcg_stmt_offset;
    fcg_offset2 = ddg->vertices[v2].fcg_stmt_offset;
    
    src_scc_id = ddg->vertices[v1].scc_id;
    dest_scc_id = ddg->vertices[v2].scc_id;

    if (options->fuse == TYPED_FUSE && (src_scc_id != dest_scc_id) &&
            (sccs[src_scc_id].is_parallel || sccs[dest_scc_id].is_parallel)) {
        check_parallel = true;
    } else {
        check_parallel = false;
    }


    /* Solve Pluto LP by setting corresponding coeffs to 0 without any objective.
     * This is the check for fusability of two dimensions */
    for (i=0; i<stmts[v1]->dim_orig; i++) {
        /* note that the vertex should not be coloured. Even if the vertex has a 
         * self edge, it must be considered during construction of the FCG. This
         * is because,even after satisfying the permute preventing dep, it might
         * still prevent fusion. */
        if (colour[fcg_offset1+i] != 0 &&
               colour[fcg_offset1+i] != current_colour) {
            continue;
        }

        if (fcg->adj->val[fcg_offset1+i][fcg_offset1+i]==1) {
            /* Do not solve LPs if a dimenion of a
             * statement is not permutable */
            for (j=0;j<stmts[v2]->dim_orig; j++) {
                fcg->adj->val[fcg_offset1+i][fcg_offset2+j]=1;
            }
            continue;
        }

        /* Set the lower bound of i^th dimension of v1 to 1 */
        conflictcst->val[row_offset + src_offset+i][CST_WIDTH-1] = -1;
        conflictcst->is_eq[row_offset + src_offset+i] = 0;

        for (j=0; j<stmts[v2]->dim_orig; j++) {
            if (colour[fcg_offset2 + j] != 0 &&
                    colour[fcg_offset2 + j] != current_colour) {
                continue;
            }
            if (fcg->adj->val[fcg_offset1+i][fcg_offset1+i]==1) {
                fcg->adj->val[fcg_offset1+i][fcg_offset2+j]=1;
                continue;

            }

            /* Set the lower bound of i^th dimension of v1 to 1 */
            conflictcst->val[row_offset + dest_offset+j][CST_WIDTH-1] = -1;
            conflictcst->is_eq[row_offset + dest_offset+j] = 0;

            /* Check if fusing ith dimesion of the source with ith dimension
             * of the target is valid */
            prog->num_lp_calls ++;
            tstart = rtclock();
            sol = pluto_fusion_constraints_feasibility_solve(conflictcst, obj);
            prog->mipTime += rtclock()-tstart;

            /* If no solutions, then dimensions are not fusable.
             * Add an edge in the conflict graph. */
            if (sol == NULL) {
                IF_DEBUG(printf("Unable to fuse dimension %d of statement %d",i,v1););
                IF_DEBUG(printf(" with dimension %d of statement %d \n",j,v2););
                IF_DEBUG(printf(" Adding edge %d to %d in fcg\n", fcg_offset1+i,
                            fcg_offset2+j););
                fcg->adj->val[fcg_offset1+i][fcg_offset2+j] = 1;
            } else {
                if (check_parallel && !is_lp_solution_parallel(sol,npar)) {
                    /* Add parallelism preventing edge */
                    fcg->adj->val[fcg_offset1+i][fcg_offset2+j] = 1;
                }
                free(sol);
            }
            /* Unset the lowerbound for the coefficient of c_j.
             * The same constraint matrix is reused for all coeffs. */
            conflictcst->val[row_offset+dest_offset+j][CST_WIDTH-1] = 0;
            conflictcst->is_eq[row_offset+dest_offset+j] = 1;
        }

        /* Unset the lowerbound for the coefficient of c_i.
         * The same constraint matrix is reused for all coeffs. */
        conflictcst->val[row_offset + src_offset+i][CST_WIDTH-1] = 0;
        conflictcst->is_eq[row_offset + src_offset+i] = 1;
    }
    conflictcst->nrows = row_offset+CST_WIDTH-1;
    return;
}

/* Returns both intra and inter dependence constraints
 * for dependences between SCC1 and SCC2 */
PlutoConstraints* get_inter_scc_dep_constraints (int scc1, int scc2,
        PlutoProg* prog)
{
    int i, ndeps, src_stmt, dest_stmt;
    int src_scc, dest_scc, nrows, ncols;
    Dep **deps;
    Dep *dep;
    Stmt **stmts;
    PlutoConstraints *inter_scc_dep_cst;

    deps = prog->deps;
    ndeps = prog->ndeps;
    stmts = prog->stmts;

    inter_scc_dep_cst = NULL;

    for (i=0;i<ndeps; i++) {
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
        src_scc = stmts[src_stmt]->scc_id;
        dest_scc = stmts[dest_stmt]->scc_id;
        if ((src_scc == scc1 || src_scc == scc2) &&
                (dest_scc == scc1 || dest_scc == scc2)) {
            IF_DEBUG(printf("Computing Inter Scc deps for SCCs %d and %d\n",
                        src_scc, dest_scc););
            if (dep->cst==NULL) {
                compute_pairwise_permutability(dep,prog);
            }
            if (inter_scc_dep_cst == NULL) {
                nrows = dep->cst->nrows *ndeps;
                ncols = dep->cst->ncols;
                inter_scc_dep_cst = pluto_constraints_alloc(nrows, ncols);
                inter_scc_dep_cst->nrows = 0;
                inter_scc_dep_cst -> ncols = dep->cst->ncols;
            }
            pluto_constraints_add(inter_scc_dep_cst, dep->cst);
        }
    }
    return inter_scc_dep_cst;
}

/* Add inter SCC edges in the FCG in the clustered approach where
 * a vertex of the FCG corresponds to a dimension of a SCC */
void fcg_scc_cluster_add_inter_scc_edges (Graph* fcg, int *colour, PlutoProg *prog,
        PlutoConstraints *conflictcst, int current_colour, PlutoMatrix* obj)
{
    int i,j,num_sccs, scc1,scc2, row_offset, npar,nstmts, nvar, dim1, dim2;
    Graph *ddg;
    Scc *sccs;
    double tstart;
    double* sol;
    int scc1_fcg_offset, scc2_fcg_offset;
    int stmt1, stmt2, stmt1_offset, stmt2_offset;
    Stmt **stmts;
    bool check_parallel = false;


    PlutoConstraints* inter_scc_cst;

    ddg = prog->ddg;
    sccs = ddg->sccs;
    num_sccs = prog->ddg->num_sccs;
    nstmts = prog->nstmts;
    npar = prog->npar;
    nvar = prog->nvar;
    stmts = prog->stmts;

    for (scc1=0; scc1<num_sccs; scc1++) {
        scc1_fcg_offset = sccs[scc1].fcg_scc_offset;
        for (scc2=scc1+1; scc2<num_sccs; scc2++) {
            scc2_fcg_offset = sccs[scc2].fcg_scc_offset;

            if (!ddg_sccs_direct_connected(ddg, prog, scc1, scc2)) {
                continue;
            }

            inter_scc_cst = get_inter_scc_dep_constraints(scc1, scc2, prog);
            if ((sccs[scc1].is_parallel || sccs[scc2].is_parallel) &&
                    options->fuse == TYPED_FUSE) {
                check_parallel = true;
            }

            /* Conflict constraints are added at the end of
             * inter_scc_cst. Hence, we have inter_scc_constraints,
             * followed by bounding constraints, followed by dimension wise
             * constraints, which are set or unset based on the sccs between
             * which edges have to be added in the fcg */
            row_offset = conflictcst->nrows-CST_WIDTH+1+inter_scc_cst->nrows;

            /* Add conflict constraints at the end of inter_scc_constraints */
            /* pluto_constraints_cplex_print (stdout,conflictcst); */
            pluto_constraints_add(inter_scc_cst, conflictcst);

            /* Set the shifting lb of coefficient for each statement in SCC1 to 0 */
            for (i=0;i<sccs[scc1].size; i++) {
                stmt1 = sccs[scc1].vertices[i];
                stmt1_offset = npar+1+(nvar+1)*stmt1;
                inter_scc_cst->is_eq[row_offset+stmt1_offset+nvar] = 0;
            }
            /* Set the shifting lb of coefficient for each statement in SCC2 to 0 */
            for (j=0;j<sccs[scc2].size; j++) {
                stmt2 = sccs[scc2].vertices[j];
                stmt2_offset = npar+1+(nvar+1)*stmt2;
                inter_scc_cst->is_eq[row_offset+stmt2_offset+nvar] = 0;
            }
            /* Check for pairwise permutability of dimensions between scc1 and scc2 */
            for (dim1=0; dim1<sccs[scc1].max_dim; dim1++) {

                if (colour[scc1_fcg_offset+dim1] !=0 &&
                       colour[scc1_fcg_offset+dim1] != current_colour) {
                    continue;
                }
                /* If there is a self edge on this vertex, then do not
                 * solve LP's. Just add edges to all dimensions of SCC2. */
                if (fcg->adj->val[scc1_fcg_offset+dim1][scc1_fcg_offset+dim1] == 1) {
                    for (dim2=0; dim2<sccs[scc2].max_dim; dim2++) {
                        fcg->adj->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 1;
                    }
                    continue;
                }
                /* Set lb of dim1 of each statement in Scc 1 */
                for (i=0; i<sccs[scc1].size; i++) {
                    stmt1 = sccs[scc1].vertices[i];
                    if (dim1 >stmts[stmt1]->dim_orig) {
                        continue;
                    }
                    stmt1_offset = npar+1+(nvar+1)*stmt1;
                    inter_scc_cst->val[row_offset+stmt1_offset+dim1][CST_WIDTH-1] = -1;
                    inter_scc_cst->is_eq[row_offset+stmt1_offset+dim1] = 0;
                }
                for (dim2=0; dim2<sccs[scc2].max_dim; dim2++) {
                    /* Set the lower bounds of dimensions of each statement in SCC2 */
                    if (colour[scc2_fcg_offset+dim2] != 0 &&
                            colour[scc2_fcg_offset+dim2] != current_colour) {
                        continue;
                    }
                    if (fcg->adj->val[scc2_fcg_offset+dim2][scc2_fcg_offset+dim2]==1) {
                        fcg->adj->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 1;
                        continue;
                    }
                    for (j=0; j<sccs[scc2].size; j++) {
                        stmt2 = sccs[scc2].vertices[j];
                        if (dim2<=stmts[stmt2]->dim_orig) {
                            stmt2_offset = npar+1+(nvar+1)*stmt2;
                            inter_scc_cst->val[row_offset+stmt2_offset+dim2][CST_WIDTH-1] = -1;
                            inter_scc_cst->is_eq[row_offset+stmt2_offset+dim2] = 0;
                        }
                    }

                    /* Check if fusing ith dimesion of the source with i^th
                     * dimension of the target is valid */
                    prog->num_lp_calls ++;
                    tstart = rtclock();
                    sol = pluto_fusion_constraints_feasibility_solve(inter_scc_cst, obj);
                    prog->mipTime += rtclock()-tstart;

                    /* If no solutions, then dimensions are not fusable.
                     * Hence add an edge in the conflict graph. */
                    if (sol == NULL) {
                        IF_DEBUG(printf("Unable to fuse dim %d of scc %d",dim1, scc1););
                        IF_DEBUG(printf(" with dim %d of scc %d \n",dim2 ,scc2););
                        IF_DEBUG(printf(" Adding edge %d to %d in fcg\n",
                                    scc1_fcg_offset+dim1,scc2_fcg_offset+dim2););
                        fcg->adj->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 1;
                    } else {
                        if (options->lpcolour) {
                            dep_distance_mat->val[i][j] = get_dep_dist_from_pluto_sol(sol, npar);
                        }
                        if (check_parallel && !is_lp_solution_parallel(sol,npar)) {
                            IF_DEBUG(printf("Adding Parallelism preventing edge"););
                            IF_DEBUG(printf("%d to %d in fcg \n", scc1_fcg_offset+dim1,
                                       scc2_fcg_offset+dim2););
                            par_preventing_adj_mat->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 1;
                        }
                        free(sol);
                    }
                    /* Reset the lower bounds of dimensions of
                     * each statement in SCC2 */
                    for (j=0; j<sccs[scc2].size; j++) {
                        stmt2 = sccs[scc2].vertices[j];
                        if (dim2<=stmts[stmt2]->dim_orig) {
                            stmt2_offset = npar+1+(nvar+1)*stmt2;
                            inter_scc_cst->val[row_offset+stmt2_offset+dim2][CST_WIDTH-1] = 0;
                            inter_scc_cst->is_eq[row_offset+stmt2_offset+dim2] = 1;
                        }
                    }
                }
                /* Reset lb of dim1 of each statement in Scc 1 */
                for (i=0; i<sccs[scc1].size; i++) {
                    stmt1 = sccs[scc1].vertices[i];
                    if (dim1<=stmts[stmt1]->dim_orig) {
                        stmt1_offset = npar+1+(nvar+1)*stmt1;
                        inter_scc_cst->val[row_offset+stmt1_offset+dim1][CST_WIDTH-1] = 0;
                        inter_scc_cst->is_eq[row_offset+stmt1_offset+dim1] = -1;
                    }
                }
            }
            pluto_constraints_free(inter_scc_cst);
        }
    }
}

/* Computes intra statement dependence constraints for every unstisfied dependence */
void compute_intra_stmt_deps(PlutoProg *prog)
{
    int ndeps, src_stmt,dest_stmt,i;
    int nrows, ncols;
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
            IF_DEBUG(printf("Computing intra stmt deps stmt: %d\n",src_stmt););
            if (dep->cst==NULL) {
                compute_pairwise_permutability(dep, prog);
            }
            if (stmt->intra_stmt_dep_cst == NULL) {
                nrows = dep->cst->nrows;
                ncols = dep->cst->ncols;
                stmt->intra_stmt_dep_cst = pluto_constraints_alloc(nrows, ncols);
                stmt->intra_stmt_dep_cst->nrows = nrows;
                stmt->intra_stmt_dep_cst->ncols = ncols;
                pluto_constraints_copy(stmt->intra_stmt_dep_cst, dep->cst);
            } else {
                pluto_constraints_add(stmt->intra_stmt_dep_cst, dep->cst);
            }
        }
    }
}

/* Computes dependence constraints for all dependences in the given SCC */
PlutoConstraints* compute_intra_scc_dep_cst(int scc_id, PlutoProg *prog)
{
    int ndeps, src_scc,dest_scc,i;
    int nrows, ncols;
    Dep **deps;
    Dep *dep;
    Stmt **stmts;
    PlutoConstraints *intra_scc_dep_cst = NULL;

    deps = prog->deps;
    ndeps = prog->ndeps;
    stmts = prog->stmts;

    for (i=0; i<ndeps; i++) {
        dep = deps[i];
        if (options->rar ==0 && IS_RAR(dep->type)) {
            continue;
        }

        if (dep_is_satisfied(dep)) {
            continue;
        }

        src_scc = stmts[dep->src]->scc_id;
        dest_scc = stmts[dep->dest]->scc_id;
        if (src_scc == scc_id && dest_scc == scc_id) {
            IF_DEBUG(printf("Computing intra scc deps for Scc: %d\n",scc_id););
            if (dep->cst==NULL) {
                compute_pairwise_permutability(dep,prog);
            }
            if (intra_scc_dep_cst == NULL) {
                nrows = dep->cst->nrows*ndeps;
                ncols = dep->cst->ncols;
                intra_scc_dep_cst = pluto_constraints_alloc(nrows, ncols);
                intra_scc_dep_cst->nrows = 0;
                intra_scc_dep_cst->ncols = dep->cst->ncols;
            } 
            pluto_constraints_add(intra_scc_dep_cst, dep->cst);
        }
    }
    return intra_scc_dep_cst;
}

/* Adds permute preventing edges for intra statement dependences.  These edges
 * are added as self loops on FCG vertices. These vertices can not be coloured 
 * until the self loops are removed by reconstruction of the FCG.
 * Assumes that there are no loop shifts */
void add_permute_preventing_edges(Graph* fcg, int *colour, PlutoProg *prog, 
        PlutoConstraints* boundcst, int current_colour, PlutoMatrix *obj)
{
    int nstmts,nvar,npar,i,j,stmt_offset,fcg_stmt_offset;
    int nrows;
    double *sol, tstart;
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
        if (stmts[i]->intra_stmt_dep_cst == NULL) {
            fcg_stmt_offset += stmts[i]->dim_orig;
            continue;
        }
        /* Constraints to check permutability are added in the first row */
        coeff_bounds = pluto_constraints_alloc(1,CST_WIDTH);
        coeff_bounds->nrows = 0;
        coeff_bounds->ncols = CST_WIDTH;

        /* Add intra statement dependence constraints
         * and bounding constraints */
        intra_stmt_dep_cst = stmts[i]->intra_stmt_dep_cst;
        pluto_constraints_add(coeff_bounds,boundcst);
        pluto_constraints_add(coeff_bounds,intra_stmt_dep_cst);

        stmt_offset = (npar+1)+ i*(nvar+1);

        for (j=0; j<stmts[i]->dim_orig; j++) {
            if (colour[fcg_stmt_offset+j] !=0 &&
                    colour[fcg_stmt_offset+j] != current_colour) {
                continue;
            }
            IF_DEBUG(printf("Checking permutability of dimension %d",j););
            IF_DEBUG(printf(" of statement %d \n",i););
            /* Not an equality constraint. Set the lower bound to 1. */
            coeff_bounds->is_eq[nrows+stmt_offset+j] = 0;
            coeff_bounds->val[nrows+stmt_offset+j][CST_WIDTH-1] = -1;
            /* coeff_bounds->val[0][stmt_offset+j] = 1; */
            prog->num_lp_calls++;

            tstart = rtclock();
            sol = pluto_fusion_constraints_feasibility_solve(coeff_bounds,obj);
            prog->mipTime += rtclock()-tstart;
            /* If the constraints are infeasible then add a self edge in the FCG */

            if (sol == NULL) {
                IF_DEBUG(printf("Dimension %d of stmt %d is not permutable",j, i););
                fcg->adj->val[fcg_stmt_offset+j][fcg_stmt_offset+j] = 1;
            } else {
                free(sol);
            }
            /* reset the coeff bound of this dimension */
            coeff_bounds->is_eq[nrows+stmt_offset+j] = 1;
            coeff_bounds->val[nrows+stmt_offset+j][CST_WIDTH-1] = 0;
        }
        pluto_constraints_free(coeff_bounds);
        fcg_stmt_offset += stmts[i]->dim_orig;
    }
}

/* Adds permute preventing edges for intra SCC dependences. These edges
 * are added as self loops on FCG vertices. These vertices can not be coloured 
 * until the self loops are removed by reconstruction of the FCG.
 * Assumes that there are no loop shifts */
void fcg_scc_cluster_add_permute_preventing_edges(Graph* fcg, int *colour, 
        PlutoProg *prog, PlutoConstraints* boundcst, int current_colour, PlutoMatrix *obj)
{
    int nstmts,nvar,npar,i,j,k, stmt_offset,fcg_scc_offset;
    int nrows, stmt_id;
    double *sol, tstart;
    Stmt **stmts;
    Scc *sccs;
    int num_sccs;
    PlutoConstraints *intra_scc_dep_cst, *coeff_bounds;

    nstmts = prog->nstmts;
    nvar = prog->nvar;
    npar = prog->npar;
    stmts = prog->stmts;
    num_sccs=prog->ddg->num_sccs;
    sccs = prog->ddg->sccs;
    nrows = boundcst->nrows-CST_WIDTH+1;

    /* Compute the intra Scc dependence constraints */
    for (i=0; i<num_sccs; i++) {
        sccs[i].is_parallel = 1;
        intra_scc_dep_cst = compute_intra_scc_dep_cst(i, prog);

        if (intra_scc_dep_cst!=NULL) {
            sccs[i].is_parallel = 0;
            /* Constraints to check permutability are added in the beginning */
            coeff_bounds = pluto_constraints_alloc(1,CST_WIDTH);
            coeff_bounds->nrows = 0;
            coeff_bounds->ncols = CST_WIDTH;

            pluto_constraints_add(coeff_bounds, boundcst);

            pluto_constraints_add(coeff_bounds,intra_scc_dep_cst);

            fcg_scc_offset = sccs[i].fcg_scc_offset;

            fcg_scc_offset = sccs[i].fcg_scc_offset;
            for (j=0; j<sccs[i].max_dim; j++) {
                /* Check for permutability only if the current scc vertex in the 
                 * fcg is not coloured or has the current colour */
                if (colour[fcg_scc_offset + j]==0 || colour[fcg_scc_offset+j] == current_colour) {
                    IF_DEBUG(printf("[Permute_preventing_edges]: Checking permutability of dimension %d of Scc %d \n",j,i););
                    /* Set the lower bounds of the jth coefficient for all the statments in scc i to 1. */
                    for (k=0; k<sccs[i].size; k++) {
                        stmt_id = sccs[i].vertices[k];
                        if (stmts[stmt_id]->is_orig_loop[j]) {
                            stmt_offset = npar+1+stmt_id*(nvar+1)+j;
                            coeff_bounds->is_eq[nrows+stmt_offset] = 0;
                            coeff_bounds->val[nrows+stmt_offset][CST_WIDTH-1] = -1;
                        }
                    }
                    prog->num_lp_calls++;

                    tstart = rtclock();
                    sol = pluto_fusion_constraints_feasibility_solve(coeff_bounds,obj);
                    prog->mipTime += rtclock()-tstart;

                    /* If the constraints are infeasible then add a self edge in the FCG */
                    if (sol == NULL) {
                        IF_DEBUG(printf("Dimension %d of scc %d is not permutable\n",j,i););
                        fcg->adj->val[fcg_scc_offset+j][fcg_scc_offset+j] = 1;
                    } else {
                        if (options->fuse == TYPED_FUSE) {
                            if (! is_lp_solution_parallel(sol, prog->npar)) {
                                par_preventing_adj_mat->val[fcg_scc_offset+j][fcg_scc_offset+j] = 1;

                            } else {
                                sccs[i].is_parallel = 1;
                            }
                        }
                        IF_DEBUG(printf("Dimension %d of scc %d is permutable\n",j,i););
                        free(sol);
                    }
                    /* reset the coeff bound of this dimension for all statements in the SCC*/
                    for (k=0; k<sccs[i].size; k++) {
                        stmt_id = sccs[i].vertices[k];
                        if (j<=stmts[stmt_id]->dim_orig) {
                            stmt_offset = npar+1+stmt_id*(nvar+1)+j;
                            coeff_bounds->is_eq[nrows+stmt_offset] = 1;
                            coeff_bounds->val[nrows+stmt_offset][CST_WIDTH-1] = 0;
                        }
                    }
                }
            }
            pluto_constraints_free(coeff_bounds);
        }
        pluto_constraints_free(intra_scc_dep_cst);
    }
}

/* Updates FCG between SCCs in the clustered approach.
 * Edges that from any vertex less than scc2 to an SCC with id 
 * greater than or equal to scc2 will be removed. 
 * Parallelism preventing edges are removed as well in 
 * case of typed fuse */
void update_scc_cluster_fcg_between_sccs(Graph *fcg, int scc1, int scc2, PlutoProg *prog)
{
    int i, j, num_sccs, dim1, dim2;
    Graph *ddg;
    Scc *sccs;
    int scc1_fcg_offset, scc2_fcg_offset;
    ddg = prog->ddg;
    sccs = ddg->sccs;
    num_sccs = ddg->num_sccs;
    assert (scc1 != scc2);

    if (options->fuse == NO_FUSE) {
        for (i=0; i<num_sccs; i++) {
            scc1_fcg_offset = sccs[i].fcg_scc_offset;
            for (dim1 = 0; dim1<=sccs[i].max_dim ; dim1++) {
                for (j=0; j<num_sccs; j++) {
                    scc2_fcg_offset = sccs[j].fcg_scc_offset;
                    for (dim2=0; dim2<=sccs[j].max_dim; dim2++) {
                        /* No fusion. Hence all sccs are cut. Therefore remove in inter scc edges in FCG */
                        if (i!=j) {
                            fcg->adj->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 0;
                        }
                    }
                }
            }
        }
    } else {
        /* Update fcg only between scc1 and scc2 */
        IF_DEBUG(printf("Updating FCG between SCCs%d and %d\n",scc1,scc2););
        for (i=0; i<scc2; i++) {
            scc1_fcg_offset = sccs[i].fcg_scc_offset;
            for (dim1 =0; dim1<sccs[i].max_dim; dim1++) {
                for (j=scc2; j<num_sccs; j++) {
                    scc2_fcg_offset = sccs[j].fcg_scc_offset;
                    for (dim2 =0; dim2<sccs[j].max_dim; dim2++) {
                        fcg->adj->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 0;
                        fcg->adj->val[scc2_fcg_offset+dim2][scc1_fcg_offset+dim1] = 0;
                        if (options->fuse == TYPED_FUSE) {
                            par_preventing_adj_mat->val[scc1_fcg_offset+dim1][scc2_fcg_offset+dim2] = 0;
                            par_preventing_adj_mat->val[scc2_fcg_offset+dim2][scc1_fcg_offset+dim1] = 0;
                        }
                    }
                }
            }
        }
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

    if (!(options->fuse == TYPED_FUSE)) {
        /* This assertion might not hold in case of typed fuse */
        assert (fcg->to_be_rebuilt == false);
    }
    tstart = rtclock();

    if (nstmts == 1) {
        return;
    }

    if (options->scc_cluster) {
        update_scc_cluster_fcg_between_sccs(fcg, scc1, scc2, prog);
        prog->fcg_update_time += rtclock() - tstart;
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
                    for (k=0; k<stmts[i]->dim_orig; k++) {
                        for (l=0; l<stmts[j]->dim_orig; l++) {
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
            for (j=0; j<nstmts; j++) {
                if ((stmts[i]->scc_id >= scc2 && stmts[j]->scc_id<scc2) ||
                        (stmts[j]->scc_id>=scc2 && stmts[i]->scc_id < scc2)) {
                    stmt_offset1 = ddg->vertices[i].fcg_stmt_offset;
                    stmt_offset2 = ddg->vertices[j].fcg_stmt_offset;
                    for (k=0; k<stmts[i]->dim_orig; k++) {
                        for (l=0; l<stmts[j]->dim_orig; l++) {
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

void fcg_add_intra_scc_edges(Graph *fcg, PlutoProg *prog)
{
    Graph *ddg;
    int i, j, k ,num_sccs, scc_offset;

    ddg = prog->ddg;
    num_sccs = ddg->num_sccs;
    scc_offset = 0;

    for (i=0; i<num_sccs; i++) {
        for (j=0;j<ddg->sccs[i].max_dim; j++) {
            for (k =j+1; k<ddg->sccs[i].max_dim; k++) {
                fcg->adj->val[scc_offset+j][scc_offset+k] = 1;
                fcg->adj->val[scc_offset+k][scc_offset+j] = 1;
            }
        }
        scc_offset += ddg->sccs[i].max_dim;
    }
    return;
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
    double t_start;
    PlutoConstraints *boundcst, **conflicts;
    PlutoMatrix *obj;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;

    ddg = prog->ddg;

    t_start = rtclock();

    fcg = graph_alloc(num_nodes);

    if (options->fuse == TYPED_FUSE) {
        par_preventing_adj_mat = pluto_matrix_alloc (num_nodes, num_nodes);
        for (i=0; i<num_nodes; i++) {
            bzero(par_preventing_adj_mat->val[i], num_nodes*sizeof(int64));
        }
    }

    if (options->lpcolour) {
        dep_distance_mat = pluto_matrix_alloc (num_nodes, num_nodes);
        for (i=0; i<num_nodes; i++) {
            bzero(dep_distance_mat->val[i], num_nodes*sizeof(int64));
        }
    }

    boundcst = get_coeff_bounding_constraints(prog);

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

    /* The last CST_WIDTH-(npar+1) number of rows, correspond to equality constraints.
     * These are changed during dimension wise computation of edges of the FCG. The 
     * equality constraints are used to set the transformation coeffs to zero*/
    for (i=npar+1; i<CST_WIDTH-1; i++) {
        (*conflicts)->is_eq[nrows+i] = 1;
        (*conflicts)->val[nrows+i][i] = 1;
    }
    
    /* Add premutation preventing intra statement dependence edges in the FCG.
     * These are self loops on vertices of the FCG. */ 
    if (options->scc_cluster) {
        fcg_scc_cluster_add_permute_preventing_edges(fcg, colour, prog, *conflicts, current_colour, obj);
    } else {
        add_permute_preventing_edges(fcg, colour, prog, *conflicts, current_colour, obj);
    }

    /* Add inter statement fusion and permute preventing edges.  */
    if (options->fuse == TYPED_FUSE && !options->scc_cluster) {
        /* The lp solutions are found and the parallel sccs are marked. 
         * However marking is only used in parallel case of typed fuse only */
        mark_parallel_sccs(colour, prog);
        IF_DEBUG(print_parallel_sccs(prog->ddg););
    }

    if (options->scc_cluster) {
        fcg_scc_cluster_add_inter_scc_edges (fcg, colour, prog, *conflicts, current_colour, obj);
    } else {
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
    }

    pluto_matrix_free(obj);

    if (options->scc_cluster) {
        fcg_add_intra_scc_edges(fcg,prog);
    } else {
        /* Add egdes between different dimensions of the same statement */
        stmt_offset = 0;
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
    }

    pluto_constraints_free(boundcst);
    pluto_constraints_free(*conflicts);
    free(conflicts);
    prog->fcg_const_time += rtclock() - t_start;

    IF_DEBUG(printf("FCG \n"););
    IF_DEBUG(pluto_matrix_print(stdout, fcg->adj));
    if(options->fuse == TYPED_FUSE) {
        IF_DEBUG(printf("Parallelism preventing edges \n"));
        IF_DEBUG(pluto_matrix_print(stdout, par_preventing_adj_mat));
    }

    IF_DEBUG(printf("[Pluto] Build FCG: Total number of LP calls in building the FCG: %ld\n",prog->num_lp_calls););
    return fcg;
}

/******************  FCG Colouring Routines **********************************/

/* Prints colour of each vertex of the FCG */
void pluto_print_colours(int *colour,PlutoProg *prog)
{
    int nstmts, i, j, stmt_offset;
    int color;
    Stmt **stmts;

    nstmts = prog->nstmts;
    stmts = prog->stmts;

    stmt_offset = 0;

    if (options->scc_cluster) {
        for (i=0; i<prog->ddg->num_sccs;i++){
            for (j=0; j<prog->ddg->sccs[i].max_dim;j++){
                color = colour[stmt_offset+j];
                printf("Colour of dimension %d of Scc %d: %d\n",j,i,color);
            }
            stmt_offset += j;
        }
        return;
    }
    for (i=0; i<nstmts; i++) {
        for (j=0; j<stmts[i]->dim_orig; j++) {
            color = colour[stmt_offset+j];
            printf("Colour of Dimension %d of Stmt %d: %d\n",j,i,color);
        }
        stmt_offset += j;
    }
}

/* Check if it is valid to give colour c to a vertex v in the fcg.
 * Colour is the array containing the colours assigned to each vertex */
bool is_valid_colour(int v, int c, Graph *fcg, int * colour, bool is_parallel)
{
    int i, fcg_nVertices;
    bool par_preventing_edge;
    fcg_nVertices = fcg->nVertices;
    for (i=0; i<fcg_nVertices; i++) {
        if ((fcg->adj->val[i][v]==1 || fcg->adj->val[v][i]==1) && colour[i]==c) {
            return false;
        }

        if (is_parallel) {
            par_preventing_edge = par_preventing_adj_mat->val[v][i] ||
                par_preventing_adj_mat->val[i][v];
            if(par_preventing_edge && colour[i] == c) {
                return false;
            }
        }
    }
    return true;
}

/* Checks if a vertex was decided to be not a suitable candidate for colouring */
bool is_discarded(int v, int list[], int num)
{
    int i;
    for (i=0; i<num; i++){
        if (list[i]==v)
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

    for (i=0; i<stmts[stmt_id]->dim_orig; i++) {
        if (!is_discarded(fcg_stmt_offset+i, list, num)) {
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
    }
    return min;
}

/* Fix: Modify this routine to handle single SCC case.
 * Routine is called with unclustered approach with typed fuse. Currently not supported. 
 * It returns the common parallel dimensions between SCC1 and SCC2 by checking the 
 * non zero component in the LP solution representing the parallel hyperplanes of scc1 and scc2*/
int* get_common_parallel_dims_for_sccs(Scc scc1, Scc scc2, PlutoProg *prog)
{
    int i, j, stmt1, stmt2, npar, nvar, stmt_offset;
    int *parallel_dims;
    Stmt **stmts;
    Graph* ddg;
    ddg = prog->ddg;

    nvar = prog->nvar;
    stmts = prog->stmts;
    npar = prog->npar;

    stmt1 = -1;
    stmt2 = -1;
    parallel_dims = NULL;

    for (i=0; (i<scc1.size) && (stmt1 == -1); i++) {
        for (j=0; j<scc2.size; j++) {
            if (is_adjecent(ddg, scc1.vertices[i], scc2.vertices[j])) {
                stmt1 = scc1.vertices[i];
                stmt2 = scc2.vertices[j];
                break;
            }
        }
    }
    assert((stmt1 >= 0) && (stmt2 >= 0));
    
    stmt_offset = npar+1;
    for (i=0; i<nvar ; i++) {
        if (stmts[stmt1]->is_orig_loop[i] && stmts[stmt2]->is_orig_loop[i]) {
            if ((scc1.sol[stmt_offset+stmt1*(nvar+1)+i] > 0.0f) 
                    && (scc2.sol[stmt_offset+stmt2*(nvar+1)+i] > 0.0f)) {
                if (parallel_dims == NULL) {
                    parallel_dims = (int*) malloc (sizeof(int)*nvar);
                    bzero(parallel_dims, nvar*sizeof(int));
                }
                parallel_dims[i] = 1;
            }
        }
    }
    return parallel_dims;
}

/* Check if there is a predecessor of SCC2 which is numbered greater than SCC1 and has not been coloured */
bool is_convex_scc(int scc1, int scc2, Graph *ddg, PlutoProg * prog)
{
    int i;
    Scc *sccs;
    sccs = prog->ddg->sccs;
    for (i=scc1+1; i<scc2; i++) {
        if (ddg_sccs_direct_connected(ddg, prog, i,scc2)) {
            /* In case of typed fuse, this scc may have already been coloured */
            if (options->fuse == TYPED_FUSE && sccs[i].is_scc_coloured) {
                continue;
            }
            return false;
        }
    }
    return true;
}

/* The following routine isn't used anymore. Can be removed in future commits if found unnecessary. */
/* Currently typedfuse is supported only with SCC clustering heuristics */
#if 0
/* In the unclustered approach this colouring heuristic chooses 
 * parallel hyperlanes to be found whenever they exist */
bool colour_scc_from_lp_solution_with_parallelism (int scc_id, int *colour, PlutoProg *prog, int c)
{
    int i, nvar,j ;
    Graph *ddg;
    bool is_successful;

    /* Parallel dims of the current SCC */
    int *parallel_dims;

    nvar = prog->nvar;
    ddg = prog->ddg;

    assert (prog->ddg->sccs[scc_id].is_parallel == 1);
    parallel_dims = NULL;

    /* Looks for SCCs that are connected that have a common dimension */
    for (i=scc_id+1; i< prog->ddg->num_sccs; i++) {
        if (ddg_sccs_direct_connected (ddg, prog, scc_id, i) &&
                ddg->sccs[i].is_parallel && is_convex_scc(scc_id,i, ddg, prog)) {
/* Todo: Add greedy heuristic here that looks for dimensions that have maximum dimensions that can be fused. */
            parallel_dims = get_common_parallel_dims_for_sccs(ddg->sccs[scc_id], ddg->sccs[i], prog);
        }
        if (parallel_dims!=NULL) {
            break;
        }
    }

    if (parallel_dims != NULL) {
        printf("Parallel dims\n");
        for (i=0; i<nvar; i++) {
            printf("%d,", parallel_dims[i]);
        }
        printf("\n");
    }
    /* exit(0); */
    if (parallel_dims == NULL) {
        printf("No parallel dims in scc %d\n", scc_id);
        is_successful = colour_scc(scc_id, colour,c,0, -1, prog);
        return is_successful;
    } else {
        for (j=0; j<prog->nvar; j++) {
            if (parallel_dims[i]==1) {

            }
            /* if (parallel_colour_scc_dimension(scc_id, colour, c, prog)) */
                break;
        }
    }
/* TODO: Fix this return statement after incorporating clustering heuristics */
   return false;

}
#endif

/* This cuts disconnects the SCC from the DDG by cutting all
 * the incoming and outgoing edges from the given scc */
void cut_around_scc (int scc_id, PlutoProg *prog)
{
    int j;
    Graph *ddg;
    ddg = prog->ddg;
    for (j=0; j<ddg->num_sccs; j++) {

        if (scc_id == j) continue;

        if ((j < scc_id) && ddg_sccs_direct_connected(ddg,prog,j,scc_id)) {
            IF_DEBUG(printf("[colour SCC]: Cutting between scc %d and %d\n",j,scc_id););
            if (options->fuse == NO_FUSE) {
                cut_all_sccs(prog,ddg);
            } else {
                cut_between_sccs(prog,ddg,j,scc_id);

                /* You also need to cut between a successor node as well */
                for (j=scc_id+1; j<ddg->num_sccs; j++) {
                    if (ddg_sccs_direct_connected(ddg, prog, scc_id, j)) {
                        IF_DEBUG(printf("[colour SCC]: Cutting between scc %d and %d\n",scc_id,j););
                        cut_between_sccs(prog, ddg, scc_id, j);
                        /* cut_all_sccs(prog, ddg); */
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
            break;
        }
    }
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
    bool is_parallel = false;

    /* ToDo: Check if this condition can really happen.  */
    if (stmt_pos >= sccs[scc_id].size){
        return true;
    }

    if (prog->coloured_dims >= sccs[scc_id].max_dim) {
        if (prog->coloured_dims > sccs[scc_id].max_dim){
            return true;
        }
        IF_DEBUG(printf("[colour SCC]: All Dimensions of statment %d in SCC %d have been coloured\n", sccs[scc_id].vertices[stmt_pos], scc_id););
        /* cut if the scc's are not already distributed and you can not colour further. */
        /* for each SCC which is greater than the current scc, if there is a dep edge 
         * between these scc's then cut between these scc's. The cut has to respect 
         * the existing dependence */

        /* This is just for experimental purposes. The following if code can 
         * be removed if the assert never fails */
        if (sccs[scc_id].size !=1 ){
            printf("SCC %d has size %d\n", scc_id,sccs[scc_id].size);
            int i;
            for (i=0;i<sccs[scc_id].size; i++){
                printf("S%d,",sccs[scc_id].vertices[i]);
            }
            printf("\n");
        }
        assert (sccs[scc_id].size ==1);

        if (sccs[scc_id].size == 1) {
            for (j=0;j<ddg->num_sccs; j++) {
                if (scc_id!=j) {
                    if ((j < scc_id) && ddg_sccs_direct_connected(ddg,prog,j,scc_id)) {
                        IF_DEBUG(printf("[colour SCC]: Cutting between scc %d and %d\n",j,scc_id););
                        if (options->fuse == NO_FUSE) { 
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

    while (num_discarded!=nvar) {
        j = get_next_min_vertex(fcg_offset, stmt_id, list, num_discarded, pv, prog);
        if (options->scc_cluster) {
            IF_DEBUG(printf("[Colour SCC] Trying Colouring dimension %d of scc %d with colour %d\n",j,scc_id,c););
        } else {
            IF_DEBUG(printf("[Colour SCC] Trying Colouring dimension %d of statement %d with colour %d\n",j,stmt_id,c););
        }

        v = fcg_offset+j;

        /* If the dimension is already coloured with a different colour. 
         * Else it tries to check if the existing colour is fine. This is done 
         * as opposed to undoing the existing colour and then redoing it in 
         * the next step once FCG is rebuilt */
        if (colour[v]>0 && colour[v]!=c) {
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

        /* This check is redundant. covered in the next condition; */
        if (pv>=0 && is_adjecent(fcg,v,pv)) {
            list[num_discarded] = v;
            num_discarded++;
            continue;
        }

        /* Check if this is a valid colour */
        if (is_valid_colour(v,c,fcg,colour, is_parallel)) {
            colour[v] = c;
            /* If this is a valid colour, then try colouring the next vertex in the SCC */
            if (colour_scc(scc_id, colour, c, stmt_pos+1, v, prog)) {
                IF_DEBUG(printf("[Colour SCC] Colouring dimension %d",j););
                IF_DEBUG(printf("of statement %d with colour %d\n",stmt_id,c););
                return true;
            } else { 
                list[num_discarded] = v;
                num_discarded++;
                IF_DEBUG(printf("[Colour SCC] Unable to colour dimension %d of",j););
                IF_DEBUG(printf("statement %d with colour %d\n",stmt_id,c););
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

/* Returns an array corresponding to convex
 * successors that have parallel dimensions */
int* get_convex_successors(int scc_id, PlutoProg *prog,
        int* num_convex_successors)
{
    int i, num_successors, num_sccs;
    Graph *ddg;
    Scc *sccs;
    int* convex_successors = NULL;
    ddg = prog->ddg;
    sccs = ddg->sccs;
    num_sccs = ddg->num_sccs;

    num_successors=0;
    for (i=scc_id+1; i<ddg->num_sccs; i++) {
        if (sccs[i].is_parallel && is_convex_scc(scc_id, i, ddg, prog)) {
            if (convex_successors == NULL) {
                convex_successors = (int*) malloc (num_sccs*sizeof(int));
            }
            convex_successors[num_successors++] = i;
        }
    }
    *num_convex_successors = num_successors;
    return convex_successors;
}

/* Returns true if there is a parallelism preventing edge between vertex v
 * and any vertex i which has been coloured with the current_colour */
bool is_colour_par_preventing (int v, int *colour, int current_colour)
{
    int i;
    for (i=0; i<par_preventing_adj_mat->nrows; i++) {
        if (colour[i]==current_colour && (par_preventing_adj_mat->val[v][i]
                    || par_preventing_adj_mat->val[i][v])) {
            return true;
        }
    }
    return false;
}

/* Returns the common parallel dims of other scc that 
 * can be coloured with a dimension of the current scc. 
 * A dimension is colourable if it is parallel and does 
 * not have a self edge in the FCG */
int* get_common_parallel_dims(int scc_id, int* convex_successors,
        int num_convex_successors, int *colour, int current_colour,
        bool *is_colourable, PlutoProg *prog)
{
    int i, j, k, v;
    Graph *ddg, *fcg;
    Scc *sccs;
    int scc_offset, succ_scc_offset, succ_scc;
    int *common_dims;
    bool is_parallel = true;

    ddg= prog->ddg;
    fcg = prog->fcg;
    sccs = ddg->sccs;

    common_dims = NULL;
    scc_offset = sccs[scc_id].fcg_scc_offset;
    for (k=0; k<sccs[scc_id].max_dim; k++) {
        if (colour[scc_offset+k]!=0 || !is_colourable[k]) continue;
        for (i=0; i<num_convex_successors; i++) {
            succ_scc = convex_successors[i];
            succ_scc_offset = sccs[succ_scc].fcg_scc_offset;
            for (j=0; j<sccs[succ_scc].max_dim; j++) {
                v = succ_scc_offset+j;
                /* The vertex can be coloured if there is no self edge on j, 
                 * no edge between dimension j and k in the fcg, and there is 
                 * no vertex adjecent to j that is already coloured. Also 
                 * vertex j must be parallel and fusing with dimension k 
                 * must not hinder parallelism */
                if (colour[v]==0 && !fcg->adj->val[v][v] && ! is_adjecent(fcg, v, scc_offset+k) 
                        && is_valid_colour (v, current_colour, fcg, colour, is_parallel) &&
                        !par_preventing_adj_mat->val[v][v] && 
                        !(par_preventing_adj_mat->val[v][scc_offset+k] || 
                            par_preventing_adj_mat->val[scc_offset+k][v])) {
                    if (common_dims==NULL) {
                        common_dims = (int*) malloc (sizeof(int)*sccs[scc_id].max_dim);
                        bzero(common_dims, sccs[scc_id].max_dim*sizeof(int));
                    }
                    /* If you have found one dimension in a convex successor then 
                     * check the next convex successor. Hence the common_dims array 
                     * indicates the number of successors that can be fused with the 
                     * kth dimension of scc scc_id */
                    common_dims[k] ++;
                    break;
                }
            }
        }
    }
    return common_dims;
}

/* Returns the dimension which has the maximum number of 
 * convex successors. If there are none, then it returns -1  */
int get_colouring_dim(int *common_dims, int max_dim)
{
    int i, dim, max;
    max = 0;
    dim = -1;
    if (common_dims == NULL) {
        return dim;
    }
    for (i=0;i<max_dim; i++) {
        if (common_dims[i] > max) {
            max = common_dims[i];
            dim = i;
        }
    }
    return dim;
}


/* colours sccs listed in convex_successors with the colour 
 * current_colour. k is the vertex that was coloured previously */
void colour_convex_successors(int k, int *convex_successors, int num_successors,
        int *colour, int current_colour, PlutoProg *prog)
{
    Graph *fcg;
    Scc *sccs;
    int i, j, v, max_dim, scc_offset, scc_id;
    bool colourable_successor, is_parallel;

    fcg = prog->fcg;
    sccs = prog->ddg->sccs;
    is_parallel = true;

    for (i=0; i<num_successors; i++) {
        scc_id = convex_successors[i];
        max_dim = sccs[scc_id].max_dim;
        scc_offset = sccs[scc_id].fcg_scc_offset;
        for (j=0;j<max_dim; j++) {
            v = scc_offset +j;
            colourable_successor = colour[v]==0 && !fcg->adj->val[v][v] &&
                !par_preventing_adj_mat->val[v][k] &&
                !par_preventing_adj_mat->val[v][v] &&
                is_valid_colour (v, current_colour, fcg, colour, is_parallel);
            if (!colourable_successor) continue;

            IF_DEBUG(printf("[pluto] Colouring dimension %d of Scc", j););
            IF_DEBUG(printf("%d with colour %d\n", scc_id, current_colour););
            colour[v] = current_colour;
            sccs[scc_id].is_scc_coloured = true;
            sccs[scc_id].has_parallel_hyperplane = true;
            break;
        }
    }
}

/* This routine implements the greedy clustering heuristic in typed fuse.
 * The scc being coloured has atleast one  parallel dimension. It looks at the
 * successors that are convex with the current SCC (scc_id) and finds the
 * parallel dimensions that are common with dimensions of the current SCC
 * and can be coloured. The dimension of the SCC scc_id with the maximum
 * common successors is chosen for colouring */
bool colour_scc_cluster_greedy(int scc_id, int *colour, int current_colour,
       PlutoProg* prog)
{
    Graph *ddg, *fcg;
    Scc *sccs;
    int i,v,max_dim, num_convex_successors, num_parallel_dims;
    bool* parallel_dims, is_dim_parallel, is_parallel;
    int *convex_successors, *common_dims, colouring_dim;

    ddg = prog->ddg;
    fcg = prog->fcg;
    sccs = ddg->sccs;

    v = sccs[scc_id].fcg_scc_offset;
    max_dim = sccs[scc_id].max_dim;
    num_convex_successors = 0;
    common_dims = NULL;

    /* Get parallel dimensions of the input scc */
    parallel_dims = (bool*) malloc(max_dim*sizeof(bool));
    bzero (parallel_dims, max_dim*sizeof(bool));

    is_parallel = true;

    num_parallel_dims = 0;
    for (i=0; i<max_dim; i++) {
        is_dim_parallel = (colour[v+i]==0 && !fcg->adj->val[v+i][v+i] &&
                !par_preventing_adj_mat->val[v+i][v+i] &&
                is_valid_colour(v+i, current_colour, fcg, colour, is_parallel));
        if (is_dim_parallel) {
            parallel_dims[i] = true;
            num_parallel_dims ++;
        }
    }

    if (num_parallel_dims == 0) {
        free (parallel_dims);
        return false;
    }

    convex_successors = get_convex_successors(scc_id, prog,
            &num_convex_successors);

    /* If there are no convex successors, colour this scc */
    if (num_convex_successors == 0) {
        for (i=0; i<max_dim; i++) {
            if (parallel_dims[i]) {
                colour[v+i] = current_colour;
                sccs[scc_id].is_scc_coloured = true;
                sccs[scc_id].has_parallel_hyperplane = true;
                IF_DEBUG(printf("Dimension %d of SCC %d ", i, scc_id););
                IF_DEBUG(printf("coloured with colour %d\n", current_colour););
                free (parallel_dims);
                return true;
            }
        }
    }

    common_dims = get_common_parallel_dims(scc_id, convex_successors, 
            num_convex_successors, colour, current_colour, parallel_dims, prog);

    colouring_dim = get_colouring_dim(common_dims,max_dim);

    if (common_dims != NULL) {
        free(common_dims);
    }

    if (colouring_dim == -1) {
        for (i=0; i<max_dim; i++) {
            if (parallel_dims[i]) {
                colour[v+i] = current_colour;
                sccs[scc_id].is_scc_coloured = true;
                sccs[scc_id].has_parallel_hyperplane = true;
                IF_DEBUG(printf("Dimension %d of SCC %d ", i, scc_id););
                IF_DEBUG(printf("coloured with colour %d\n", current_colour););
                free (parallel_dims);
                free (convex_successors);
                return true;
            }
        }
    }
    colour[v+colouring_dim] = current_colour;
    sccs[scc_id].is_scc_coloured = true;
    sccs[scc_id].has_parallel_hyperplane = true;
    IF_DEBUG(printf("Dimension %d of SCC %d ", colouring_dim, scc_id););
    IF_DEBUG(printf("coloured with colour %d\n", current_colour););

    colour_convex_successors(v+colouring_dim, convex_successors,
            num_convex_successors, colour, current_colour, prog);

    free(convex_successors);
    free(parallel_dims);
    return true;
}

/* Cuts the SCC given by scc_id from its immediate predecessor. */
void cut_from_predecessor(int scc_id, PlutoProg* prog)
{
    int i;
    Graph *ddg;
    ddg = prog->ddg;
    for (i=scc_id-1; i>=0; i--) {
        if (ddg_sccs_direct_connected(ddg, prog, i, scc_id)) {
            cut_between_sccs(prog, ddg, i, scc_id);
        }
    }
}

/* Returns a dimension of an SCC (vetex in the FCG)
 * that can be colored next. */
int get_next_min_vertex_scc_cluster(int scc_id, PlutoProg *prog,
        int num_discarded, int* discarded_list)
{
    int i, v, max_dim, scc_offset;
    Scc *sccs;

    sccs = prog->ddg->sccs;
    max_dim = sccs[scc_id].max_dim;
    scc_offset = sccs[scc_id].fcg_scc_offset;
    for (i=0; i<max_dim; i++) {
        v = scc_offset + i;
        if (!is_discarded(v, discarded_list, num_discarded)) {
            return v;
        }
    }
    /* If none of the vertices are colourable */
    return -1;
}

/* Colours an SCC of the FCG in the clustered approach */
bool colour_scc_cluster (int scc_id, int *colour,
        int current_colour, PlutoProg* prog)
{
    int max_dim, scc_offset;
    int i, v;
    Graph *ddg, *fcg;
    Scc *sccs;
    int *disc_list, num_discarded;

    bool check_parallel;
    ddg = prog->ddg;
    fcg = prog->fcg;
    sccs = ddg->sccs;

    max_dim = ddg->sccs[scc_id].max_dim;

    /* All dimensions of the current SCC have already been coloured */
    if (prog->coloured_dims > max_dim) return true;

    if (prog->coloured_dims == max_dim) {
        if (!options->delayed_cut) {
            cut_around_scc (scc_id, prog);
        }
        return true;
    }

    bool hybrid_cut = options->hybridcut && sccs[scc_id].has_parallel_hyperplane;
    /* If the SCC has a parallel hyperplane and the fusion strategy is hybrid, 
     * then look max_fuse instead of greedy typed fuse heuristic */
    if (options->fuse == TYPED_FUSE &&
            sccs[scc_id].is_parallel && !hybrid_cut ) {
        if (colour_scc_cluster_greedy(scc_id, colour, current_colour, prog)) {
            sccs[scc_id].has_parallel_hyperplane = true;
            return true;
        } else {
            /* Colouring might fail because of a parallelism preventing egde. 
             * Hence cut between SCCs and try colouring again. */
            cut_from_predecessor(scc_id, prog);
            /* pluto_print_colours(colour, prog); */
            IF_DEBUG(printf("Updating FCG between SCCs %d and %d\n", scc_id, scc_id-1););
            update_fcg_between_sccs (fcg, scc_id-1, scc_id, prog);

            if (colour_scc_cluster_greedy(scc_id, colour, current_colour, prog)) {
                sccs[scc_id].has_parallel_hyperplane = true;
                return true;
            }
            /* Colouring the cluster has failed due to a permute 
             * or fusion preventing edge. Hence FCG has to be rebuilt. */
            IF_DEBUG(printf ("Unable to colour scc %d\n", scc_id););
            return false;
        }
    } 

    scc_offset = sccs[scc_id].fcg_scc_offset;
    check_parallel = false;
    disc_list = (int*)malloc(max_dim*sizeof(int));
    num_discarded = 0;
    /* for (i =0; i<max_dim; i++) { */
    do {
        v = get_next_min_vertex_scc_cluster(scc_id, prog, num_discarded, disc_list);
        assert (v != -1);
        /* v = scc_offset + i; */
        /* Used for debugging purposes */
        i = scc_offset - max_dim;
        if (colour[v]>0 && colour[v]!=current_colour) {
            IF_DEBUG(printf("Dimension %d of SCC %d ", i, scc_id););
            IF_DEBUG(printf("already coloured with colour %d\n",colour[v]););
            disc_list[num_discarded] = v;
            num_discarded ++;
            continue;
        }

        /* Cannot colour a vertex with a self edge. This is will not be 
         * covered in the next case as the vertex v is not coloured in 
         * the first place */
        if (fcg->adj->val[v][v] == 1) {
            disc_list[num_discarded] = v;
            num_discarded ++;
            continue;
        }

        /* Check if there is an adjecent vertex with the same colour. 
         * In case of typed fuse it checks if there is an adjecent vertex 
         * with the same colour and has a parallelism preventing edge  */
        hybrid_cut = options->hybridcut && sccs[scc_id].has_parallel_hyperplane;
        if (is_valid_colour(v, current_colour, fcg, colour, check_parallel)) {
            if (!hybrid_cut && options->fuse == TYPED_FUSE && sccs[scc_id].is_parallel &&
                    is_colour_par_preventing(v, colour, current_colour)) {
                IF_DEBUG(printf("Cannot colour dimension %d of SCC %d\n", i, scc_id););
                disc_list[num_discarded] = v;
                num_discarded ++;
                continue;
            }
            colour[v] = current_colour; 
            IF_DEBUG(printf("Colouring dimension %d of SCC %d ", i, scc_id););
            IF_DEBUG(printf("with colour %d\n", colour[v]););
            return true;
        }

        disc_list[num_discarded] = v;
        num_discarded ++;
    } while (num_discarded != max_dim);
    return false;
}

/* Returns colours corresponding vertices of the original FCG 
 * from the colours of vertices of scc clustered FCG */
int* get_vertex_colour_from_scc_colour (PlutoProg *prog, int *colour,
        int* has_parallel_hyperplane)
{
    int i, j, scc_offset, scc_id;
    int nvar, nstmts;
    int *stmt_colour;
    Scc *sccs;
    Stmt **stmts;

    nvar = prog->nvar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;
    sccs = prog->ddg->sccs;

    stmt_colour = (int*) malloc (nstmts*(nvar)*sizeof(int));
    bzero(stmt_colour,nvar*nstmts*sizeof(int));
    for (i=0; i<nstmts; i++) {
        scc_id = stmts[i]->scc_id;
        scc_offset = sccs[scc_id].fcg_scc_offset;
        for (j=0; j<stmts[i]->dim_orig; j++) {
            stmt_colour[i*(nvar)+j] = colour[scc_offset+j];
        }
        if (options->fuse == TYPED_FUSE) {
            has_parallel_hyperplane[i] = sccs[scc_id].has_parallel_hyperplane;
        }
    }
    return stmt_colour;
}

/* Returns colours corresponding to clustered FCG from the colours of
 * the statement. The routine picks the colour of the statement whose
 * dimensionality is same as the dimensionality of the SCC as the colour
 * of the SCC */
int* get_scc_colours_from_vertex_colours (PlutoProg *prog, int *stmt_colour,
        int current_colour, int nvertices, int* has_parallel_hyperplane)
{
    int i, j, scc_offset, stmt_id;
    int nvar, num_sccs;
    int *scc_colour;
    Scc *sccs;
    Stmt **stmts;

    nvar = prog->nvar;
    stmts = prog->stmts;
    num_sccs = prog->ddg->num_sccs;
    sccs = prog->ddg->sccs;

    scc_colour = (int*) malloc (nvertices*sizeof(int));
    bzero(scc_colour, nvertices*sizeof(int));

    scc_offset = 0;
    stmt_id = -1;

    for (i=0; i<num_sccs; i++) {
        sccs[i].is_scc_coloured = false;
        for (j=0; j<sccs[i].size; j++) {
            stmt_id = sccs[i].vertices[j];
            if (sccs[i].max_dim == stmts[stmt_id]->dim) break;
        }
        assert (stmt_id >= 0);
        for (j=0; j<sccs[i].max_dim; j++) {
            if (stmt_colour[(stmt_id*nvar)+j] == current_colour) {
                sccs[i].is_scc_coloured = true;
            }
            scc_colour[scc_offset+j] = stmt_colour[stmt_id*(nvar)+j];
        }
        sccs[i].fcg_scc_offset = scc_offset;
        scc_offset += sccs[i].max_dim;
        if (options->fuse == TYPED_FUSE) {
            sccs[i].has_parallel_hyperplane = has_parallel_hyperplane[stmt_id];
        }
    }
    return scc_colour;
}

/* Routine to rebuild FCG in the clustered approach. 
 * As SCC ids may change when SCC's are recomputed, 
 * the colours for each statement is computed and then 
 * the colours for the updated sccs are obtained using statement colours. */
int* rebuild_scc_cluster_fcg (PlutoProg *prog, int *colour, int c)
{
    int *stmt_colour, nvertices, i, num_sccs;
    int *scc_colour;
    Graph *ddg, *fcg;
    int* has_parallel_hyperplane = NULL;

    fcg = prog->fcg;
    ddg = prog->ddg;


    if (options->fuse == TYPED_FUSE) {
        has_parallel_hyperplane = (int*)malloc((prog->nstmts)*sizeof(int));
    }
    stmt_colour = get_vertex_colour_from_scc_colour(prog, colour,
                                        has_parallel_hyperplane);
    free_scc_vertices(ddg);

    /* You can update the DDG but do not update the FCG.  Doing otherwise will remove 
     * edges wich prevents permutation which is unsound */
    ddg_update(ddg, prog);
    IF_DEBUG(printf("DDG after colouring with colour %d\n",c););
    IF_DEBUG(pluto_matrix_print(stdout, ddg->adj););
    ddg_compute_scc(prog);
    compute_scc_vertices(ddg);
    num_sccs = prog->ddg->num_sccs;

    nvertices = 0;

    for (i=0;i<num_sccs; i++) {
        nvertices += prog->ddg->sccs[i].max_dim;
    }

    scc_colour = get_scc_colours_from_vertex_colours (prog, stmt_colour, c,
            nvertices, has_parallel_hyperplane);
    if(options->fuse == TYPED_FUSE) {
        pluto_matrix_free(par_preventing_adj_mat);
    }
    prog->fcg = build_fusion_conflict_graph(prog, scc_colour, nvertices, c);

    /* These two have to be reset in the clustered apporoach as 
     * Scc's will change when FCG is rebuilt and 
     * they will be revisited during colouring */
    prog->fcg->num_coloured_vertices = 0;
    prog->total_coloured_stmts[c-1] = 0;
    prog->fcg->to_be_rebuilt = false;

    graph_free(fcg);
    free (colour);
    free(stmt_colour);
    free (has_parallel_hyperplane);
    return scc_colour;
}

/* The routine checks if the two SCCs are fused in the till the current level.
 * If yes it returns true else returns false */
bool are_sccs_fused(PlutoProg *prog, int scc1, int scc2)
{
    int i, num_hyperplanes, stmt1, stmt2, nvar, npar;
    Scc *sccs;
    Stmt **stmts;
    bool sccs_fused = true;
    int cut_id1, cut_id2;

    num_hyperplanes = prog->num_hyperplanes;
    sccs = prog->ddg->sccs;
    stmts = prog->stmts;
    nvar = prog->nvar;
    npar = prog->npar;

    stmt1 = sccs[scc1].vertices[0];
    stmt2 = sccs[scc2].vertices[0];

    for (i=0; i<num_hyperplanes; i++) {
        if (prog->hProps[i].type == H_LOOP) {
            continue;
        }

        cut_id1 = stmts[stmt1]->trans->val[i][nvar+npar];
        cut_id2 = stmts[stmt2]->trans->val[i][nvar+npar];

        if (cut_id1 != cut_id2) {
            sccs_fused = false;
            break;
        }
    }
    return sccs_fused;
}

/* Routine adds a scalar hyperplane between two SCCs.
 * Note that SCCs might not be connected. Ideally this routine
 * can be moved to pluto.c (or framework.c) and approprite 
 * changes can be made in DDG cut routines */
void pluto_add_scalar_hyperplanes_between_sccs(PlutoProg *prog, int scc1, int scc2)
{
    int i, j, nstmts, nvar, npar;
    Stmt **stmts;

    stmts = prog->stmts;
    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;

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
}

/* Adds a hyperplane (H_LOOP) from the solution an ILP solution. 
 * This method can also be moved to framework.c (or pluto.c) and 
 * appropriate changes can be made. */
void add_hyperplane_from_ilp_solution (int64 *sol, PlutoProg *prog)
{
    int j, k;
    int nvar, npar,nstmts;
    Stmt **stmts;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;

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
}

/* Colours all scc's with a colour c. Returns the current colouring of the fcg. */
int* colour_fcg_scc_based(int c, int *colour, PlutoProg *prog)
{
    int i,j,nsccs,prev_scc;
    bool is_distributed, is_successful;
    Graph *ddg, *fcg;
    double t_start;
    bool is_parallel_scc_coloured, hybrid_cut;

    ddg = prog->ddg;
    fcg = prog->fcg;
    nsccs = ddg->num_sccs;

    is_distributed = false;
    is_parallel_scc_coloured = false;
    prev_scc = -1;

    for (i=0; i<nsccs; i++) {
        t_start = rtclock();

        /* In clustering approach, when FCG is rebuilt, DDG is upadated.
         * However, if some Sccs were previously coloured before the rebuilding
         * the FCG, we dont have to re-colour those SCCs again. If fcg has to
         * be rebuilt, then SCC ids would not have changed from previous
         * colouring */
        if (options->scc_cluster && prog->ddg->sccs[i].is_scc_coloured) {
            fcg->num_coloured_vertices += ddg->sccs[i].max_dim;
            prog->total_coloured_stmts[c-1] += ddg->sccs[i].size;
            prev_scc = i;
            prog->fcg_colour_time += rtclock() - t_start;
            /* pluto_print_colours(colour, prog); */
            continue;
        }

        IF_DEBUG(printf("Trying Colouring Scc %d of Size %d ",i,ddg->sccs[i].size););
        IF_DEBUG(printf("with colour %d\n", c););
        if (options->scc_cluster) {
            hybrid_cut = options->hybridcut && ddg->sccs[i].has_parallel_hyperplane;
            if (options->fuse == TYPED_FUSE) {
                /* case when the previous scc that was coloured was parallel
                 * and the current one is seqential */
                if (!ddg->sccs[i].is_parallel && is_parallel_scc_coloured
                        && prev_scc != -1) {
                    if (are_sccs_fused(prog, prev_scc, i)) {
                        /* distribute the loops here. Note that
                         * sccs may not be connected at all. However we still
                         * need to cut to preserve parallelism */
                        if (ddg_sccs_direct_connected(ddg, prog, prev_scc, i)) {
                            cut_between_sccs(prog, ddg, prev_scc, i);
                            update_fcg_between_sccs(fcg, prev_scc, i, prog);
                        } else {
                            pluto_add_scalar_hyperplanes_between_sccs(prog, prev_scc, i);
                        }
                    }
                } else if (ddg->sccs[i].is_parallel && !hybrid_cut 
                        && !is_parallel_scc_coloured && prev_scc != -1) {
                    if (are_sccs_fused(prog, prev_scc, i)) {
                        /* distribute the loops here. Note that
                         * sccs may not be connected at all. However we still
                         * need to cut to preserve parallelism */
                        if (ddg_sccs_direct_connected(ddg, prog, prev_scc, i)) {
                            cut_between_sccs(prog, ddg, prev_scc, i);
                            update_fcg_between_sccs(fcg, prev_scc, i, prog);
                        } else {
                            IF_DEBUG(printf("Adding Scalar hyperplanes without cut\n"););
                            pluto_add_scalar_hyperplanes_between_sccs(prog, prev_scc, i);
                        }
                    }
                }
                /* Set that a parallel SCC is being coloured */
                if (ddg->sccs[i].is_parallel) {
                    is_parallel_scc_coloured = true;
                } else {
                    is_parallel_scc_coloured = false;
                }
            }
            is_successful = colour_scc_cluster (i, colour, c, prog);
        /* } else if (options->fuse ==TYPED_FUSE && ddg->sccs[i].is_parallel) { */
            /* is_successful = colour_scc_from_lp_solution_with_parallelism(i, colour, prog, c); */
        } else {
            is_successful = colour_scc(i, colour, c, 0, -1, prog);
        }
        /* If colouring fails in the fist SCC */
        if (!is_successful) {
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
                if (options->scc_cluster) {
                    colour = rebuild_scc_cluster_fcg (prog, colour,c);
                    /* rebuliding the cluster_fcg will update ddg, hence number of sccs can increase */
                    nsccs = prog->ddg->num_sccs;
                    fcg = prog->fcg;

                    /* Sccs will be renumbered; hence all sccs have to be revisited; */
                    i=-1;
                    prev_scc = -1;
                    continue;
                } else {
                    prog->fcg = build_fusion_conflict_graph(prog, colour, fcg->nVertices, c);
                }

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
                    if (options->scc_cluster) {
                        is_distributed = colour_scc_cluster(i, colour, c, prog);
                    } else {
                        is_distributed = colour_scc(i, colour, c, 0, -1, prog);
                    }
                    if (!is_distributed){
                        /* Colouring was prevented by a fusion preventing dependence. 
                         * Therefore cut DDG then update FCG and then colour */
                        IF_DEBUG(printf("FCG Before Updating\n"););
                        IF_DEBUG(pluto_matrix_print(stdout, fcg->adj););

                        if (options->fuse == NO_FUSE) {
                            cut_all_sccs(prog,ddg);
                            /* TODO: Update this call testing is done */
                            update_fcg_between_sccs(fcg, 0, 0, prog);
                        } else {
                            for (j=prev_scc; j>=0; j--) {
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
                        if (options->scc_cluster) {
                            is_distributed = colour_scc_cluster(i, colour, c, prog);
                        } else {
                            is_distributed = colour_scc(i,colour,c,0, -1, prog);
                        }
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
                    for (j=prev_scc; j>=0; j--) {
                        if (ddg_sccs_direct_connected(ddg, prog, j, i)) {
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
                if (options->scc_cluster) {
                    is_distributed = colour_scc_cluster (i, colour, c, prog);
                } else {
                    is_distributed = colour_scc(i,colour,c,0, -1, prog);
                }
            }

            /* Needed in case of partial satisfaction */
            if (is_distributed == false) {
                printf("Num Deps satisfied with precise check %d\n",
                        pluto_compute_dep_satisfaction_precise(prog));

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
                if (options->scc_cluster) {
                    is_distributed = colour_scc_cluster (i, colour, c, prog);
                } else {
                    is_distributed = colour_scc(i,colour,c,0, -1, prog);
                }

            }
            assert (is_distributed == true);
        }

        prog->ddg->sccs[i].is_scc_coloured = true;
        if (options->scc_cluster) {
            fcg->num_coloured_vertices += ddg->sccs[i].max_dim;
        } else {
            fcg->num_coloured_vertices += ddg->sccs[i].size;
        }
        prog->total_coloured_stmts[c-1] += ddg->sccs[i].size;
        prev_scc = i;
        prog->fcg_colour_time += rtclock() - t_start;
    }
    return colour;
}

/* Resets that an scc has coloured. This is called after permutation
 * found at the current level is scaled and shifted, thus enabling
 * colouring at the next level. */
void reset_scc_colouring (Graph *ddg)
{
    int nsccs, i;
    nsccs = ddg->num_sccs;

    for(i=0; i<nsccs; i++) {
        ddg->sccs[i].is_scc_coloured = false;
    }
    return;
}

/* Routine to find permutable hyperplanes in the FCG based approach. 
 * Colouring is done on a per SCC basis. Natural number ordering of the 
 * scc ids is used to ensure convexity. */
void find_permutable_dimensions_scc_based(int *colour, PlutoProg *prog)
{
    int i,j,num_coloured_dims,max_colours;
    Stmt **stmts;
    Graph *ddg;
    double t_start;

    max_colours = prog->nvar;
    stmts = prog->stmts;

    if (options->fuse == TYPED_FUSE && options->scc_cluster) {
        for (i=0;i<prog->ddg->num_sccs; i++) {
            prog->ddg->sccs[i].has_parallel_hyperplane = 0;
        }
    }

    for (i=1; i<=max_colours; i++) {
        for (j=0; j<prog->ddg->num_sccs; j++) {
            prog->ddg->sccs[j].is_scc_coloured = false;
        }
        if (options->lpcolour && !options->scc_cluster) {
            mark_parallel_sccs(colour, prog);
        }
        colour = colour_fcg_scc_based(i, colour, prog);

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

        ddg = prog->ddg;

        if (options->lpcolour) {
            for (j=0; j<ddg->num_sccs; j++) {
                if (ddg->sccs[j].sol != NULL) {
                    free(ddg->sccs[j].sol);
                    ddg->sccs[j].sol = NULL;
                }
            }
        }

        /* Do not update ddg or sccs if sccs are clustered. It will be updated when FCG is rebuilt */
        if (!options->scc_cluster) {
            free_scc_vertices(ddg);

            /* You can update the DDG but do not update the FCG.  Doing otherwise will remove 
             * edges wich prevents permutation which is unsound */
            ddg_update(ddg, prog);
            IF_DEBUG(printf("DDG after colouring with colour %d\n",i););
            IF_DEBUG(pluto_matrix_print(stdout, ddg->adj););
            IF_DEBUG(printf("[Find_permutable_dims_scc_based]: Updating SCCs \n"););
            ddg_compute_scc(prog);
            compute_scc_vertices(ddg);
        }else{
            /* Reset SCC colouring to enable colouring at the next level*/
            reset_scc_colouring(ddg);
        }

        IF_DEBUG2(pluto_transformations_pretty_print(prog););
        IF_DEBUG2(pluto_compute_dep_directions(prog););
        IF_DEBUG2(pluto_compute_dep_satisfaction(prog););
        IF_DEBUG2(pluto_print_dep_directions(prog););
    }
    /* This can be delayed further. Introduce this cut after diamond tiling */

#if 0
    /* All dimensions have been coloured but still there are some deps that 
     * need to be satisfied at the innermost level by distribution.  */
    if (i == max_colours+1 && !deps_satisfaction_check(prog)) {
        cut_all_sccs(prog, prog->ddg);
    }
#endif

    IF_DEBUG(printf("[Pluto] Colouring Successful\n"););
    IF_DEBUG(pluto_print_colours(colour,prog););

    free(colour);
    if (options->fuse == TYPED_FUSE) {
        pluto_matrix_free(par_preventing_adj_mat);
    }
    if (options->lpcolour) {
        pluto_matrix_free(dep_distance_mat);
    }
    return;
}

/*************************** Scaling Routines ******************/
/* Set the lower bounds for statements that are coloured with colour c. 
 * These are used to find the scaling and shifting factors. 
 * When vertices of FCG are clustered lower bound of dimensions of all 
 * statements in the SCC is set to 1 if and only if the corresponding 
 * vertex of the FCG is coloured with colour c */
void add_coeff_constraints_from_scc_clustered_fcg_colouring (PlutoConstraints *coeffcst, 
        int *colour, int c, PlutoProg *prog) 
{
    int i, j, k, num_sccs, stmt_id, nvar ,npar, scc_offset;
    Stmt **stmts;
    Graph *ddg;
    Scc *sccs;

    nvar = prog->nvar;
    npar = prog->npar;
    ddg = prog->ddg;
    num_sccs = ddg->num_sccs;
    sccs = ddg->sccs;
    stmts = prog->stmts;
    scc_offset = 0;

    for (j=0; j<num_sccs; j++) {
        /* If the Scc's dimesionality is greater less c, then all vertices 
         * of this scc are coloured. Hence set the lower bound of the 
         * transformation coefficients of this SCC to zero.*/
        if (sccs[j].max_dim < c) {
            for (i=0; i<sccs[j].size; i++) {
                stmt_id = sccs[j].vertices[i];
                for(k=0; k<sccs[j].max_dim; k++) {
                    if (stmts[stmt_id]->is_orig_loop[k]) {
                        pluto_constraints_add_lb(coeffcst, 
                                npar+1+stmt_id*(nvar+1)+k, 0);
                    } else {
                        pluto_constraints_add_equality(coeffcst);
                        coeffcst->val[coeffcst->nrows-1]
                            [npar+1+stmt_id*(nvar+1)+k] = 1;
                    }
                }
            }
            scc_offset += sccs[j].max_dim;
            continue;
        }


        for (i=0; i<sccs[j].size; i++) {
            stmt_id = sccs[j].vertices[i];
            for (k=0; k<sccs[j].max_dim; k++) {
                if (colour[scc_offset+k]==c && stmts[stmt_id]->is_orig_loop[k]) {
                    pluto_constraints_add_lb(coeffcst,npar+1+stmt_id*(nvar+1)+k,1);
                } else {
                    pluto_constraints_add_equality(coeffcst);
                    coeffcst->val[coeffcst->nrows-1][npar+1+stmt_id*(nvar+1)+k] = 1;
                }
            }
        }
        scc_offset += sccs[j].max_dim;
    }
    return;
}

/* Set the lower bounds for statements that are coloured with colour c. 
 * These are used to find the scaling and shifting factors. */
void add_coeff_constraints_from_fcg_colouring (PlutoConstraints *coeffcst,
        int *colour, int c, PlutoProg *prog)
{
    int j, k, nstmts, nvar ,npar, stmt_offset;
    Stmt **stmts;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;
    stmt_offset = 0;

    for (j=0; j<nstmts; j++) {
        /* If a statements's dimesionality is less than c, then all 
         * dimensions of this statement are coloured. Hence set the
         *  lower bound of the transformation coefficients of this 
         *  statement to zero.*/
        if (stmts[j]->dim_orig < c) {
            for (k =0; k<stmts[j]->dim_orig; k++) {
                pluto_constraints_add_lb(coeffcst,npar+1+j*(nvar+1)+k, 1);
            }
            stmt_offset += stmts[j]->dim_orig;
            continue;
        }
        for (k=0; k<nvar; k++) {
            if (stmts[j]->is_orig_loop[k] && colour[stmt_offset+k] == c) {
                pluto_constraints_add_lb(coeffcst,npar+1+j*(nvar+1)+k, 1);
            } else {
                pluto_constraints_add_equality(coeffcst);
                coeffcst->val[coeffcst->nrows-1][npar+1+j*(nvar+1)+k] = 1;
            }
        }
        stmt_offset += stmts[j]->dim_orig;
    }
}

/* Once the permutation is found, it finds the scling and shifting factors for the permtation
 * Scales the dimensions in the with colour c+1. Returns 1 if scaling
 * was successful. Else returns 0. */
int scale_shift_permutations(PlutoProg *prog, int *colour, int c)
{
    int select;
    int nvar, npar, nrows;
    int nstmts;
    double t_start;
    PlutoConstraints *basecst,*coeffcst, *boundcst;
    int64 *sol;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;

    basecst = get_permutability_constraints(prog);
    assert (basecst->ncols == CST_WIDTH);

    boundcst = get_coeff_bounding_constraints(prog);
    nrows = basecst->nrows + boundcst->nrows + (nstmts*nvar);
    coeffcst = pluto_constraints_alloc(nrows, basecst->ncols);
    coeffcst->nrows = 0;
    coeffcst->ncols = basecst->ncols;
    assert (coeffcst->ncols == CST_WIDTH);

    IF_DEBUG(printf("Num stmts coloured with colour %d: %d\n", c+1,
                prog->total_coloured_stmts[c]););

    if (prog->total_coloured_stmts[c] != nstmts) {
        IF_DEBUG(printf("Not All statements have been coloured\n"););
        pluto_constraints_free(coeffcst);
        return 0;
    }

    coeffcst = pluto_constraints_copy(coeffcst,boundcst);
    pluto_constraints_free(boundcst);

    /* Pick a colour that you would start */
    select = c+1;
    IF_DEBUG(printf("[pluto] Finding Scaling factors for colour %d\n",select););

    /* Add CST_WIDTH number of cols and set appropriate constraints to 1 and set the rest to 0
     * These redundant cols are then removed. */

    if (options->scc_cluster) {
        add_coeff_constraints_from_scc_clustered_fcg_colouring(coeffcst, colour, select, prog);
    } else {
        add_coeff_constraints_from_fcg_colouring(coeffcst, colour, select, prog);
    }
    coeffcst = pluto_constraints_add(coeffcst,basecst);

    /* Solve the constraints to find the hyperplane at this level */
    t_start = rtclock();

    sol = pluto_prog_constraints_lexmin(coeffcst, prog);

    if (sol != NULL) {
        IF_DEBUG(printf("[pluto]: found a hyperplane\n"));
        /* num_sols_found++; */

        /* if (options->varliberalize) { */
        /*     for (j=0; j<ndeps; j++) { */
        /*         #<{(| Check if it has to be c or c+1 |)}># */
        /*         if (deps[j]->temp_across && c < deps[j]->fuse_depth  */
        /*                 && pluto_domain_equality(stmts[deps[j]->src],stmts[deps[j]->dest])) { */
        /*             for (k=0;k<nvar;k++) { */
        /*                 if (sol[npar+1+(deps[j]->src)*(nvar+1)+k] != sol[npar+1+(deps[j]->dest)*(nvar+1)+k]) */
        /*                     break; */
        /*             } */
        /*             if (k!=nvar || (sol[npar+1+(deps[j]->src)*(nvar+1)+nvar]!=sol[npar+1+(deps[j]->dest)*(nvar+1)+nvar])) { */
        /*                 printf("Cutting between SCCs to prevent illegal transformation with var-lib"); */
        /*                 cut_between_sccs(prog,ddg,stmts[deps[j]->src]->scc_id, stmts[deps[j]->dest]->scc_id); */
        /*             } */
        /*         } */
        /*     } */
        /* } */
        add_hyperplane_from_ilp_solution(sol, prog);
        prog->scaling_cst_sol_time += rtclock()-t_start;
        free(sol);
        IF_DEBUG(pluto_transformation_print_level(prog, prog->num_hyperplanes-1););
        pluto_constraints_free(coeffcst);
        return 1;
    } else {
        printf("[pluto] No Hyperplane found\n");
        if (options->delayed_cut) {
            coeffcst->nrows  = coeffcst->nrows - basecst->nrows;
            cut_smart(prog, prog->ddg);
            basecst = get_permutability_constraints(prog);
            coeffcst = pluto_constraints_add(coeffcst, basecst);
            sol = pluto_prog_constraints_lexmin(coeffcst, prog);
            if (sol != NULL) {
                add_hyperplane_from_ilp_solution(sol, prog);
                prog->scaling_cst_sol_time += rtclock()-t_start;
                free(sol);
                IF_DEBUG(pluto_transformation_print_level(prog,
                            prog->num_hyperplanes-1););
                pluto_constraints_free(coeffcst);
                return 1;
            }
        }
        pluto_constraints_free(coeffcst);
        prog->scaling_cst_sol_time += rtclock()-t_start;
        return 0;
    }
}

/* Routines that introduce loop skewing after loop permutations, loop skewing
 * and loop shifting transfomations have been found. */
bool get_negative_components(Dep *dep, bool *dims_with_neg_components,
       PlutoProg *prog, int level)
{
    int i;
    bool has_negative_comp;
    HyperplaneProperties *hProps;
    int loop_dims;

    hProps = prog->hProps;
    has_negative_comp = false;
    loop_dims = 0;
    for (i=0; i<prog->num_hyperplanes; i++){
        if (hProps[i].type == H_SCALAR && i < level){
            continue;
        }
        if (hProps[i].type == H_LOOP && i < level){
            loop_dims++;
            continue;
        }
        if (hProps[i].type == H_SCALAR && i >= level){
            continue;
        }
        if (dep->dirvec[i] == DEP_MINUS || dep->dirvec[i] == DEP_STAR ) {
            dims_with_neg_components[loop_dims] = 1;
            has_negative_comp = true;
            break;
        }
        loop_dims++;
    }
    return has_negative_comp;
}

/* Checks if there are any non constant dependences in the loop nest.
 * This is done by computing \vec{d}.\vec{h} and the resultant \vec{u}
 * has to be zero in pluto's cost function. Constraints to ensure u
 * is zero is added and if the resulting LP formulation is unsat, then
 * there are non constant deps in the SCC  */
bool constant_deps_in_scc(int scc_id, int level, PlutoConstraints *basecst,
       PlutoProg *prog) {
    int i, j, ndeps, nstmts, nvar, npar;
    int src_id, dest_id;
    Stmt **stmts;
    Dep **deps;

    ndeps = prog->ndeps;
    nstmts = prog->nstmts;
    npar = prog->npar;
    nvar = prog->nvar;
    deps = prog->deps;
    stmts = prog->stmts;

    basecst->nrows = 0;
    basecst->ncols = CST_WIDTH;
    for (i=0; i<ndeps; i++) {
        Dep *dep = deps[i];
        src_id = stmts[dep->src]->scc_id;
        dest_id = stmts[dep->dest]->scc_id;
        if ((src_id== scc_id) && (dest_id == scc_id) ) {
            pluto_constraints_add(basecst, dep->cst);
        }
    }
    /* If deps are constant then \vec{u} is zero vector.
     * Hence set the value of u to be zero */
    for (i=0;i<npar; i++){
        pluto_constraints_add_equality(basecst);
        basecst->val[basecst->nrows-1][i] = 1;
    }
    pluto_constraints_add_lb(basecst, npar, 0);
    /* Set the transformation coefficients to be equalities with
     * same values as in the input permutation */
    for (i=0; i<nstmts; i++) {
        for(j=0; j<nvar+1; j++) {
            pluto_constraints_add_equality(basecst);
            basecst->val[basecst->nrows-1][npar+1+i*(nvar+1)+j] = 1;
            basecst->val[basecst->nrows-1][basecst->ncols -1] = -stmts[i]->trans->val[level][j];
        }
    }

    IF_DEBUG(pluto_constraints_cplex_print(stdout, basecst););

    /* Replace this with LP call */
    int64* sol = pluto_prog_constraints_lexmin(basecst, prog);
    if (sol == NULL) {
        return false;
    }

    if(sol[npar]>10) {
        free(sol);
        return false;
    }
    free(sol);
    return true;
}

/* Returns a boolean array in which the vales that are set represent 
 * the dimensions of the current SCC that have to be skewed in order 
 * the make the loop nest tileable */
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
      
        if (dep_is_satisfied(dep)) {
            if (get_negative_components(dep,dims_with_neg_components, prog, level)){
                tile_preventing_deps[i] = 1;
            }
        }
    }
    return dims_with_neg_components;
}

/* Returns the dimensions at which tiling preventing dependences are satisfied */
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
        if (!tile_preventing_deps[i]) {
            continue;
        }
        for (j=0; j<prog->num_hyperplanes; j++) {
            if (j == dep->satisfaction_level) {
                break;
            }
            else if (hProps[j].type == H_LOOP) {
                loop_dims++;
            }
        }
        sat_dim[loop_dims] = 1;
    }
    return sat_dim;
}

/* Returns skewing constraints. The lower bounds of the 
 * dimensions that satisfy tiling preventing dependences are set to 1 */
PlutoConstraints *get_skewing_constraints(bool *src_dims, bool* skew_dims,
        int scc_id, PlutoProg* prog, int level, PlutoConstraints *skewCst)
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
        for (j=0; j<stmts[i]->dim_orig; j++){
            if (src_dims[j] && stmts[i]->scc_id == scc_id) {
                pluto_constraints_add_lb(skewCst, npar+1+ i*(nvar+1)+j, 1);
            }
            else {
                pluto_constraints_add_equality(skewCst);
                skewCst->val[skewCst->nrows-1][npar+1+i*(nvar+1)+j] = 1;
                /* Set the value of the current coeff to the one that you have already found */
                skewCst->val[skewCst->nrows-1][CST_WIDTH-1] = -stmts[i]->trans->val[level][j];
            }
        }
        pluto_constraints_add_lb(skewCst, npar+1+i*(nvar+1)+nvar, 0); 
    }
    return skewCst;
}

/*  Introduce loop skewing transformations if necessary.Returns true if
 *  skew was introuduced at some level for some SCC  */
bool introduce_skew(PlutoProg *prog)
{
    int i, j, k, num_sccs, nvar, npar, nstmts, level,ndeps;
    int initial_cuts, nrows, stmt_offset;
    Graph *orig_ddg;
    int64* sol;
    PlutoConstraints *skewingCst, *basecst, *const_dep_check_cst;
    HyperplaneProperties *hProps;
    Stmt **stmts;
    bool *src_dims, *skew_dims, tile_preventing_deps[prog->ndeps];
    double tstart;
    bool is_skew_introduced = false;

    nvar = prog->nvar;
    npar = prog->npar;
    nstmts = prog->nstmts;
    stmts = prog->stmts;
    ndeps = prog->ndeps;

    /* If there are zero or one hyperpane then you dont need to skew */
    if (prog->num_hyperplanes <=1) {
        return is_skew_introduced;
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

    pluto_dep_satisfaction_reset(prog);

    Graph *newDDG = ddg_create (prog);
    orig_ddg = prog->ddg;
    prog->ddg = newDDG;

    for (i=0; i<prog->ndeps; i++) {
        tile_preventing_deps[i] = 0;
    }

    initial_cuts = 0;
    for (level = 0; level< prog->num_hyperplanes; level++) {
        if (hProps[level].type == H_LOOP) {
            break;
        }
        initial_cuts ++; 
        dep_satisfaction_update(prog, level);
    }

    /* Needed to handle the case when there are no loops  */
    if (initial_cuts == prog->num_hyperplanes) {
        return is_skew_introduced;
    }
    basecst = get_permutability_constraints(prog);
    ddg_update(newDDG, prog);

    assert(level == initial_cuts);
    ddg_compute_scc(prog);
    num_sccs = newDDG->num_sccs;

    const_dep_check_cst = pluto_constraints_alloc(ndeps*nvar+1, CST_WIDTH);
    nrows = basecst->nrows+nstmts*(nvar+1);
    skewingCst = pluto_constraints_alloc(nrows, basecst->ncols);
    skewingCst->nrows = 0;
    skewingCst->ncols = CST_WIDTH;
    pluto_constraints_add(skewingCst, basecst);
    dep_satisfaction_update(prog,level);
    for (i =0; i<num_sccs; i++) {
        IF_DEBUG(printf("Looking for skews in SCC %d \n", i););
        /* dep_satisfaction_update(prog,level); */
        /* if (!constant_deps_in_scc(i, level, const_dep_check_cst, prog)) { */
        /*  */
        /*     IF_DEBUG(printf("Scc %d has atleast one non constant dep \n",i);); */
        /*     continue; */
        /* } */
        skew_dims = dims_to_be_skewed(prog, i, tile_preventing_deps, level);
        src_dims = innermost_dep_satisfaction_dims(prog, tile_preventing_deps);
        level++;

        for (; level <prog->num_hyperplanes; level++) {
            if (hProps[level].type != H_LOOP) {
                continue;
            }

            int skew_dim = 0;
            for (j=initial_cuts; j<prog->num_hyperplanes; j++) {
                if (prog->hProps[j].type == H_LOOP && skew_dims[skew_dim] == 1) {
                    level = j;
                    break;
                }
                else if (prog->hProps[j].type == H_LOOP){
                    skew_dim++;
                }
            }

            /* Skewing has to be done at level j+1 */
            if (j==prog->num_hyperplanes) {
                break;
            }

            skewingCst->nrows = basecst->nrows;
            get_skewing_constraints(src_dims, skew_dims, i, prog, level, skewingCst);

            sol = pluto_prog_constraints_lexmin(skewingCst, prog);

            if (!sol) {
                /* The loop nest is not tileable */
                free(skew_dims);
                break;
            }

            /* Set the Appropriate coeffs in the transformation matrix */
            for (j=0; j<nstmts; j++) {
                stmt_offset = npar+1+j*(nvar+1);
                for (k = 0; k<nvar; k++){
                    stmts[j]->trans->val[level][k] = sol[stmt_offset+k];
                }
                /* No parametric Shifts */
                for (k=nvar; k<nvar+npar; k++)    {
                    stmts[j]->trans->val[level][k] = 0;
                }
                /* The constant Shift */
                stmts[j]->trans->val[level][nvar+npar] = sol[stmt_offset+nvar];

            }
            free(sol);
            is_skew_introduced = true;

            dep_satisfaction_update(prog, level);
            pluto_compute_dep_directions(prog);

            free(skew_dims);

            if (level >= prog->num_hyperplanes-1) {
                continue;
            }

            skew_dims = dims_to_be_skewed(prog, i, tile_preventing_deps, level+1);
            free(src_dims);
            src_dims = innermost_dep_satisfaction_dims(prog, tile_preventing_deps);
        }
        free(src_dims);
        level = initial_cuts;
    }

    pluto_constraints_free(const_dep_check_cst);
    pluto_constraints_free(skewingCst);
    prog->ddg = orig_ddg;
    graph_free(newDDG);
    prog->skew_time += rtclock()-tstart;
    if (!options->silent) {
        printf("[Pluto]: Post processing skewing complete\n");
    }
    return is_skew_introduced;
}
#endif
