/*
 * PLUTO: An automatic parallelizer and locality optimizer
 *
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
#ifdef GUROBI
#include "constraints.h"
#include "pluto.h"
#include <math.h>
#include <assert.h>
#include <gurobi_c.h>

/* The file provides an interface to solve Pluto's constraints with Gurobi as 
 * the (I)LP solver. The file contains routines that creates a (I)LP problem
 * from Pluto's constraint matrix. The objective for the optimization problem 
 * is given as a PlutoMatrix. The constraints are solved a non-null solution 
 * is returned if a solution exists. The choice of solving ILP/LP can be 
 * given as a command line option to polycc.  */

/* Prints the error message corresponding to the error value passed */
inline void check_error_gurobi(GRBmodel *lp, int error)
{
    if(error) {
        GRBenv *env;
        env = GRBgetenv(lp);
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);

    }
}

/* Creates a double array for gurobi objective from the objective matrix.
 * The input PlutoMatrix should have a single row and the column elements 
 * correspond to the coefficients of the variables in the objective function. */
double* get_gurobi_objective_from_pluto_matrix(PlutoMatrix *obj)
{
    int j;
    double *grb_obj;
    assert (obj->nrows == 1);
    grb_obj = (double*) malloc(sizeof(double)*(obj->ncols));
    for (j=0; j<obj->ncols; j++) {
        grb_obj[j] = (double)(obj->val[0][j]);
    }
    return grb_obj;
}

/* Constructs constraints for gurobi problem from pluto_constraints.
 * Assumes that there are no rows or cols in the input model lp */
void set_gurobi_constraints_from_pluto_constraints (GRBmodel *lp, const PlutoConstraints *cst)
{
    int i, j, k, nrows, ncols;
    int *index;
    double *value, rhs;

    nrows = cst->nrows;
    ncols = cst->ncols-1;

    index = (int*) malloc ((ncols)*sizeof(int));
    value = (double*) malloc ((ncols)*sizeof(double));

    /*Populate the constraints row by row. */
    for (i=0; i<nrows; i++) {
        k = 0;
        for (j=0; j<ncols; j++) {
            if (cst->val[i][j] != 0) {
                index[k] = j;
                value[k] = (double) cst->val[i][j];
                k++;
            }
        }
        rhs = (double) -(cst->val[i][ncols]);
        /* Trivial constraint. can be skipped */
        if(k == 0) {
            continue;
        }
        if (cst->is_eq[i]) {
            /* For equality constraints the bounds are fixed */
            GRBaddconstr(lp, k, index, value, GRB_EQUAL, rhs, NULL);
        } else {
            /* Add a lower bound on each row of the inequality */
            GRBaddconstr(lp, k, index, value, GRB_GREATER_EQUAL, rhs, NULL);
        }
    }
    free(index);
    free(value);
}

/* Retrives ilp solution from the input gurobi problem. 
 * Assumes that the optimal solution exists and has been found*/
int64* get_ilp_solution_from_gurobi_problem(GRBmodel *lp)
{
    int num_cols, i, error;
    int64 *sol;
    double x;

    error = GRBgetintattr(lp, GRB_INT_ATTR_NUMVARS, &num_cols);
    check_error_gurobi(lp, error);

    sol = (int64*) malloc (num_cols*sizeof(int64));

    for (i=0; i<num_cols; i++) {
        error = GRBgetdblattrelement(lp, GRB_DBL_ATTR_X, i, &x);
        IF_DEBUG(printf("c%d = %lld, ", i, (int64) round(x)););
        sol[i] = (int64) round(x);

    }
    return sol;
}

/* Retrives lp solution from the input gurobi problem. 
 * Assumes that the optimal solution exists and has been found*/
double* get_lp_solution_from_gurobi_problem(GRBmodel *lp)
{
    int num_cols, error, i;
    double *sol;

    error = GRBgetintattr(lp, GRB_INT_ATTR_NUMVARS, &num_cols);
    check_error_gurobi(lp, error);

    sol = (double*)malloc(num_cols*sizeof(double));

    GRBgetdblattrarray(lp, GRB_DBL_ATTR_X, 0, num_cols, sol);

    if (options->debug) {
        for (i=0; i<num_cols; i++) {
            printf("c_%d=%0.6f,", i, sol[i]);
        }
        printf("\n");
    }

    return sol;
}

inline void find_optimal_solution_gurobi(GRBmodel *lp, double tol)
{
     GRBsetdblparam(GRBgetenv(lp), "IntFeasTol", tol);

     GRBoptimize(lp);
}

/* Solve the gurobi problem lp. If optimal solution is found, then it returns 0.
 * The caller can retrive the funtion from the gurobi model object lp. If the 
 * problem is infeasible then the routine returns 1. If the problem is 
 * unbounded, program terminates with the corresponding error message. */
int pluto_constraints_solve_gurobi(GRBmodel *lp, double tol)
{
    int optim_status;
    double objval;
    GRBenv *env = NULL;

     find_optimal_solution_gurobi(lp, tol);

     GRBgetintattr(lp, GRB_INT_ATTR_STATUS, &optim_status);

     if (optim_status == GRB_INFEASIBLE) {
         /* env = GRBgetenv(lp); */
         /* GRBfreemodel(lp); */
         /* GRBfreeenv(env); */
         return 1;
     }

     if (optim_status == GRB_UNBOUNDED) {
         env = GRBgetenv(lp);
         GRBfreemodel(lp);
         GRBfreeenv(env);
         printf("[Pluto]: Unbounded model. This appears to be a problem in Pluto's ILP/LP construction. Aborting\n");
         exit(1);
     }

     GRBgetdblattr(lp, GRB_DBL_ATTR_OBJVAL, &objval);
     IF_DEBUG(printf("Objective value: %0.6f\n",objval););

     return 0;
}

/* Returns the objective array for scaling MIP. The objective is to reduce the 
 * sum of the scaling factors of each connected component. */
double* get_gurobi_scaling_obj(int num_ccs, int num_sols_to_be_scaled)
{
    double *obj;
    int i;
    int nvars;
    nvars = num_sols_to_be_scaled + num_ccs;
    obj = (double*) malloc(nvars*sizeof(double));
    for (i=0; i<num_ccs; i++) {
        obj[i] = 1.0f;
    }
    for (i=num_ccs; i<nvars; i++) {
        obj[i] = 0.0f;
    }
    return obj;
}

/* Sets the lower bound of first lb1 variables to 1 
 * and the next lb0 variables to zero */
double* get_lower_bounds_for_variables(int lb1, int lb0)
{
    double *lb;
    int i;

    lb = (double*)malloc((lb0+lb1)*sizeof(double));

    /* Elements with lower bound one */
    for (i=0; i<lb1; i++) {
        lb[i] = 1.0f;
    }

    /* Elements with lower bound zero */
    for (i=lb1; i<lb1+lb0; i++) {
        lb[i] = 0.0f;
    }
    return lb;
}

/* Solves scaling MIP with an integer tolorence of 0.01 */
double* pluto_mip_scale_solutions_gurobi(GRBmodel *lp)
{
    double *scale_sols;
    int is_unsat;

    is_unsat = pluto_constraints_solve_gurobi(lp, 1e-2);
    if (is_unsat) {
        printf("[Pluto]: Error in scaling MIP, Aborting.\n");
        exit(1);
    }

    scale_sols = get_lp_solution_from_gurobi_problem(lp);
    return scale_sols;

}

/* Constructs a gurobi model for scaling MIP. Assumes that the CSR matrices are indexed from 1 instead of 0 */
GRBmodel *get_scaling_lp_gurobi(double *fpsol, int num_sols, double **val, int **index, int npar, int num_ccs, GRBenv *env)
{
    int i, num_sols_to_be_scaled, col_num;
    GRBmodel *lp = NULL;
    double *grb_obj, *lb;
    char *vtype;

    num_sols_to_be_scaled = num_sols-npar-1;

    vtype = (char*)malloc(num_sols_to_be_scaled + num_ccs);
    for (i=0; i<num_ccs; i++) {
        vtype[i] = GRB_CONTINUOUS;
    }
    for(i=num_ccs; i<num_ccs+num_sols_to_be_scaled; i++) {
        vtype[i] = GRB_INTEGER;
    }

    grb_obj = get_gurobi_scaling_obj(num_ccs, num_sols_to_be_scaled);

    lb = get_lower_bounds_for_variables(num_ccs, num_sols_to_be_scaled);

    /* GRBnewmodel(env, &lp, NULL, num_ccs+num_sols_to_be_scaled, grb_obj, NULL, NULL, vtype, NULL); */
    GRBnewmodel(env, &lp, NULL, num_ccs+num_sols_to_be_scaled, grb_obj, lb, NULL, vtype, NULL);

    for (i=npar+1; i<num_sols; i++) {
        col_num = i-npar-1+num_ccs;
        if (fabs(fpsol[i]) > 1e-7) {
            val[i-npar-1][1] = fpsol[i];

            /* This is specific to gurobi as it assumes that 
             * constraint matrix is indexed from 0 (as opposed to 1 in glpk) */
            index[i-npar-1][1] --;
            index[i-npar-1][2] --;
            int error = GRBaddconstr(lp, 2, &(index[i-npar-1][1]), &(val[i-npar-1][1]), GRB_EQUAL, 0.0, NULL);
            if(error) {
                printf("ERROR: %s\n", GRBgeterrormsg(env));
                exit(1);
            }
            GRBsetdblattrelement(lp, GRB_DBL_ATTR_LB ,col_num, 1.0);
        } else {
            /* Set the lower and upper bounds of the variable to zero */
            GRBsetdblattrelement(lp, GRB_DBL_ATTR_UB ,col_num, 0.0);
        }
    }

    return lp;
}

/* The rational solutions of pluto-lp are scaled on a per connected component basis. 
 * The following routine returns the maximum scaling factor among the scaling 
 * factors of each connected component; the largest among the first num_ccs 
 * elements of the array sol*/
int64 get_max_scale_factor(double *sol, int num_ccs)
{
    int i, max;
    max = 0;

    for (i=0; i<num_ccs; i++) {
        if (sol[i] > max) {
            max = (int64) round(sol[i]);
        }
    }
    return max;
}

/* Finds the lexmin solution for the constraints */
int64 *pluto_prog_constraints_lexmin_gurobi(const PlutoConstraints *cst, 
        PlutoMatrix *obj, double **val, int**index, int npar, int num_ccs)
{
    int i, is_unsat, num_sols;
    int num_vars, col_iter;
    int64 max_scale_factor;
    double *grb_obj, *fpsol, *scale_sols;
    char* vtype;

    int64 *sol;

    num_vars = cst->ncols-1;

    GRBenv *env = NULL;
    GRBmodel *lp = NULL;

    IF_DEBUG(printf("[pluto] pluto_prog_constraints_lexmin_gurobi (%d variables, %d constraints)\n",
                cst->ncols-1, cst->nrows););
    GRBloadenv(&env, NULL);
    GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);

    grb_obj = get_gurobi_objective_from_pluto_matrix(obj);

    /* Default type is GRB_CONTINIOUS - that is reals */
    vtype = NULL;

    /* Set integer constraints on all variables of the ILP */
    if (!options->lp) {
        vtype = (char*)malloc(sizeof(char)*num_vars);
        for (i = 0; i<num_vars; i++) {
            vtype[i] = GRB_INTEGER;
        }
    }

    /* Create gurobi model. Add objective and variable types during creation of the object itself */
    GRBnewmodel(env, &lp, NULL, num_vars, grb_obj, NULL, NULL, vtype, NULL);

    /* Set objective direction to minimization */
    GRBsetdblattr(lp, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);

    set_gurobi_constraints_from_pluto_constraints(lp, cst);


    if (options->debug) {
        GRBwrite(lp, "pluto.lp");
    }

    is_unsat = pluto_constraints_solve_gurobi(lp, 1e-7);

    free(grb_obj);

    if (is_unsat) {
        GRBfreemodel(lp);
        GRBfreeenv(env);
        return NULL;
    }

    if (options->lp) {
        fpsol = get_lp_solution_from_gurobi_problem(lp);
        num_sols = num_vars;
        GRBfreemodel(lp);
        /* GRBfreeenv(env); */

        /* GRBloadenv(&env, NULL); */
        lp = get_scaling_lp_gurobi(fpsol, num_sols, val, index, npar, num_ccs, env);

        if (options->debug) {
            GRBwrite(lp, "pluto-scaling-mip.lp");
        }

        scale_sols = pluto_mip_scale_solutions_gurobi(lp);

        max_scale_factor = get_max_scale_factor(scale_sols, num_ccs);

        sol =(int64*)malloc(sizeof(int64)*num_sols);
        for (i=0; i<npar+1; i++) {
            sol[i] = max_scale_factor;
        }

        col_iter = num_ccs;

        for (i=npar+1; i<cst->ncols-1; i++) {
            double x = scale_sols[col_iter++];
            IF_DEBUG(printf("c%d = %lld, ", i, (int64) round(x)););
            sol[i]=(int64)round(x);
        }

        IF_DEBUG(printf("\n"););
        free(fpsol);
        free(scale_sols);
    } else {

        sol = get_ilp_solution_from_gurobi_problem(lp);
    }

    GRBfreemodel(lp);
    GRBfreeenv(env);

    return sol;

}

/* The routine is called during the construction of FCG in pluto-lp-dfp. 
 * It returns a rational solution if the optimal solution exists else it returns NULL. */
double *pluto_fcg_constraints_lexmin_gurobi(const PlutoConstraints* cst, PlutoMatrix* obj) 
{
    double *grb_obj, *sol;
    int i, num_vars, is_unsat;
    char *vtype;

    GRBenv *env = NULL;
    GRBmodel *lp = NULL;
    GRBloadenv(&env, NULL);
    GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);

    num_vars = cst->ncols-1;

    grb_obj = get_gurobi_objective_from_pluto_matrix(obj);

    vtype = NULL;

    /* Set integer constraints on all variables of the ILP */
    if (!options->lp) {
        vtype = (char*)malloc(sizeof(char)*num_vars);
        for (i = 0; i<num_vars; i++) {
            vtype[i] = GRB_INTEGER;
        }
    }
    
    /* Create gurobi model. Add objective and variable types during creation of the object itself */
    GRBnewmodel(env, &lp, NULL, num_vars, grb_obj, NULL, NULL, vtype, NULL);

    set_gurobi_constraints_from_pluto_constraints(lp, cst);

    if (options->debug) {
        GRBwrite(lp, "pluto-pairwise-constraints-gurobi.lp");
    }

    is_unsat = pluto_constraints_solve_gurobi(lp, 1e-7);

    free(grb_obj);

    if (is_unsat) {
        return NULL;
    } else {
        sol = get_lp_solution_from_gurobi_problem(lp);
        return sol;
    }
}
#endif
