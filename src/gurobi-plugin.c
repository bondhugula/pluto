#include "constraints.h"
#include "pluto.h"
#include <math.h>
#include <assert.h>
#include <gurobi_c.h>

inline void check_error_gurobi(GRBmodel *lp, int error)
{
    if(error) {
        GRBenv *env;
        env = GRBgetenv(lp);
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);

    }
}

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

    /* non_zero = pluto_constraints_get_num_non_zero_coeffs(cst); */
    /* row = (int*)malloc((non_zero)*sizeof(int)); */
    /* col = (int*)malloc((non_zero)*sizeof(int)); */
    /* coeff = (double*)malloc((non_zero)*sizeof(double)); */

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

/* Retrives ilp solution from the input glpk problem. 
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

int64 *pluto_prog_constraints_lexmin_gurobi(const PlutoConstraints *cst, 
        PlutoMatrix *obj, double **val, int**index, int npar, int num_ccs)
{
    int i, j, is_unsat, num_sols;
    int num_vars;
    double *grb_obj;
    char* vtype;

    int64 *sol;

    double objval;
    int optim_status;

    
    num_vars = cst->ncols-1;
    
    GRBenv *env = NULL;
    GRBmodel *lp = NULL;

     IF_DEBUG(printf("[pluto] pluto_prog_constraints_lexmin_gurobi (%d variables, %d constraints)\n",
                                 cst->ncols-1, cst->nrows););
     GRBloadenv(&env, NULL);
     GRBsetintparam(env, GRB_INT_PAR_OUTPUTFLAG, 0);

     grb_obj = get_gurobi_objective_from_pluto_matrix(obj);

     vtype = (char*)malloc(sizeof(char)*num_vars);

     /* Set integer constraints on all variables of the ILP */
     if (!options->lp) {
         for (i = 0; i<num_vars; i++) {
             vtype[i] = GRB_INTEGER;
         }
     }

     /* Create gurobi model. Add objective during creation itself */
     GRBnewmodel(env, &lp, NULL, num_vars, grb_obj, NULL, NULL, vtype, NULL);

     /* Set objective direction to minimization */
     GRBsetdblattr(lp, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);

     set_gurobi_constraints_from_pluto_constraints(lp, cst);


     if (options->debug) {
         GRBwrite(lp, "pluto.lp");
     }

     GRBoptimize(lp);

     free(grb_obj);

     GRBgetintattr(lp, GRB_INT_ATTR_STATUS, &optim_status);

     if (optim_status == GRB_INFEASIBLE) {
         GRBfreemodel(lp);
         GRBfreeenv(env);
         return NULL;
     }

     if (optim_status == GRB_UNBOUNDED) {
         GRBfreemodel(lp);
         GRBfreemodel(env);
         printf("[Pluto]: Unbounded model. This appears to be a problem in Pluto's ILP/LP construction. Aborting\n");
         exit(1);
     }

     GRBgetdblattr(lp, GRB_DBL_ATTR_OBJVAL, &objval);
     IF_DEBUG(printf("Objective value: %0.6f\n",objval););
     
     sol = get_ilp_solution_from_gurobi_problem(lp);

         GRBfreemodel(lp);
         GRBfreeenv(env);

     return sol;


}
