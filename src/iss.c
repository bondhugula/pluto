/* Index set splitting */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "pluto.h"
#include "math_support.h"
#include "constraints.h"

/*
 * The affine form of the Farkas lemma
 * cst is the domain on which the affine form described in phi is non-negative
 *
 * Returns: constraints on the variables that define the rows of \phi (each
 * row of phi is an affine function these variables). The output is thus a
 * constraint set with phi->ncols-1 variables (farkas multipliers are eliminated
 * by FM). In effect, this allows one to linearize the constraint describing
 * \phi.
 *
 * The rows of phi's correspond to the coefficients of the variables in dom in
 * that order with the last row of phi representing the translation part of
 * the affine form. The number of rows in phi is thus the same as the number
 * of columns in dom (number of dom dimensions + 1).
 *
 * Each row of phi itself is an affine function of a set of variables
 * (phi->ncols-1 variables); the last column corresponds to the constant.
 *
 * (c_1 + c_2)*i + (c_2 - c3)*j + (c1 + c2 + c3 + 1) >= 0 over a domain (dom) on
 * (i,j)
 *
 * Here phi would be
 *
 * [1 1 0  0]
 * [0 1 -1 0]
 * [1 1 1  1]
 *
 * Let cst have faces (inequalities representing non-negative half-spaces) f1,
 * f2, ..., fn
 *
 * The affine form of the farkas lemma now says that
 *
 * (c_1 + c_2)*i + (c_2 - c3)*j + (c1 + c2 + c3 + 1) = \lambda_1*f1 +
 * \lambda_2*f2 + ... + \lambda_n*fn + \lambda_0,
 * with all \lambda_i >= 0
 *
 * Eliminate Farkas multipliers by FM and return constraints in c_1, c_2, c_3
 *
 * */
PlutoConstraints *farkas_affine(const PlutoConstraints *dom, const PlutoMatrix *phi)
{
    int i, j;

    /* Convert everything into inequalities of >= 0 form */
    PlutoConstraints *idom = pluto_constraints_to_pure_inequalities_single(dom);

    // printf("Initial constraints\n");
    // pluto_constraints_pretty_print(stdout, idom);

    // printf("phi matrix\n");
    // pluto_matrix_print(stdout, phi);
    
    /* Add a trivial row (1 >= 0) for the translation Farkas
     * multiplier (\lambda_0) so that the non-negative linear
     * combination of the faces is modeled naturally below */
    pluto_constraints_add_inequality(idom);
    idom->val[idom->nrows-1][idom->ncols-1] = 1;

    assert(phi->nrows == idom->ncols);

    /*
     * Farkas space
     * idom->ncols equalities (one for each of the idom->ncols-1 variables)
     * and one for the constant part, followed by idom->nrows constraints for
     * the farkas multipliers *
     * idom->nrows is the number of Farkas multipliers
     * format: [phi->ncols-1 vars, idom->nrows farkas multipliers, const]
     * the translation farkas multiplier appears last
     *
     *    Eg: [c_1, c_2, c_3, l_1, l_2, ..., l_n, l_0, 1]
     */
    PlutoConstraints *farkas = pluto_constraints_alloc(
            idom->ncols+idom->nrows, phi->ncols+idom->nrows);
    farkas->nrows = idom->ncols+idom->nrows;

    int farkas_offset = phi->ncols-1;
    /* First idom->ncols equalities */
    for (i=0; i<idom->ncols; i++) {
        farkas->is_eq[i] = 1;
        for (j=0; j<phi->ncols-1; j++) {
            farkas->val[i][j] = phi->val[i][j];
        }
        for (j=0; j<idom->nrows; j++) {
            farkas->val[i][farkas_offset+j] = -idom->val[j][i];
        }
        farkas->val[i][farkas_offset + idom->nrows+1] = phi->val[i][phi->ncols-1];
    }

    /* All farkas multipliers are non-negative */
    for (j=0; j<idom->nrows; j++) {
        farkas->is_eq[idom->ncols+j] = 0;
        farkas->val[idom->ncols+j][farkas_offset + j] = 1;
    }

    // printf("After equating both sides\n");
    // pluto_constraints_pretty_print(stdout, farkas);
    for (i=0; i<idom->nrows; i++) {
        int best_elim = pluto_constraints_best_elim_candidate(farkas, idom->nrows-i);
        fourier_motzkin_eliminate(farkas, best_elim);
        // printf("After eliminating %d multiplier\n", i);
        // printf("%d rows\n", farkas->nrows);
        // pluto_constraints_pretty_print(stdout, farkas);
    }
    assert(farkas->ncols == phi->ncols);
    
    // printf("After farkas multiplier elimination\n");
    // pluto_constraints_pretty_print(stdout, farkas);

    pluto_constraints_free(idom);

    return farkas;
}


PlutoConstraints **get_lin_ind_constraints(PlutoMatrix *mat, int *orthonum)
{
    int i, j, k;
    PlutoConstraints **orthcst;
    isl_ctx *ctx;
    isl_mat *h;

    assert(mat != NULL);

    int ndim = mat->ncols;

    ctx = isl_ctx_alloc();
    assert(ctx);

    // printf("Input matrix\n");
    // pluto_matrix_print(stdout, mat);

    h = isl_mat_alloc(ctx, mat->nrows, ndim);

    for (i=0; i<mat->ncols; i++) {
        for (j=0; j<mat->nrows; j++) {
            h = isl_mat_set_element_si(h, j, i, mat->val[j][i]);
        }
    }

    h = isl_mat_right_kernel(h);

    PlutoMatrix *ortho = pluto_matrix_from_isl_mat(h);

    isl_mat_free(h);

    orthcst = (PlutoConstraints **) malloc((ndim+1)*sizeof(PlutoConstraints *)); 

    for (i=0; i<ndim+1; i++)  {
        orthcst[i] = pluto_constraints_alloc(1, ndim+1);
        orthcst[i]->ncols = ndim+1;
    }

    /* All non-negative orthant only */
    /* An optimized version where the constraints are added as
     * c_1 >= 0, c_2 >= 0, ..., c_n >= 0, c_1+c_2+..+c_n >= 1
     *
     * basically only look in the orthogonal space where everything is
     * non-negative
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
    // printf("Ortho matrix\n");
    // pluto_matrix_print(stdout, ortho); 

    for (i=0; i<ortho->ncols; i++) {
        for (j=0; j<ndim; j++) {
            orthcst[i]->val[0][j] = ortho->val[j][i];
        }
        orthcst[i]->nrows = 1;
        orthcst[i]->val[0][ndim] = -1;
        orthcst[i]->val[0][ndim] = 0;
    }

    // pluto_matrix_print(stdout, stmt->trans);

    if (ortho->ncols >= 1)  {
        /* Sum of all of the above is the last constraint */
        for(j=0; j<ndim+1; j++)  {
            for (i=0; i<ortho->ncols; i++) {
                orthcst[ortho->ncols]->val[0][j] += orthcst[i]->val[0][j];
            }
        }
        orthcst[ortho->ncols]->nrows = 1;
        orthcst[ortho->ncols]->val[0][ndim] = -1;
        *orthonum = ortho->ncols+1;
    }else *orthonum = 0;


    // printf("Ortho constraints: %d set(s)\n", *orthonum);
    // for (i=0; i<*orthonum; i++) {
    	// print_polylib_visual_sets("li", orthcst[i]);
        // pluto_constraints_print(stdout, orthcst[i]);
    // }

    /* Free the unnecessary ones */
    for (i=*orthonum; i<ndim+1; i++)    {
        pluto_constraints_free(orthcst[i]);
    }

    pluto_matrix_free(ortho);
    isl_ctx_free(ctx);

    return orthcst;
}



/*
 * Index Set Splitting for close-to mid-point cutting
 *
 * Refer to PACT'14 paper on tiling periodic domains for
 * the formulation
 */
PlutoConstraints *pluto_iss(const PlutoConstraints **doms, int ndoms,
        int npar, PlutoConstraints *indcst)
{
    int i, j, k, ndim;

    if (ndoms == 0) return NULL;

    assert(ndoms >= 0);

    const PlutoConstraints *dom0 = doms[0];

    for (i=0; i<ndoms; i++) {
        assert((doms[i]->ncols - 1 - npar)%2 == 0);
    }
    ndim = (doms[0]->ncols - 1 - npar)/2;

    /* 
     * [m | sigma(h) | h | P r| const ] 
     *
     * The ISS is given by h.i = P.p + r, i \in dom
     * m = bound on distance from mid-point
     *
     * */
    int iss_cst_width = 2+ndim+npar+1+1;
    PlutoConstraints *cst = pluto_constraints_alloc(10, iss_cst_width);

    for (k=0; k<ndoms; k++) {
        // printf("[iss] Domain %d\n", k);
        const PlutoConstraints *dom = doms[k];

        /* Linearize m + 2v(p) - h.s - h.t >= 0 */
        PlutoMatrix *mat = pluto_matrix_alloc(dom->ncols, iss_cst_width);
        pluto_matrix_initialize(mat, 0);

        for (i=0; i<ndim; i++) {
            mat->val[i][2+i] = -1;
        }
        for (i=0; i<ndim; i++) {
            mat->val[ndim+i][2+i] = -1;
        }
        for (i=0; i<npar; i++) {
            mat->val[2*ndim+i][2+ndim+i] = 2;
        }
        /* m + 2*r for the const part */
        mat->val[2*ndim+npar][2+ndim+npar] = 2;
        mat->val[2*ndim+npar][0] = 1;

        /*
         * cst is of the following form
         * [m | sigma(h) | h | P r| const ] 
         */
        PlutoConstraints *cst1 = farkas_affine(dom, mat);

        /* Linearize: m - 2v(p) + h.s + h.t >= 0 */
        pluto_matrix_initialize(mat, 0);

        for (i=0; i<ndim; i++) {
            mat->val[i][2+i] = 1;
        }
        for (i=0; i<ndim; i++) {
            mat->val[ndim+i][2+i] = 1;
        }
        for (i=0; i<npar; i++) {
            mat->val[2*ndim+i][2+ndim+i] = -2;
        }
        /* m + 2*r for the const part */
        mat->val[2*ndim+npar][2+ndim+npar] = -2;
        mat->val[2*ndim+npar][0] = 1;

        PlutoConstraints *cst2 = farkas_affine(dom, mat);
        pluto_matrix_free(mat);

        pluto_constraints_add(cst, cst1);
        pluto_constraints_add(cst, cst2);
        pluto_constraints_free(cst1);
        pluto_constraints_free(cst2);
    }

    /* Set sigma h */
    pluto_constraints_add_equality(cst);
    /* \sigma h_i >= 1 */
    for (j=2; j<2+ndim; j++) {
        cst->val[cst->nrows-1][j] = -1;
    }
    cst->val[cst->nrows-1][1] = 1;

    /* Avoid trivial zero solution of h */
    PlutoConstraints *nz = pluto_constraints_alloc(1, iss_cst_width);
    nz->nrows = 1;
    for (j=0; j<ndim; j++) {
        nz->val[0][2+j] = 1;
    }
    nz->val[0][nz->ncols-1] = -1;
    pluto_constraints_add(cst, nz);

    /* Add linear independence constraints */
    if (indcst) {
        /* for m and sigma h */
        pluto_constraints_add_dim(indcst, 0, NULL);
        pluto_constraints_add_dim(indcst, 0, NULL);

        for (j=0; j<npar+1; j++) {
            pluto_constraints_add_dim(indcst, ndim+2, NULL);
        }
        pluto_constraints_add(cst, indcst);
    }

    // pluto_constraints_pretty_print(stdout, cst);
    int64 *sol = pluto_constraints_solve(cst, 0);

    pluto_constraints_free(cst);
    pluto_constraints_free(nz);

    if (sol) {
        PlutoConstraints *h = pluto_constraints_alloc(1, ndim+npar+1);
        h->nrows = 1;

        h->is_eq[0] = 1;
        for (j=0; j<ndim; j++) {
            h->val[0][j] = sol[2+j];
        }
        for (j=0; j<npar; j++) {
            h->val[0][ndim+j] = -sol[2+ndim+j];
        }
        h->val[0][ndim+npar] = sol[2+ndim+npar];
        pluto_constraints_set_names_range(h, dom0->names, 0, 0, ndim);
        pluto_constraints_set_names_range(h, dom0->names, ndim, 2*ndim, npar);

        printf("[iss] m = %lld\n", sol[0]);
        printf("[iss] h (cut) is ");
        pluto_constraints_compact_print(stdout, h);
        free(sol);
        return h;
    }else{
        printf("[iss] No solution (close to mid-point cut)\n");
        return NULL;
    }
}


#if 0
void pluto_iss(const PlutoConstraints **doms, int ndoms, int npar)
{
    int j;
    int orthonum;
    PlutoConstraints *h = NULL;
    PlutoConstraints *indcst_i = NULL;

    const PlutoConstraints *dom0 = doms[0];

    PlutoMatrix *H = pluto_matrix_alloc(3, (dom0->ncols-1-npar)/2);
    H->nrows = 0;

    //do{
        h = iss(doms, ndoms, npar, indcst_i);
        if (h != NULL) {
            pluto_matrix_add_row(H, H->nrows);

            for (j=0; j<h->ncols-1-npar; j++) {
                H->val[H->nrows-1][j] = h->val[0][j];
            }

            PlutoConstraints **indcst = get_lin_ind_constraints(H, &orthonum);
            /* There will be at least two or zero */
            assert(orthonum != 1);
            if (orthonum >= 2) {
                // pluto_constraints_set_names_range(indcst[orthonum-1], dom->names, 0, 0, h->ncols-1-npar);
                // printf("Linear independence constraints\n");
                // pluto_constraints_pretty_print(stdout, indcst[orthonum-1]);
                indcst_i = pluto_constraints_dup(indcst[orthonum-1]);
            }else indcst_i = NULL;
        }
    //}while (h != NULL && indcst_i != NULL);
}
#endif

int is_long_bidirectional_dep(const Dep *dep, int dim, int npar)
{
    assert(dep->src == dep->dest);

    int ndim = (dep->dpolytope->ncols-1-npar)/2;

    assert(dim >= 0);
    assert(dim <= ndim-1);

    PlutoConstraints *dpolyc = pluto_constraints_dup(dep->dpolytope);
    pluto_constraints_add_dim(dpolyc, 0, NULL);
    pluto_constraints_add_equality(dpolyc);
    dpolyc->val[dpolyc->nrows-1][0] = 1;
    dpolyc->val[dpolyc->nrows-1][1+dim] = 1;
    dpolyc->val[dpolyc->nrows-1][1+ndim+dim] = -1;

    int64 lb, ub;
    int retval1, retval2;

    retval1 = pluto_constraints_get_const_lb(dpolyc, 0, &lb);
    retval2 = pluto_constraints_get_const_ub(dpolyc, 0, &ub);

    // if (retval1)  printf("lb = %lld\n", lb);
    // if (retval2)  printf("ub = %lld\n", ub);

    return !(retval1 && retval2 && ub - lb <= 5);
}


void pluto_iss_dep(PlutoProg *prog)
{
    int ndeps = prog->ndeps;
    int npar = prog->npar;

    if (prog->nstmts == 0 || prog->nstmts >= 2) return;

    if (ndeps == 0) return;

    int ndim = prog->stmts[0]->dim;

    int is_long[prog->ndeps][ndim];
    int num_long_deps[ndim];

    bzero(num_long_deps, ndim*sizeof(int));

    int i, j;

    for (i=0; i<ndeps; i++) {
        if (prog->deps[i]->src != prog->deps[i]->dest) continue;
        int ndim = (prog->deps[i]->dpolytope->ncols-1-npar)/2;
        for (j=0; j<ndim; j++) {
            is_long[i][j] = is_long_bidirectional_dep(prog->deps[i], j, npar);
            num_long_deps[j] += is_long[i][j];
        }
    }

    for (j=0; j<ndim; j++) {
        printf("[iss] Dimension %d has %d long deps\n", j, num_long_deps[j]);
    }

    PlutoConstraints ***long_dep_doms = (PlutoConstraints ***) malloc(sizeof(PlutoConstraints **)*ndim);
    for (i=0; i<ndim; i++) {
        if (num_long_deps[i] >= 1) {
        long_dep_doms[i] = malloc(num_long_deps[i]*sizeof(PlutoConstraints *));
        }else long_dep_doms[i] = NULL;
    }

    for (j=0; j<ndim; j++) {
        // printf("Along dimension %d\n", j);
        int q = 0;
        for (i=0; i<ndeps; i++) {
            if (is_long[i][j]) {
                assert(q <= num_long_deps[j]-1);
                long_dep_doms[j][q++] = prog->deps[i]->dpolytope;
                // pluto_constraints_compact_print(stdout, prog->deps[i]->dpolytope);
            }
        }
    }

    for (i=0; i<ndim; i++) {
        if (num_long_deps[i] >= 1) {
            pluto_iss((const PlutoConstraints **) long_dep_doms[i], num_long_deps[i], npar, NULL);
        }
    }
}

