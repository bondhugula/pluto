/* This module implements "Variable Liberalization" proposed by Mehta et al., TACO 2016. This
 * implementation has been adapted from Sanyam's Implementation with some minor changes */

#include <strings.h>
#include "constraints.h"
#include "pluto.h"
#include "ddg.h"

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

struct dist{
    PlutoConstraints* value;
    int dep;
    struct dist* next;
};

bool pluto_domain_equality(Stmt* stmt1, Stmt* stmt2)
{

    int i, j, depth, scalar;
    depth = 0;

    PlutoMatrix* mat1 = stmt1->trans,* mat2 = stmt2->trans;

    if (stmt1->dim_orig != stmt2->dim_orig || mat1->nrows!=mat2->nrows || mat1->ncols!=mat2->ncols)
        return false;


    for (i=0; i<mat1->nrows; i++) {
        scalar = 1;
        if (depth == stmt1->dim_orig) 
            break;
        for (j=0; j<mat1->ncols; j++) {
            if (j<mat1->ncols-1 && mat1->val[i][j])
                scalar = 0;
            if (mat1->val[i][j] != mat2->val[i][j])
                return false;
        }
        if (scalar == 0)
            depth++;
    }
    return true;
}


/* Returns true if the two domains are the same */
bool pluto_domain_equality1(PlutoConstraints* mat1, PlutoConstraints* mat2)
{

    int i,j;
    if (mat1->nrows!=mat2->nrows || mat1->ncols!=mat2->ncols)
        return false;
    for (i=0; i<mat1->nrows; i++) {
        for (j=0; j<mat2->ncols; j++) {
            if(mat1->val[i][j]!=mat2->val[i][j])
                return false;
        }
    }
    return true;
}


int is_on_loop(PlutoProg* prog, PlutoConstraints* dpolytope, int j)
{
    int i, count;
    count = 0;
    for (i=0; i<prog->nvar; i++) {
        if (dpolytope->val[j][i]) {
            count++;
            break;
        }
    }
    for (i=prog->nvar; i<2*prog->nvar; i++) {
        if (dpolytope->val[j][i]) {
            count++;
            break;
        }
    }

    if (count == 2)
        return true;
    else
        return false;
}

bool is_loop_independent(Dep *dep, int k, PlutoProg *prog)
{
    int j;
    for(j=0; j<dep->dpolytope->nrows; j++) {
        if ((dep->dpolytope->val[j][k]) && (is_on_loop(prog, dep->dpolytope,j)) &&
                !(dep->dpolytope->val[j][dep->dpolytope->ncols-1])) {
            return false;
        }
    }
    return true;
}

void get_unique_deps_between_statements (struct dist ***dist_, PlutoProg *prog)
{
    int ndeps, npar, i, j, k, count1, count2;
    Dep **deps;
    struct dist *d;

    ndeps = prog->ndeps;
    npar = prog->npar;

    deps = prog->deps;
    /* first collect all unique dependences between a pair of statements.
     * If the dependence is not unique, then it can be removed */

    for(i=0; i<ndeps; i++) {
        if(deps[i]->src == deps[i]->dest && IS_WAW(deps[i]->type))
            continue;
        d = dist_[deps[i]->src][deps[i]->dest];

        if(d==NULL) {

            dist_[deps[i]->src][deps[i]->dest] = (struct dist*) malloc(sizeof(struct dist));
            d = dist_[deps[i]->src][deps[i]->dest];
            d->dep = i;
            d->next = NULL;
            d->value = NULL;

            count2 = 0;
            for(j=0; j<deps[i]->dpolytope->nrows; j++) {
                count1 = 0;
                for(k=0; k<deps[i]->dpolytope->ncols-npar-1; k++) {
                    /* if the jth row has a non zero component in kth dimension */
                    /* TODO: Is there a loop enclosing the statement? */
                    if(deps[i]->dpolytope->val[j][k])
                        count1++;
                }
                /* If there is a loop enclosing the statement or if it is an equality constraint */
                if((count1>1) || prog->deps[i]->dpolytope->is_eq[j])
                    count2++;
            }

            d->value = pluto_constraints_alloc(count2, deps[i]->dpolytope->ncols);
            d->value->nrows = count2;

            count2 = 0;
            for(j=0; j<deps[i]->dpolytope->nrows; j++) {
                count1 = 0;
                for(k=0; k<deps[i]->dpolytope->ncols-npar-1; k++) {
                    if(deps[i]->dpolytope->val[j][k])
                        count1++;
                }
                if((count1>1) || deps[i]->dpolytope->is_eq[j]) {
                    for(k=0; k<deps[i]->dpolytope->ncols; k++) {
                        d->value->val[count2][k]=prog->deps[i]->dpolytope->val[j][k];
                    }
                    count2++;
                }
            }

            IF_DEBUG(printf("dep: %d (%d, %d)", i, deps[i]->src, deps[i]->dest););
        } else {

            while(d->next != NULL)
                d = d->next;

            d->next = (struct dist*) malloc(sizeof(struct dist));
            d->next->dep = i;
            d->next->value = NULL;
            d->next->next = NULL;

            count2 = 0;
            for (j=0; j<deps[i]->dpolytope->nrows; j++) {
                count1 = 0;
                for (k=0; k<deps[i]->dpolytope->ncols-npar-1; k++) {
                    if (deps[i]->dpolytope->val[j][k])
                        count1++;
                }
                if ((count1>1) || deps[i]->dpolytope->is_eq[j])
                    count2++;
            }

            d->next->value = pluto_constraints_alloc(count2,prog->deps[i]->dpolytope->ncols);
            d->next->value->nrows = count2;

            count2 = 0;
            for (j=0; j<deps[i]->dpolytope->nrows; j++) {
                count1 = 0;
                for (k=0; k<deps[i]->dpolytope->ncols-npar-1; k++) {
                    if(deps[i]->dpolytope->val[j][k])
                        count1++;
                }
                if ((count1>1) || deps[i]->dpolytope->is_eq[j]) {
                    for (k=0;k<deps[i]->dpolytope->ncols;k++) {
                        d->next->value->val[count2][k]=deps[i]->dpolytope->val[j][k];
                    }
                    count2++;
                }
            }

            struct dist* d1 = dist_[deps[i]->src][deps[i]->dest];

            while(d1 != d->next) {
                if(pluto_domain_equality1(d1->value,d->next->value))
                    break;
                d1 = d1->next;
            }


            if(d1 != d->next) {
                free(d->next);
                d->next = NULL;
            }
        }
    }

}

/*
 * Marks the dependencies that can be skipped in case of false dependencies
 * involving loop temporary variables. This implemation has been adopted 
 * from Sanyam Mehta's implementation of Variable liberalization with some 
 * minor changes 
 */
void pluto_variable_liberalize(PlutoProg *prog)
{
    int nvar, nstmts, i, j, k, ndeps;
    Graph *ddg;
    Stmt **stmts;
    Dep **deps;
    int* iter_priv_sccs_depth;

      /* dist_[i][j] points to the list of all unique dependences between statements i and j */
    struct dist ***dist_;
    struct dist *d;

    nvar = prog->nvar;
    nstmts = prog->nstmts;
    ndeps = prog->ndeps;

    stmts = prog->stmts;
    deps = prog->deps;

    ddg = prog->ddg;

    ddg_compute_scc(prog);


    dist_ = (struct dist***) malloc(sizeof(struct dist**)*nstmts);

    for (i=0; i<nstmts; i++) {
        dist_[i] = (struct dist**)malloc(sizeof(struct dist*)*nstmts);
        for (j=0; j<nstmts; j++) {
            dist_[i][j] = NULL;
        }
    }

    for (i=0; i<prog->nstmts; i++) {
        prog->stmts[i]->orig_scc_id = (int*)malloc(prog->num_hyperplanes*sizeof(int));
    }

    for (j=0; j<prog->num_hyperplanes; j++) {
        dep_satisfaction_update(prog,j);
        ddg_update(ddg,prog);
        ddg_compute_scc(prog);
        for(i=0;i<prog->nstmts;i++)
            prog->stmts[i]->orig_scc_id[j] = prog->stmts[i]->scc_id;
    }
    pluto_compute_dep_directions(prog);
    /* pluto_print_dep_directions(prog); */

    for (j=0; j<prog->ndeps; j++) {
        deps[j]->satisfied = false;
        deps[j]->skipdep = false;
    }
    ddg_update(prog->ddg,prog);
    ddg_compute_scc(prog);
    get_unique_deps_between_statements(dist_, prog);

 

    /* Variable Liberalization begin */
    iter_priv_sccs_depth = (int*) malloc(sizeof(int)*ddg->num_sccs);

    for (i=0; i<ddg->num_sccs; i++) {
        iter_priv_sccs_depth[i] = -1;
    }

    /* Get innermost loop with iteration private live ranges.
     * A live range is iteration private if the corresonding
     * RAW depenendence is not carried by that loop */
    for (i=0; i<prog->ndeps; i++) {
        if ((stmts[deps[i]->src]->scc_id == stmts[deps[i]->dest]->scc_id) && IS_RAW(deps[i]->type)) {
            IF_DEBUG(printf("Analyzing dep %d\n",i););
            /* check if the dep i is loop independent in level k. The current check seems inappropriate.
             * Ideally one should see if d.\vec{k} = 0  */
            /* This following section finds the outermost iteration private live range
             * as opposed to the innermost. This was directly taken from Sanyam's code */

            for (j=0,k=0; (j<deps[i]->dpolytope->nrows && k<nvar); j++) {
                if ((deps[i]->dpolytope->val[j][k]) && (is_on_loop(prog, deps[i]->dpolytope,j)) &&
                        !(deps[i]->dpolytope->val[j][deps[i]->dpolytope->ncols-1])) {
                    IF_DEBUG(printf("Dep %d is not loop independent in level %d\n",i,k););
                    k++;
                    j =-1;
                    continue;
                }
            }

            IF_DEBUG(printf("Dep %d, k: %d, ipscc: %d\n", i,k, iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id]););
            if ((k < iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id]) ||
                    (iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id] == -1)){
                iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id] = k;
            IF_DEBUG(printf(" Dep %d: depth of iteration private live range %d\n", i, iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id]););
            }
        }
    }

    bool *unique_deps;
    unique_deps = (bool*)malloc(ndeps * sizeof (bool));
    bzero(unique_deps, ndeps*sizeof(bool));



    for(i=0; i<nstmts; i++)
    {
        for(j=0; j<nstmts; j++)
        {
            d = dist_[i][j];
            while(d != NULL)
            {
                unique_deps[d->dep] = true;
                d=d->next;
            }
        }
    }


    int* order = (int*) calloc(nvar,sizeof(int));
    for(i=0; i<ndeps; i++)
    {
        IF_DEBUG(printf("Dep being considered for skipping : %d\n",i););
        deps[i]->temp_across = false;
        deps[i]->fuse_depth = 0;


        if (!unique_deps[i]) {
            deps[i]->skipdep = true;
        }
        /* The case where the dependence is in the same SCC */
        else if (stmts[prog->deps[i]->src]->scc_id == stmts[prog->deps[i]->dest]->scc_id) {

            for (k=0; k<nvar; k++) {
                if(is_loop_independent(deps[i], k, prog)) {
                    break;
                }
            }
            /* finds the outermost dimension in whch the dependence is loop independent */
            for (j=0,k=0; (j<deps[i]->dpolytope->nrows && k<nvar); j++) {
                if ((deps[i]->dpolytope->val[j][k]) && (is_on_loop(prog,deps[i]->dpolytope,j)) &&
                        !(deps[i]->dpolytope->val[j][deps[i]->dpolytope->ncols-1])) {
                    IF_DEBUG(printf("Dep %d is not loop independent in level %d\n",i,k););
                    k++;
                    j = -1;
                    continue;
                }
                IF_DEBUG(printf("Dep %d is loop independent in level %d\n",i,k););
            }
            IF_DEBUG(printf("depth of iteration private live range %d\n", iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id]););
            /* k will point to the outermost dimension that is loop independent */
            if(k<iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id]) {
                /* these dependences should be removed if possible as these may restrict transformations */
                int l = 0;
                for(j=0; (j<prog->num_hyperplanes && (l<min(stmts[deps[i]->src]->dim_orig, stmts[deps[i]->dest]->dim_orig))); j++) {
                    if(get_loop_type(stmts[deps[i]->src],j) == H_LOOP) { /* its a loop-hyperplane */
                        l++;
                        if((stmts[deps[i]->src]->orig_scc_id[j] != stmts[deps[i]->dest]->orig_scc_id[j]) && (k<l)) {
                            deps[i]->skipdep = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    for(i=0; i<ndeps; i++) {
        deps[i]->temp_across = false;
        deps[i]->fuse_depth = 0;

        if (stmts[deps[i]->src]->scc_id != stmts[deps[i]->dest]->scc_id &&
                (iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id]>0) &&
                (iter_priv_sccs_depth[stmts[deps[i]->dest]->scc_id]>0)) {
            PlutoConstraints* polytope;
            polytope = deps[i]->dpolytope;

            for(j=0; j<polytope->nrows; j++) {
                if(polytope->is_eq[j] && !is_on_loop(prog,polytope,j)) {
                    IF_DEBUG(printf("Removing row for dependence %d\n",i););
                    pluto_constraints_remove_row(polytope,j); // absence of a row is fine; an equality may be inserted later to relax scheduling
                    j=-1;
                    continue;
                }
            }
            /* The following is to particularly handle cases when the loop order
             * in the candidate nests is different */
            for (j=0; j<nvar; j++) {
                order[j] = 0;
            }

            for(j=0; j<polytope->nrows; j++) {
                if(is_on_loop(prog,polytope,j)) {
                    for(k=0; k<nvar; k++) {
                        if(polytope->val[j][k]) order[k] = -1;
                    }
                    for(k=0; k<nvar; k++) {
                        if(polytope->val[j][k+nvar]) order[k] = -1;
                    }
                }
            }

            for(k=0; k<min(iter_priv_sccs_depth[stmts[deps[i]->src]->scc_id],
                        iter_priv_sccs_depth[stmts[deps[i]->dest]->scc_id]) && order[k]!=-1; k++) {
                IF_DEBUG(printf("Adding equality constraints for dep %d\n",i););
                pluto_constraints_add_equality(polytope);
                polytope->val[j][k] = 1;
                polytope->val[j][k+nvar] = -1;
                j++;
            }
            deps[i]->temp_across = true;
            deps[i]->fuse_depth = k;
        } // if
    }

    /* Cleanup dist */
    for (i=0; i<nstmts; i++){
        for (j=0; j<nstmts; j++) {
            struct dist *tmp;
            d = dist_[i][j];
            while (d != NULL) {
                tmp = d;
                d = d->next;
                if (tmp->value != NULL) {
                    pluto_constraints_free(tmp->value);
                }
                free(tmp);
            }
        }
        free(dist_[i]);
    }

    for (j=0; j<nstmts; j++) {
        free(stmts[j]->orig_scc_id);
    }
    free(dist_);

    if (options->debug) {
        printf("[Pluto] Unique dependences \n");
        for (i=0; i<ndeps ; i++) {
            if (unique_deps[i]) {
                printf("%d,",i);
            }
        }
        printf("\n");
        printf("[Pluto] Dependences skipped by variable liberalization \n");
        for (i=0; i<ndeps; i++) {
            if (deps[i]->skipdep)
                printf("%d,", i);
        }
        printf("\n");
    }

    free (unique_deps);


    free(iter_priv_sccs_depth);
    free(order);

    pluto_dep_satisfaction_reset(prog);
}

