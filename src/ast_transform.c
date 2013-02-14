#include <stdlib.h>
#include <stdio.h>

#include "pluto.h"
#include "program.h"
#include "ast_transform.h"

#include "cloog/cloog.h"

/*
 * Clast-based parallel loop marking */
void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions)
{
    int i, j, nloops, nstmts, nploops;
    struct clast_for **loops;
    int *stmts;
    assert(root != NULL);

    // int filter[1] = {1};

    Ploop **ploops = pluto_get_dom_parallel_loops(prog, &nploops);

    // pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);

    IF_DEBUG(printf("[pluto_mark_parallel] parallel loops\n"););
    IF_DEBUG(pluto_loops_print(ploops, nploops););

    // clast_pprint(stdout, root, 0, cloogOptions);

    for (i=0; i<nploops; i++) {
        char iter[5];
        sprintf(iter, "t%d", ploops[i]->depth+1);
        int *stmtids = malloc(ploops[i]->nstmts*sizeof(int));
        for (j=0; j<ploops[i]->nstmts; j++) {
            stmtids[j] = ploops[i]->stmts[j]->id+1;
        }

        // IF_DEBUG(printf("Looking for loop\n"););
        // IF_DEBUG(pluto_loop_print(ploops[i]););
        // IF_DEBUG(printf("%s in \n", iter););
        // IF_DEBUG(clast_pprint(stdout, root, 0, cloogOptions););

        ClastFilter filter = {iter, stmtids, ploops[i]->nstmts, subset};
        clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);

        /* There should be at least one */
        if (nloops==0) {
            /* Sometimes loops may disappear (1) tile size larger than trip count
             * 2) it's a scalar dimension but can't be determined from the
             * trans matrix */
            printf("Warning: parallel poly loop not found in AST\n");
            continue;
        }else{
            for (j=0; j<nloops; j++) {
                // printf("Marking %s parallel\n", loops[j]->iterator);
                if (options->parallel) {
                    loops[j]->parallel = CLAST_PARALLEL_OMP;
                }
                loops[j]->private_vars = strdup("lbv,ubv");
            }
        }
        free(stmtids);
        free(loops);
        free(stmts);
    }

    pluto_loops_free(ploops, nploops);
}


/*
 * Clast-based vector loop marking */
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog,
        CloogOptions *cloogOptions)
{
    int i, j, nloops, nstmts, nploops;
    struct clast_for **loops;
    int *stmts;
    assert(root != NULL);

    Ploop **ploops = pluto_get_parallel_loops(prog, &nploops);

    // pluto_print_depsat_vectors(prog->deps, prog->ndeps, prog->num_hyperplanes);
    // clast_pprint(stdout, root, 0, cloogOptions);

    for (i=0; i<nploops; i++) {
        /* Only the innermost ones */
        if (!pluto_is_loop_innermost(ploops[i], prog)) continue;

        IF_DEBUG(printf("[pluto_mark_vector] marking loop\n"););
        IF_DEBUG(pluto_loop_print(ploops[i]););
        char iter[5];
        sprintf(iter, "t%d", ploops[i]->depth+1);
        int *stmtids = malloc(ploops[i]->nstmts*sizeof(int));
        for (j=0; j<ploops[i]->nstmts; j++) {
            stmtids[j] = ploops[i]->stmts[j]->id+1;
        }

        // IF_DEBUG(printf("Looking for loop\n"););
        // IF_DEBUG(pluto_loop_print(ploops[i]););
        // IF_DEBUG(printf("%s in \n", iter););
        // IF_DEBUG(clast_pprint(stdout, root, 0, cloogOptions););

        ClastFilter filter = {iter, stmtids, ploops[i]->nstmts, subset};
        clast_filter(root, filter, &loops, &nloops, &stmts, &nstmts);

        /* There should be at least one */
        if (nloops==0) {
            /* Sometimes loops may disappear (1) tile size larger than trip count
             * 2) it's a scalar dimension but can't be determined from the
             * trans matrix */
            printf("Warning: parallel poly loop not found in AST\n");
            continue;
        }
        for (j=0; j<nloops; j++) {
            // printf("Marking %s ivdep\n", loops[j]->iterator);
            loops[j]->parallel += CLAST_PARALLEL_VEC;
        }
        free(stmtids);
        free(loops);
        free(stmts);
    }

    pluto_loops_free(ploops, nploops);
}
