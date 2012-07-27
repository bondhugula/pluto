#include "pluto.h"
#include "constraints.h"
#include "candl/candl.h"
#include "program.h"
#include "pluto/libpluto.h"
#include "isl/map.h"

PlutoOptions *options;

/*
 * Output schedules are isl relations that have dims in the order
 * isl_dim_out, isl_dim_in, div, param, const
 */
__isl_give isl_union_map *pluto_schedule(isl_union_set *domains, 
        isl_union_map *dependences, 
        PlutoOptions *options_l)
{
    int i;
    isl_ctx *ctx;

    ctx = isl_ctx_alloc();

    options = options_l;

    PlutoProg *prog = pluto_prog_alloc();

    prog->nvar = -1;

    prog->nstmts = isl_union_set_n_set(domains);

    if (prog->nstmts >= 1) {
        prog->stmts = (Stmt **)malloc(prog->nstmts * sizeof(Stmt *));
    }else{
        prog->stmts = NULL;
    }

    for (i=0; i<prog->nstmts; i++) {
        prog->stmts[i] = NULL;
    }

    extract_stmt_domains(domains, prog->stmts);

    for (i=0; i<prog->nstmts; i++) {
        prog->nvar = PLMAX(prog->nvar, prog->stmts[i]->dim);
    }

    if (prog->nstmts >= 1) {
        Stmt *stmt = prog->stmts[0];
        prog->npar = stmt->domain->ncols - stmt->dim - 1;
    }else prog->npar = 0;

    prog->ndeps = isl_union_map_n_map(dependences);

    prog->deps = (Dep **)malloc(prog->ndeps * sizeof(Dep *));
    for (i=0; i<prog->ndeps; i++) {
        prog->deps[i] = pluto_dep_alloc();
    }
    extract_deps(prog->deps, prog->ndeps, prog->stmts,
            dependences, CANDL_RAW);

    pluto_auto_transform(prog);
    
    /* Extract schedules */
    isl_space *space = isl_space_alloc(ctx, prog->npar, 0, 0);

    isl_union_map *schedules = isl_union_map_empty(space);

    for (i=0; i<prog->nstmts; i++) {
        Stmt *stmt = prog->stmts[i];
        PlutoConstraints *sched = pluto_stmt_get_schedule(stmt);

        isl_basic_map *bmap;
        isl_map *map;

        bmap = isl_basic_map_from_pluto_constraints(ctx, sched, 
                stmt->domain->ncols-1, stmt->trans->nrows, prog->npar);
        map = isl_map_from_basic_map(bmap);
        schedules = isl_union_map_union(schedules, isl_union_map_from_map(map));
    }

    pluto_prog_free(prog);
    isl_ctx_free(ctx);

    return schedules;
}
