/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2008 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
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
 * program.c
 *
 * This file contains functions that do the job interfacing the PLUTO 
 * core to the frontend and related matters
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "pluto.h"
#include "math_support.h"
#include "program.h"

#include "scoplib/statement.h"

/* Return a copy of the statement */
Stmt *stmt_copy(Stmt *src)
{
    Stmt *dest = (Stmt *) malloc(sizeof(Stmt));

    *dest = *src;

    dest->domain = pluto_constraints_copy(dest->domain, src->domain);
    dest->trans = pluto_matrix_copy(src->trans);

    return dest;
}


PlutoConstraints *clan_matrix_to_pluto_constraints(scoplib_matrix_p clanMatrix)
{
    // candl_matrix_print(stdout, candlMatrix);
    int has_equalities = 0;

    int i, j;

    /* Does it have any equalities at all? */
    for (i=0; i<clanMatrix->NbRows; i++)   {
        if (clanMatrix->p[i][0] == 0) {
            has_equalities = 1;
            break;
        }
    }

    PlutoConstraints *pmat;
    if (has_equalities) {
        /* An extra inequality will be added to capture the inequalities */
        pmat = pluto_constraints_alloc(clanMatrix->NbRows+1, clanMatrix->NbColumns-1);
        pmat->nrows = clanMatrix->NbRows+1;
    }else{
        pmat = pluto_constraints_alloc(clanMatrix->NbRows, clanMatrix->NbColumns-1);
        pmat->nrows = clanMatrix->NbRows;
    }

    pmat->ncols = clanMatrix->NbColumns-1;

    for (i=0; i<clanMatrix->NbRows; i++)   {
        for (j=0; j<pmat->ncols; j++)   {
#ifdef PIP_WIDTH_MP
            pmat->val[i][j] = mpz_get_si(clanMatrix->p[i][j+1]);
#else
            pmat->val[i][j] = (int) clanMatrix->p[i][j+1];
#endif
        }
    }

    if (has_equalities) {
        /* Last row is sigma (equalities) <= 0 */
        for (j=0; j<pmat->ncols; j++)   {
            pmat->val[pmat->nrows-1][j] = 0;
        }

        for (i=0; i<clanMatrix->NbRows; i++)   {
#ifdef PIP_WIDTH_MP
            if (mpz_get_si(clanMatrix->p[i][0]) == 0) {
#else
            if (clanMatrix->p[i][0] == 0) {
#endif
                for (j=0; j<pmat->ncols; j++)   {
#ifdef PIP_WIDTH_MP
                    pmat->val[pmat->nrows-1][j] -= mpz_get_si(clanMatrix->p[i][j+1]);
#else
                    pmat->val[pmat->nrows-1][j] -= clanMatrix->p[i][j+1];
#endif
                }
            }
        }
    }

    return pmat;
}


PlutoConstraints *candl_matrix_to_pluto_constraints(CandlMatrix *candlMatrix)
{
    // candl_matrix_print(stdout, candlMatrix);
    int has_equalities = 0;

    int i, j;

    /* Does it have any equalities at all? */
    for (i=0; i<candlMatrix->NbRows; i++)   {
#ifdef PIP_WIDTH_MP
        if (mpz_get_si(candlMatrix->p[i][0]) == 0) {
#else
        if (candlMatrix->p[i][0] == 0) {
#endif
            has_equalities = 1;
            break;
        }
    }

    PlutoConstraints *pmat;
    if (has_equalities) {
        /* An extra inequality will be added to capture the inequalities */
        pmat = pluto_constraints_alloc(candlMatrix->NbRows+1, candlMatrix->NbColumns-1);
        pmat->nrows = candlMatrix->NbRows+1;
    }else{
        pmat = pluto_constraints_alloc(candlMatrix->NbRows, candlMatrix->NbColumns-1);
        pmat->nrows = candlMatrix->NbRows;
    }

    pmat->ncols = candlMatrix->NbColumns-1;

    for (i=0; i<candlMatrix->NbRows; i++)   {
        for (j=0; j<pmat->ncols; j++)   {
#ifdef PIP_WIDTH_MP
            pmat->val[i][j] = mpz_get_si(candlMatrix->p[i][j+1]);
#else
            pmat->val[i][j] = (int) candlMatrix->p[i][j+1];
#endif
        }
    }

    if (has_equalities) {
        /* Last row is sigma (equalities) <= 0 */
        for (j=0; j<pmat->ncols; j++)   {
            pmat->val[pmat->nrows-1][j] = 0;
        }

        for (i=0; i<candlMatrix->NbRows; i++)   {
#ifdef PIP_WIDTH_MP
            if (mpz_get_si(candlMatrix->p[i][0]) == 0) {
#else
            if (candlMatrix->p[i][0] == 0) {
#endif
                for (j=0; j<pmat->ncols; j++)   {
#ifdef PIP_WIDTH_MP
                    pmat->val[pmat->nrows-1][j] -= mpz_get_si(candlMatrix->p[i][j+1]);
#else
                    pmat->val[pmat->nrows-1][j] -= candlMatrix->p[i][j+1];
#endif
                }
            }
        }
    }
    // pluto_matrix_print(stdout, pmat);

    return pmat;
}

/* Read dependences from candl structures */
static Dep *deps_read(CandlDependence *candlDeps, PlutoProg *prog)
{
    int i, ndeps;
    Dep *deps;
    int npar = prog->npar;
    Stmt *stmts = prog->stmts;

    ndeps = candl_num_dependences(candlDeps);

    deps = (Dep *) malloc(ndeps*sizeof(Dep));

    CandlDependence *candl_dep = candlDeps;

    candl_dep = candlDeps;

    IF_DEBUG(candl_dependence_pprint(stdout, candl_dep));

    /* Dependence polyhedra information */
    for (i=0; i<ndeps; i++)  {

        Dep *dep = &deps[i];

        dep->id = i;

        // candl_matrix_print(stdout, candl_dep->domain);
        dep->dpolytope = candl_matrix_to_pluto_constraints(candl_dep->domain);

        /* Get rid of rows that are all zero */
        int r, c;
        bool *remove = (bool *) malloc(sizeof(bool)*dep->dpolytope->nrows);
        for (r=0; r<dep->dpolytope->nrows; r++) {
            for (c=0; c<dep->dpolytope->ncols; c++) {
                if (dep->dpolytope->val[r][c] != 0) {
                    break;
                }
            }
            if (c == dep->dpolytope->ncols) {
                remove[r] = true;
            }else{
                remove[r] = false;
            }
        }
        int orig_nrows = dep->dpolytope->nrows;
        int del_count = 0;
        for (r=0; r<orig_nrows; r++) {
            if (remove[r])  {
                pluto_constraints_remove_row(dep->dpolytope, r-del_count);
                del_count++;
            }
        }
        free(remove);

        dep->type = candl_dep->type;

        int src_stmt_id = candl_dep->source->label;
        int target_stmt_id = candl_dep->target->label;

        dep->src = src_stmt_id;
        dep->dest = target_stmt_id;

        int src_dim = stmts[src_stmt_id].dim;
        int target_dim = stmts[target_stmt_id].dim;

        assert(candl_dep->domain->NbColumns-1 == src_dim+target_dim+npar+1);

        /* Initialize other fields used for auto transform */
        dep->satisfied = false;
        dep->satisfaction_level = -1;

        candl_dep = candl_dep->next;
    }

    return deps;
}

void dep_print(FILE *fp, Dep *dep)
{
    fprintf(fp, "--- Dep %d from S%d to S%d, Type: ",
            dep->id+1, dep->src+1, dep->dest+1);

    switch (dep->type) {
        case CANDL_UNSET : fprintf(fp, "UNSET"); break;
        case CANDL_RAW   : fprintf(fp, "RAW")  ; break;
        case CANDL_WAR   : fprintf(fp, "WAR")  ; break;
        case CANDL_WAW   : fprintf(fp, "WAW")  ; break;
        case CANDL_RAR   : fprintf(fp, "RAR")  ; break;
        default : fprintf(fp, "unknown"); break;
    }

    fprintf(fp, "\n\n");

    fprintf(fp, "Dependence polyhedron\n");
    pluto_constraints_pretty_print(fp, dep->dpolytope);
}


void deps_print(FILE *fp, Dep *deps, int ndeps)
{
    int i;
    for (i=0; i<ndeps; i++) {
        dep_print(fp, &deps[i]);
    }
}


/* Read statement info from Clan's structures */
static Stmt *stmts_read(scoplib_scop_p scop, int npar, int nvar)
{
    int i, j;
    Stmt *stmts;
    Stmt *stmt;

    int nstmts = scoplib_statement_number(scop->statement);

    /* Allocate more to account for unroll/jamming later on */
    stmts = (Stmt *) malloc(nstmts*sizeof(Stmt));

    scoplib_statement_p clan_stmt = scop->statement;

    for(i=0; i<nstmts; i++)  {
        stmt = &stmts[i];

        stmt->id = i;

        stmt->dim = clan_stmt->nb_iterators;

        assert(clan_stmt->domain->elt->NbColumns-1 == stmt->dim + npar + 1);

        stmt->domain = clan_matrix_to_pluto_constraints(clan_stmt->domain->elt);

        /* Initialization */
        stmt->num_tiled_loops = 0;

        for (j=0; j<nvar; j++)  {
            stmt->is_supernode[j] = false;
        }

        for (j=0; j<stmt->dim; j++)  {
            stmt->is_outer_loop[j] = true;
        }

        stmt->trans = pluto_matrix_alloc(MAX_TRANS_ROWS, 
                MAX_TILING_LEVELS*nvar+nvar+1);

        stmt->trans->nrows = 0;
        stmt->trans->ncols = nvar+1;

        stmt->num_ind_sols = 0;

        /* Tile it if it's tilable unless turned off by .fst/.precut file */
        stmt->tile = 1;

        clan_stmt = clan_stmt->next;
    }

    return stmts;
}


void stmts_print(FILE *fp, Stmt *stmts, int nstmts)
{
    int i;

    for(i=0; i<nstmts; i++)  {
        Stmt stmt = stmts[i];
        fprintf(fp, "S%d %d-d index set\n", stmt.id+1, stmt.dim);
        pluto_constraints_pretty_print(fp, stmt.domain);
    }
}


void stmt_free(Stmt *stmt)
{
    pluto_matrix_free(stmt->trans);
    pluto_constraints_free(stmt->domain);
}


void dep_free(Dep *dep)
{
    pluto_constraints_free(dep->dpolytope);
}


/* 
 * Extract necessary information from clan_scop to create PlutoProg - a
 * representation of the program sufficient to be used throughout Pluto. 
 * PlutoProg also includes dependences; so candl is run here.
 */
PlutoProg *scop_to_pluto_prog(scoplib_scop_p scop, PlutoOptions *options)
{
    PlutoProg *prog = (PlutoProg *) malloc(sizeof(PlutoProg));

    prog->nstmts = scoplib_statement_number(scop->statement);
    prog->options = options;

    /* Set global variables first (they are used in stmts_read, deps_read too) */
    prog->npar = scop->nb_parameters;
    scoplib_statement_p clan_stmt = scop->statement;

    prog->nvar = clan_stmt->nb_iterators;

    int i;
    for (i=0; i<prog->nstmts; i++) {
        prog->nvar = PLMAX(prog->nvar, clan_stmt->nb_iterators);
        clan_stmt = clan_stmt->next;
    }

    /* Calculate dependences using Candl */

    candl_program_p candl_program = candl_program_convert_scop(scop, NULL);

    CandlOptions *candlOptions = candl_options_malloc();
    if (options->rar)   {
        candlOptions->rar = 1;
    }
    candlOptions->lastwriter = options->lastwriter;
    candlOptions->scalar_privatization = options->scalpriv;
    // candlOptions->verbose = 1;


    CandlDependence *candl_deps = candl_dependence(candl_program, candlOptions);

    prog->stmts = stmts_read(scop, prog->npar, prog->nvar);
    prog->deps = deps_read(candl_deps, prog);
    prog->ndeps = candl_num_dependences(candl_deps);

    candl_options_free(candlOptions);
    candl_dependence_free(candl_deps);
    candl_program_free(candl_program);


    /* Allocate and initialize hProps */
    prog->hProps = (HyperplaneProperties *) 
        malloc(MAX_TRANS_ROWS*sizeof(HyperplaneProperties));

    for (i=0; i<MAX_TRANS_ROWS; i++)    {
        prog->hProps[i].unroll = NO_UNROLL;
    }

    /* Parameter names */
    prog->params = (char **) malloc(sizeof(char *)*prog->npar);
    for (i=0; i<prog->npar; i++)  {
        prog->params[i] = (char *) malloc(sizeof(char)*64);
        strcpy(prog->params[i], scop->parameters[i]);
    }

    /* Iterator names and statement text */
    clan_stmt = scop->statement;
    for (i=0; i<prog->nstmts; i++)    {
        prog->stmts[i].iterators = (char **) malloc(sizeof(char *)*prog->stmts[i].dim);
        int j;
        for (j=0; j<prog->stmts[i].dim; j++)    {
            prog->stmts[i].iterators[j] = (char *) malloc(sizeof(char)*64);
            strcpy(prog->stmts[i].iterators[j], clan_stmt->iterators[j]);
        }
        /* Statement text */
        prog->stmts[i].text = (char *) malloc(sizeof(char)*(strlen(clan_stmt->body)+1));
        strcpy(prog->stmts[i].text, clan_stmt->body);
        clan_stmt = clan_stmt->next;
    }

	// hack for linearized accesses
	FILE *lfp = fopen(".linearized", "r");
	FILE *nlfp = fopen(".nonlinearized", "r");
	char tmpstr[256];
	char linearized[256];
       if (lfp && nlfp) {
               for (i=0; i<prog->nstmts; i++)    {
                       rewind(lfp);
                       rewind(nlfp);
                       while (!feof(lfp))      {
                               fgets(tmpstr, 256, nlfp);
                               fgets(linearized, 256, lfp);
                               // printf("%s\n", tmpstr);
                               // printf("%s\n", prog->stmts[i].text);
                               if (strstr(tmpstr, prog->stmts[i].text))        {
                                       // printf("Found substring\n");
                                       // printf("%s\n", linearized);
                                       prog->stmts[i].text = (char *) realloc(prog->stmts[i].text, sizeof(char)*(strlen(linearized)+1));
                                       strcpy(prog->stmts[i].text, linearized);
                               }
			}
		}
               fclose(lfp);
               fclose(nlfp);
	}

    return prog;
}


void pluto_prog_free(PlutoProg *prog)
{
    int i;

    /* Free the dependences */
    for (i=0; i<prog->ndeps; i++) {
        dep_free(&prog->deps[i]);
    }
    free(prog->deps);

    /* Free the DDG */
    graph_free(prog->ddg);

    free(prog->hProps);

    for (i=0; i<prog->npar; i++)  {
        free(prog->params[i]);
    }
    free(prog->params);

    /* Iterator names and statement text */
    for (i=0; i<prog->nstmts; i++)    {
        int j;
        for (j=0; j<prog->stmts[i].dim; j++)    {
            /* TODO: increase iterators while tiling */
            // free(prog->stmts[i].iterators[j]);
        }
        free(prog->stmts[i].iterators);

        /* Statement text */
        free(prog->stmts[i].text);
    }

    /* Statements */
    for (i=0; i<prog->nstmts; i++) {
        stmt_free(&prog->stmts[i]);
    }
    free(prog->stmts);

    free(prog);
}


PlutoOptions *pluto_options_alloc()
{
    PlutoOptions *options;

    options  = (PlutoOptions *) malloc(sizeof(PlutoOptions));

    /* Initialize to default */
    options->tile = 0;
    options->debug = 0;
    options->moredebug = 0;
    options->scancount = 0;
    options->parallel = 0;
    options->unroll = 0;

    /* Unroll/jam factor */
    options->ufactor = 8;

    /* Ignore input deps */
    options->rar = 0;

    /* Override for first and last levels to tile */
    options->ft = -1;
    options->lt = -1;

    /* Override for first and last cloog options */
    options->cloogf = -1;
    options->cloogl = -1;

    options->multipipe = 0;
    options->l2tile = 0;
    options->prevector = 1;
    options->fuse = SMART_FUSE;

    /* Experimental */
    options->polyunroll = 0;

    /* Default context is no context */
    options->context = -1;

    options->bee = 0;

    options->lastwriter = 0;

    options->nobound = 0;

    options->scalpriv = 0;

    options->silent = 0;

    return options;
}


void pluto_options_free(PlutoOptions *options)
{
    free(options);
}
