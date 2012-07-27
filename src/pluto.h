/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007 Uday Kumar Bondhugula
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * A copy of the GNU General Public Licence can be found in the 
 * top-level directory of this program (`COPYING') 
 *
 */
#ifndef _PLUTO_H_
#define _PLUTO_H_

#include <stdbool.h>

#include "scoplib/symbol.h"

#include "math_support.h"
#include "constraints.h"
#include "ddg.h"
#include "pluto/libpluto.h"

#define IF_DEBUG(foo) {if (options->debug || options->moredebug) { foo; }}
#define IF_DEBUG2(foo) {if (options->moredebug) {foo; }}

/* Do not fuse across SCCs */
#define NO_FUSE 0
/* Geared towards maximal fusion, but not really maximal fusion */
#define MAXIMAL_FUSE 1
/* Something in between the above two */
#define SMART_FUSE 2

#define MAX_CONSTRAINTS 10000
#define MAX_FARKAS_CST  2000

#define MAX_TILING_LEVELS 2

#define DEFAULT_L1_TILE_SIZE 32

#define CST_WIDTH (npar+1+nstmts*(nvar+1)+1)

typedef enum dirvec_type {DEP_MINUS='-', DEP_ZERO='0', DEP_PLUS='+', DEP_STAR='*'} DepDir;

#define H_UNKNOWN 0
#define H_LOOP 1
#define H_TILE_SPACE_LOOP 2
#define H_SCALAR 3

/* Candl dependences are not marked uniform/non-uniform */
#define IS_UNIFORM(type) (0)
#define IS_RAR(type) (type == CANDL_RAR)


typedef enum looptype {UNKNOWN=0, PARALLEL, PIPE_PARALLEL, SEQ, 
    PIPE_PARALLEL_INNER_PARALLEL} PlutoLoopType;

typedef struct pluto_access{
    int sym_id;
    char *name;

    scoplib_symbol_p symbol;

    PlutoMatrix *mat;
} PlutoAccess;


struct statement{
    int id;

    PlutoConstraints *domain;

    /* Original iterator names */
    char **iterators;

    /* Statement text */
    char *text;

    /* Does this dimension appear in the statement's original domain? 
     * dummy dimensions added (sinking) will have is_orig_loop as false
     */
    bool *is_orig_loop;

    /* Dimensionality of statement's domain */
    int dim;

    /* Original dimensionality of statement's domain */
    int dim_orig;

    /* Should you tile even if it's tilable? */
    int tile;

    /* Affine transformation matrix that will be completed step by step */
    /* this captures the A/B/G notation of the INRIA group - except that
     * we don't have parametric shifts
     * Size: nvar+npar+1 columns inside auto_transform
     *       stmt->dim + npar + 1 after auto_transform (after denormalize)
     */
    PlutoMatrix *trans;

    /* Num of scattering dimensions tiled */
    int num_tiled_loops;

    PlutoAccess **reads;
    int nreads;

    PlutoAccess **writes;
    int nwrites;


    /***/
    /* Used by scheduling algo */
    /***/

    /* Num of independent soln's found */
    int num_ind_sols;

    /* ID of the SCC in the DDG this statement belongs to */
    int scc_id;

};
typedef struct statement Stmt;


struct dependence{

    /* Unique number of the dependence: starts with 0  */
    int id;

    /* Source statement ID -- can be used to directly index into Stmt *stmts */
    int src;

    /* Dest statement ID */
    int dest;

    /* Points into statements's read or write access */
    PlutoAccess *src_acc;
    PlutoAccess *dest_acc;

    /* 
     * Dependence polyhedra (both src & dest) 
     * (product space)
     *
     * [src|dest|par|const] >= 0
     * [nvar|nvar|npar|1]
     */
    PlutoConstraints *dpolytope;

    /* Dependence type from Candl (raw, war, or rar) */
    int type;

    /* Has this dependence been satisfied? */
    bool satisfied;

    /* Level at which the dependence is satisfied */
    int satisfaction_level;

    /* Dependence direction in transformed space */
    DepDir *dirvec;
};
typedef struct dependence Dep;


typedef enum unrollType {NO_UNROLL, UNROLL, UNROLLJAM} UnrollType;


/* Properties of the new hyperplanes found. These are common across all
 * statements or apply at a level across all statements 
 */
struct hyperplane_properties{

    /* Hyperplane property: see looptype enum definition */
    PlutoLoopType dep_prop;

    /* Hyperplane type: scalar, loop, or tile-space loop (H_SCALAR,
     * H_LOOP, H_TILE... */
    int type;

    /* The band number this hyperplane belongs to. Note that everything is a
     * hierarchy of permutable loop nests (it's not a tree, but a straight
     * hierarchy) */
    int band_num;

    /* Unroll or Unroll-jam this dimension? */
    UnrollType unroll;

    /* Mark for icc vectorization */
    int prevec;
};
typedef struct hyperplane_properties HyperplaneProperties;

struct plutoProg{
    /* Array of statements */
    Stmt **stmts;
    int nstmts;

    /* Array of dependences */
    Dep **deps;
    int ndeps;

    /* Parameters */
    char **params;

    /* Number of hyperplanes that represent the transformed space
     * same as stmts[i].trans->nrows, for all i */
    int num_hyperplanes;

    /* Data dependence graph of the program */
    Graph *ddg;

    /* Options for Pluto */
    PlutoOptions *options;

    /* Hyperplane properties */
    HyperplaneProperties *hProps;

    /* Max (original) domain dimensionality */
    int nvar;

    /* Number of program parameters */
    int npar;

    /* Param context */
    PlutoConstraints *context;

};
typedef struct plutoProg PlutoProg;


/* Globally visible, easily accessible data */
extern PlutoOptions *options;

void dep_alloc_members(Dep *);
void dep_free(Dep *);

void pluto_compute_dep_satisfaction(PlutoProg *prog);
bool dep_is_satisfied(Dep *dep);
void pluto_detect_transformation_properties(PlutoProg *prog);

PlutoConstraints *get_permutability_constraints(Dep **, int, const PlutoProg *);
PlutoConstraints **get_stmt_ortho_constraints(Stmt *stmt, const PlutoProg *prog,
        const PlutoConstraints *currcst,
        int *orthonum);
PlutoConstraints *get_non_trivial_sol_constraints(const PlutoProg *);

int pluto_auto_transform(PlutoProg *prog);
int  pluto_multicore_codegen(FILE *fp, FILE *outfp, const PlutoProg *prog);

int  find_permutable_hyperplanes(PlutoProg *prog, int max_sols);
void detect_hyperplane_type(Stmt *stmts, int nstmts, Dep *deps, int ndeps, int, int, int);
int  get_dep_direction(const Dep *dep, const PlutoProg *prog, int level);

void getInnermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd);
void getOutermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd);

void pluto_gen_cloog_file(FILE *fp, const PlutoProg *prog);
void cut_lightest_edge(Stmt *stmts, int nstmts, Dep *deps, int ndeps, int);
void pluto_tile(PlutoProg *);
bool create_tile_schedule(PlutoProg *prog, int firstD, int lastD);

int pluto_omp_parallelize(PlutoProg *prog);

void   ddg_update(Graph *g, PlutoProg *prog);
void   ddg_compute_scc(PlutoProg *prog);
Graph *ddg_create(PlutoProg *prog);
int    ddg_sccs_direct_conn(Graph *g, PlutoProg *prog, int scc1, int scc2);

void unroll_phis(PlutoProg *prog, int unroll_dim, int ufactor);

void pretty_print_affine_function(FILE *fp, const Stmt *stmt, int level);
void pluto_transformations_print(const PlutoProg *prog);
void pluto_transformations_pretty_print(const PlutoProg *prog);
void pluto_print_hyperplane_properties(const PlutoProg *prog);
void pluto_print_dep_directions(Dep **deps, int ndeps, int levels);
PlutoConstraints *pluto_stmt_get_schedule(const Stmt *stmt);
void pluto_update_deps(Stmt *stmt, PlutoConstraints *cst, PlutoProg *prog);

int generate_declarations(const PlutoProg *prog, FILE *outfp);
int pluto_gen_cloog_code(const PlutoProg *prog, FILE *cloogfp, FILE *outfp);

int gen_vecloop_file(PlutoProg *prog);

#endif
