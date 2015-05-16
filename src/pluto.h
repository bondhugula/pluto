/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007 Uday Bondhugula
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
 * A copy of the GNU General Public Licence can be found in the 
 * top-level directory of this program (`COPYING') 
 *
 */
#ifndef _PLUTO_H_
#define _PLUTO_H_

#include <stdbool.h>

#include "osl/scop.h"

#include "math_support.h"
#include "constraints.h"
#include "ddg.h"
#include "pluto/libpluto.h"

#include "osl/extensions/dependence.h"

/* Check out which piplib we are linking with */
/* Candl/piplib_wrapper converts relation to matrices */
#ifdef SCOPLIB_INT_T_IS_LONGLONG // Defined in src/Makefile.am
#define PLUTO_OSL_PRECISION 64
#elif  SCOPLIB_INT_T_IS_MP
#define PLUTO_OSL_PRECISION 0
#endif

#define IF_DEBUG(foo) {if (options->debug || options->moredebug) { foo; }}
#define IF_DEBUG2(foo) {if (options->moredebug) {foo; }}
#define IF_MORE_DEBUG(foo) {if (options->moredebug) {foo; }}
#define PLUTO_MESSAGE(foo) {if (!options->silent) { foo; }}

#define MAX_TILING_LEVELS 2

#define DEFAULT_L1_TILE_SIZE 32

#define CST_WIDTH (npar+1+nstmts*(nvar+1)+1)

#define ALLOW_NEGATIVE_COEFF 1 
#define DO_NOT_ALLOW_NEGATIVE_COEFF 0 

typedef enum dirvec_type {DEP_MINUS='-', DEP_ZERO='0', DEP_PLUS='+', DEP_STAR='*'} DepDir;

/* H_TILE_SPACE_LOOP may not always be distinguished from H_LOOP */
typedef enum hyptype {H_UNKNOWN=0, H_LOOP, H_TILE_SPACE_LOOP,
    H_SCALAR} PlutoHypType;

/* Candl dependences are not marked uniform/non-uniform */
#define IS_UNIFORM(type) (0)
#define IS_RAR(type) (type == OSL_DEPENDENCE_RAR)
#define IS_RAW(type) (type == OSL_DEPENDENCE_RAW)
#define IS_WAR(type) (type == OSL_DEPENDENCE_WAR)
#define IS_WAW(type) (type == OSL_DEPENDENCE_WAW)

typedef enum looptype {UNKNOWN=0, PARALLEL, PIPE_PARALLEL, SEQ, 
    PIPE_PARALLEL_INNER_PARALLEL} PlutoLoopType;


/* ORIG is an original compute statement provided by a polyhedral extractor */
typedef enum stmttype {ORIG=0, STMT_UNKNOWN} PlutoStmtType;

typedef struct pluto_access{
    int sym_id;
    char *name;

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

    /* The hyperplane evicted by the hyperplane enabling
     * concurrent start (diamond tiling) */
    PlutoMatrix *evicted_hyp;

    /* H_LOOP, H_SCALAR, .. */
    PlutoHypType *hyp_types;

    /* Num of scattering dimensions tiled */
    int num_tiled_loops;

    PlutoAccess **reads;
    int nreads;

    PlutoAccess **writes;
    int nwrites;

    /***/
    /* Used by scheduling algo */
    /***/

    /* ID of the SCC in the DDG this statement belongs to */
    int scc_id;

    int first_tile_dim;
    int last_tile_dim;

    PlutoStmtType type;

    /* Compute statement associated with distmem copy/sigma stmt */
    const struct statement *parent_compute_stmt;
};
typedef struct statement Stmt;

struct stmt_access_pair{
    PlutoAccess *acc;
    Stmt *stmt;
};


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
     *        [src iterators|dest iterators | params |const] >= 0
     * sizes: [src_dim      |dest_dim       |npar    |1    ]
     */
    PlutoConstraints *dpolytope;

    /*
     * Polyhedra used to store source unique dependence polyhedra
     * in FOP scheme
     */
    PlutoConstraints *src_unique_dpolytope;

    PlutoConstraints *depsat_poly;

    int *satvec;

    /* Dependence type from Candl (raw, war, or rar) */
    int type;

    /* Has this dependence been satisfied? */
    bool satisfied;

    /* Level at which the dependence is satisfied */
    int satisfaction_level;

    /* Constraints for preserving this dependence while bounding 
     * its distance */
    PlutoConstraints *cst;

    /* Constraints for bounding dependence distance */
    PlutoConstraints *bounding_cst;

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
    PlutoHypType type;

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
    /* Does not contain transitive dependences if options->lastwriter */
    Dep **deps;
    int ndeps;

    /* Array of dependences */
    /* Used for calculating write-out set only if options->lastwriter */
    /* May contain transitive WAR dependences */
    Dep **transdeps;
    int ntransdeps;

    /* Array of data variable names */
    // required only by commopt_foifi option
    char **data_names;
    int num_data;

    /* Parameters */
    char **params;

    /* Number of hyperplanes that represent the transformed space
     * same as stmts[i].trans->nrows, for all i; even if statements 
     * are allowed to have different number of rows, at codegen time,
     * they all have to be padded to the maximum across all statements;
     * in that case, num_hyperplanes is intended to be max across
     * stmt->trans->nrows */
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

    /* Codegen context */
    PlutoConstraints *codegen_context;
    /* Temp autotransform data */
    PlutoConstraints *globcst;

    /* Hyperplane that was replaced in case concurrent start 
     * had been found*/
    int evicted_hyp_pos;

    osl_scop_p scop;
};
typedef struct plutoProg PlutoProg;

/*
 * A Ploop is NOT an AST loop; this is a dimension in the scattering tree
 * which is not a scalar one. Ploop exists in the polyhedral representation
 * and corresponds to one or more loops at the *same* depth (same t<num>) in
 * the final generated AST
 */
typedef struct pLoop{
    int depth;
    Stmt **stmts;
    int nstmts;
} Ploop;

struct pluto_dep_list {
    Dep *dep;
    struct pluto_dep_list *next;
};
typedef struct pluto_dep_list PlutoDepList;

struct pluto_constraints_list {
    PlutoConstraints *constraints;
    PlutoConstraints *iterations;
    PlutoDepList *deps;
    Stmt *stmt;
    PlutoAccess *access;
    struct pluto_constraints_list *next;
};

typedef struct pluto_constraints_list PlutoConstraintsList;

PlutoConstraintsList* pluto_constraints_list_alloc(PlutoConstraints *cst);

void pluto_constraints_list_free(PlutoConstraintsList *cstlist);

void pluto_constraints_list_add(PlutoConstraintsList *list,const PlutoConstraints *cst,
    Dep *dep, int copyDep);

void pluto_constraints_list_replace(PlutoConstraintsList *list, PlutoConstraints *cst);

typedef struct band{
    Ploop *loop;
    int width;
    /* Not used yet */
    struct band **children;
} Band;

/* Globally visible, easily accessible data */
/* It's declared in main.c (for pluto binary) and in libpluto.c for 
 * the library */
extern PlutoOptions *options;

void dep_alloc_members(Dep *);
void dep_free(Dep *);

void pluto_compute_dep_satisfaction(PlutoProg *prog);
void pluto_compute_dep_satisfaction_complex(PlutoProg *prog);
bool dep_is_satisfied(Dep *dep);
void pluto_detect_transformation_properties(PlutoProg *prog);
void pluto_detect_hyperplane_types_stmtwise(PlutoProg *prog);

void pluto_compute_satisfaction_vectors(PlutoProg *prog);
void pluto_compute_dep_directions(PlutoProg *prog);

PlutoConstraints *get_permutability_constraints(PlutoProg *);
PlutoConstraints **get_stmt_ortho_constraints(Stmt *stmt, const PlutoProg *prog,
        const PlutoConstraints *currcst, int *orthonum);
PlutoConstraints *get_global_independence_cst(
        PlutoConstraints ***ortho_cst, int *orthonum, 
        const PlutoProg *prog);
PlutoConstraints *get_non_trivial_sol_constraints(const PlutoProg *, bool);

int pluto_auto_transform(PlutoProg *prog);
int  pluto_multicore_codegen(FILE *fp, FILE *outfp, const PlutoProg *prog);

int  find_permutable_hyperplanes(PlutoProg *prog, bool lin_ind_mode, 
       bool loop_search_mode, int max_sols);

void detect_hyperplane_type(Stmt *stmts, int nstmts, Dep *deps, int ndeps, int, int, int);
DepDir  get_dep_direction(const Dep *dep, const PlutoProg *prog, int level);

void getInnermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd);
void getOutermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd);

void pluto_gen_cloog_file(FILE *fp, const PlutoProg *prog);
void cut_lightest_edge(Stmt *stmts, int nstmts, Dep *deps, int ndeps, int);
void pluto_tile(PlutoProg *);
bool pluto_create_tile_schedule(PlutoProg *prog, Band **bands, int nbands);
int pluto_detect_mark_unrollable_loops(PlutoProg *prog);

int pluto_omp_parallelize(PlutoProg *prog);

void   ddg_update(Graph *g, PlutoProg *prog);
void   ddg_compute_scc(PlutoProg *prog);
Graph *ddg_create(PlutoProg *prog);
int    ddg_sccs_direct_conn(Graph *g, PlutoProg *prog, int scc1, int scc2);

void unroll_phis(PlutoProg *prog, int unroll_dim, int ufactor);

void pluto_print_dep_directions(PlutoProg *prog);
void pluto_print_depsat_vectors(PlutoProg *prog, int levels);
PlutoConstraints *pluto_stmt_get_schedule(const Stmt *stmt);
void pluto_update_deps(Stmt *stmt, PlutoConstraints *cst, PlutoProg *prog);

PlutoMatrix *get_new_access_func(const Stmt *stmt, const PlutoMatrix *acc, const PlutoProg *prog);
PlutoConstraints *pluto_get_new_domain(const Stmt *stmt);
PlutoConstraints *pluto_compute_region_data(const Stmt *stmt, const PlutoConstraints *dom,
        const PlutoAccess *acc, int copy_level, const PlutoProg *prog);

int generate_declarations(const PlutoProg *prog, FILE *outfp);
int pluto_gen_cloog_code(const PlutoProg *prog, int cloogf, int cloogl, FILE *cloogfp, FILE *outfp);
void pluto_add_given_stmt(PlutoProg *prog, Stmt *stmt);

int is_loop_dominated(Ploop *loop1, Ploop *loop2, const PlutoProg *prog);
Ploop **pluto_get_parallel_loops(const PlutoProg *prog, int *nploops);
Ploop **pluto_get_all_loops(const PlutoProg *prog, int *num);
Ploop **pluto_get_dom_parallel_loops(const PlutoProg *prog, int *nploops);
Band **pluto_get_dom_parallel_bands(PlutoProg *prog, int *nbands, int **comm_placement_levels);
void pluto_loop_print(const Ploop *loop);
void pluto_loops_print(Ploop **loops, int num);
void pluto_loops_free(Ploop **loops, int nloops);
int pluto_loop_compar(const void *_l1, const void *_l2);
Band *pluto_band_alloc(Ploop *loop, int width);
void pluto_bands_print(Band **bands, int num);
void pluto_band_print(const Band *band);

Band **pluto_get_outermost_permutable_bands(PlutoProg *prog, int *ndbands);
Ploop *pluto_loop_dup(Ploop *l);
int pluto_loop_is_parallel(const PlutoProg *prog, Ploop *loop);
int pluto_loop_is_parallel_for_stmt(const PlutoProg *prog, const Ploop *loop, 
        const Stmt *stmt);
int pluto_loop_has_satisfied_dep_with_component(const PlutoProg *prog, 
        const Ploop *loop);
void pluto_bands_free(Band **bands, int nbands);
int pluto_is_hyperplane_loop(const Stmt *stmt, int level);
void pluto_detect_hyperplane_types(PlutoProg *prog);
void pluto_tile_band(PlutoProg *prog, Band *band, int *tile_sizes);

Ploop **pluto_get_loops_under(Stmt **stmts, int nstmts, int depth,
        const PlutoProg *prog, int *num);
Ploop **pluto_get_loops_immediately_inner(Ploop *ploop, PlutoProg *prog, int *num);
int pluto_intra_tile_optimize(PlutoProg *prog,  int is_tiled);
int pluto_intra_tile_optimize_band(Band *band, int is_tiled, PlutoProg *prog);

int pluto_pre_vectorize_band(Band *band, int is_tiled, PlutoProg *prog);
int pluto_is_band_innermost(const Band *band, int is_tiled);
Band **pluto_get_innermost_permutable_bands(PlutoProg *prog, int *ndbands);
int pluto_is_loop_innermost(const Ploop *loop, const PlutoProg *prog);

PlutoConstraints *pluto_get_transformed_dpoly(const Dep *dep, Stmt *src, Stmt *dest);

void pluto_detect_scalar_dimensions(PlutoProg *prog);
int pluto_detect_mark_unrollable_loops(PlutoProg *prog);
int pluto_are_stmts_fused(Stmt **stmts, int nstmts, const PlutoProg *prog);

void pluto_iss_dep(PlutoProg *prog);
PlutoConstraints *pluto_find_iss(const PlutoConstraints **doms, int ndoms, int npar, PlutoConstraints *);
void pluto_iss(Stmt *stmt, PlutoConstraints **cuts, int num_cuts, PlutoProg *prog);


#endif
