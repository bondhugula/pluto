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

#include "constraints.h"
#include "ddg.h"
#include "math_support.h"
#include "pluto/libpluto.h"

#include "osl/extensions/dependence.h"

#ifdef GLPK
#include <glpk.h>
#endif

/* Check out which piplib we are linking with */
/* Candl/piplib_wrapper converts relation to matrices */
#ifdef SCOPLIB_INT_T_IS_LONGLONG // Defined in src/Makefile.am
#define PLUTO_OSL_PRECISION 64
#elif SCOPLIB_INT_T_IS_MP
#define PLUTO_OSL_PRECISION 0
#endif

#define IF_DEBUG(foo)                                                          \
  {                                                                            \
    if (options->debug || options->moredebug) {                                \
      foo;                                                                     \
    }                                                                          \
  }
#define IF_DEBUG2(foo)                                                         \
  {                                                                            \
    if (options->moredebug) {                                                  \
      foo;                                                                     \
    }                                                                          \
  }
#define IF_MORE_DEBUG(foo)                                                     \
  {                                                                            \
    if (options->moredebug) {                                                  \
      foo;                                                                     \
    }                                                                          \
  }
#define PLUTO_MESSAGE(foo)                                                     \
  {                                                                            \
    if (!options->silent) {                                                    \
      foo;                                                                     \
    }                                                                          \
  }

#define MAX_TILING_LEVELS 2

#define DEFAULT_L1_TILE_SIZE 32

#define PI_TABLE_SIZE 256

#define CST_WIDTH (npar + 1 + nstmts * (nvar + 1) + 1)

#define ALLOW_NEGATIVE_COEFF 1
#define DO_NOT_ALLOW_NEGATIVE_COEFF 0

/* Iterative search modes */
#define EAGER 0
#define LAZY 1

typedef enum dirvec_type {
  DEP_MINUS = '-',
  DEP_ZERO = '0',
  DEP_PLUS = '+',
  DEP_STAR = '*'
} DepDir;

/* H_TILE_SPACE_LOOP may not always be distinguished from H_LOOP */
typedef enum hyptype {
  H_UNKNOWN = 0,
  H_LOOP,
  H_TILE_SPACE_LOOP,
  H_SCALAR
} PlutoHypType;

/* Candl dependences are not marked uniform/non-uniform */
#define IS_UNIFORM(type) (0)
#define IS_RAR(type) (type == OSL_DEPENDENCE_RAR)
#define IS_RAW(type) (type == OSL_DEPENDENCE_RAW)
#define IS_WAR(type) (type == OSL_DEPENDENCE_WAR)
#define IS_WAW(type) (type == OSL_DEPENDENCE_WAW)

typedef enum looptype {
  UNKNOWN = 0,
  PARALLEL,
  PIPE_PARALLEL,
  SEQ,
  PIPE_PARALLEL_INNER_PARALLEL
} PlutoLoopType;

/* ORIG is an original compute statement provided by a polyhedral extractor */
typedef enum stmttype {
  ORIG = 0,
  ORIG_IN_FUNCTION,
  IN_FUNCTION,
  COMM_CALL,
  LW_COMM_CALL,
  SIGMA,
  TAU,
  MISC,
  DATA_DIST_INIT,
  DATA_SETUP,
  DATA_DIST_MANG,
  DATA_DIST_COPY,
  COPY_OUT,
  COPY_IN,
  LW_COPY_OUT,
  LW_COPY_IN,
  ALL_TASKS,
  STMT_UNKNOWN
} PlutoStmtType;

typedef struct pluto_access {
  int sym_id;
  char *name;

  /* scoplib_symbol_p symbol; */

  PlutoMatrix *mat;
} PlutoAccess;

struct pet_stmt;

struct statement {
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
  unsigned dim;

  /* Original dimensionality of statement's domain */
  unsigned dim_orig;

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

  /* Hyperplane that was replaced in case concurrent start
   * had been found*/
  int evicted_hyp_pos;

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

  /* ID of the CC in the DDG this statement belongs to */
  int cc_id;

  int first_tile_dim;
  int last_tile_dim;

  PlutoStmtType type;

  /* ID of the domain parallel loop that the statement belongs to */
  int ploop_id;

  /* indicates whether inner loop is vectorized or not*/
  int inner_loop_vec;

  int pos_peel_offset;

  int neg_peel_offset;

  /* Compute statement associated with distmem copy/sigma stmt */
  const struct statement *parent_compute_stmt;

  /* Intra statement depndence constraints.  Used to construct the fusion
   * conflict graph */
  PlutoConstraints *intra_stmt_dep_cst;

  struct pet_stmt *pstmt;
};
typedef struct statement Stmt;

struct stmt_access_pair {
  PlutoAccess *acc;
  Stmt *stmt;
};

struct TiledHyperplane {

  // Is the corresponding loop hyperplane tiled?
  int loop_tiled;

  // already tiled?
  int is_tiled;

  // Used to store the tile size tiled dim
  int tile_sizes;

  // Indicates the tile loop dim corresponding to original array dim
  int tiled_dim;

  // current col number of orignal dim
  int orig_dim;

  int depth;

  int firstD;

  int lastD;

  char *size_buffer;
};

struct array {

  int id;

  /*name of the array*/
  char *text;

  // datatype
  char *data_type;

  PlutoConstraints *array_bounds;

  PlutoConstraints **parmetric_domain;

  int dim_updated;

  /*Dimensionality of array*/
  int dim;

  /*Original dimensionality of array*/
  int dim_orig;

  //	Should the data layout of the array be changed?
  int change_layout;

  // Affine data transformation matrix of array
  // Size: nvar+npar+1
  PlutoMatrix *trans;
  PlutoMatrix *trans_orig;

  // Read access of the array
  PlutoAccess **reads;
  int nreads;

  // Write access of the array
  PlutoAccess **writes;
  int nwrites;

  // Statements that have this array reference
  int nstmts;
  Stmt **stmts;

  // Num of hyperplanes found
  int num_hyperplanes_found;

  // iterator names
  char **iterators;

  struct TiledHyperplane *tiled_hyperplane;

  int num_tiled_loops;

  int first_tile_dim;

  int last_tile_dim;

  int num_cur_tiled_loops;

  // dim upto copy level should be treated as parameters
  int *copy_level;

  int copy_level_used;

  // npar in arr parametric domain
  int npar;

  int *hyperplane_mapping;
};
typedef struct array Array;

struct dependence {

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

  /* Used for bounding function constraints */
  PlutoConstraints *bounding_poly;

  /*
   * Polyhedra used to store source unique dependence polyhedra
   * in FOP scheme
   */
  PlutoConstraints *src_unique_dpolytope;

  PlutoConstraints *depsat_poly;

  /* Dependence type from Candl (raw, war, or rar) */
  int type;

  /* Has this dependence been satisfied? */
  bool satisfied;

  /* Does this dependence need to be skipped? Set to true during variable
   * liberalization */
  bool skipdep;

  /* Level at which this dependence is completely satisfied (when doing
   * conservative computation) or level *by* which the dependence is
   * completely satisfied (if doing a complex/accurate check) */
  int satisfaction_level;

  /*
   * A vector which provides the levels at which a
   * dependence has been satisfied. Depending on how this is
   * computed, it may be conservative or accurate
   */
  int *satvec;

  /* Constraints for preserving this dependence while bounding
   * its distance */
  PlutoConstraints *cst;

  /* Constraints for bounding dependence distance */
  PlutoConstraints *bounding_cst;

  /* Dependence direction in transformed space */
  DepDir *dirvec;
};
typedef struct dependence Dep;

typedef enum unrollType { NO_UNROLL, UNROLL, UNROLLJAM } UnrollType;

/* Properties of the new hyperplanes found. These are common across all
 * statements or apply at a level across all statements
 */
struct hyperplane_properties {

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

struct plutoProg {
  /* Array of statements */
  Stmt **stmts;
  int nstmts;

  /* Array of dependences */
  /* Does not contain transitive dependences if options->lastwriter */
  Dep **deps;
  int ndeps;

  // Array of arrays in the prog
  Array **arrays;
  int narrays;

  /* Array of dependences */
  /* Used for calculating write-out set only if options->lastwriter */
  /* May contain transitive WAR dependences */
  Dep **transdeps;
  int ntransdeps;

  /* Array of data variable names */
  /* used by distmem to set different tags for different data */
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

  /* Fusion conflict graph of the program */
  Graph *fcg;

  /* Options for Pluto */
  PlutoOptions *options;

  /* Hyperplane properties */
  HyperplaneProperties *hProps;

  /* Max (original) domain dimensionality */
  int nvar;

  /* Max dim of all the arrays*/
  int max_array_dim;

  /* Number of program parameters */
  int npar;

  /* Param context */
  PlutoConstraints *context;

  // Declarations to be emitted.
  char *decls;

  /* Codegen context */
  PlutoConstraints *codegen_context;

  PlutoAccess **orig_accesses;
  int num_accesses;

  /* Temp autotransform data */
  PlutoConstraints *globcst;

  /* Hyperplane that was replaced in case concurrent start
   * had been found*/
  int evicted_hyp_pos;

  osl_scop_p scop;

  /* number of outermost parallel dimensions to be parameterized */
  /* used by dynschedule */
  int num_parameterized_loops;

  /* used to store decls var for vec optimization in data dist*/
  char **decl_indices;
  int num_decl_indices;
  char *vec_decls;
  int num_stmts_to_be_coloured;
  // Boolean Array indicating whether a dimension is scaled
  int *scaled_dims;

  /* Total number of statements coloured per dimension in the FCG */
  int *total_coloured_stmts;

  /* Total number of coloured dimensions that have been coloured for all
   * statements */
  int coloured_dims;

  /* Used to store constraint solving times */
  double mipTime, ilpTime, cst_solve_time, cst_const_time, cst_write_time,
      scaling_cst_sol_time, skew_time;
  double fcg_const_time, fcg_colour_time, fcg_dims_scale_time, fcg_update_time,
      fcg_cst_alloc_time;
  long int num_lp_calls;
};
typedef struct plutoProg PlutoProg;

/*
 * A Ploop is NOT an AST loop; this is a dimension in the scattering tree
 * which is not a scalar one. Ploop exists in the polyhedral representation
 * and corresponds to one or more loops at the *same* depth (same t<num>) in
 * the final generated AST
 */
typedef struct pLoop {
  unsigned depth;
  Stmt **stmts;
  unsigned nstmts;
} Ploop;

struct pluto_dep_list {
  Dep *dep;
  struct pluto_dep_list *next;
};
typedef struct pluto_dep_list PlutoDepList;

#if defined(__cplusplus)
extern "C" {
#endif

PlutoDepList *pluto_dep_list_alloc(Dep *dep);

void pluto_deps_list_free(PlutoDepList *deplist);

PlutoDepList *pluto_deps_list_dup(PlutoDepList *src);

void pluto_deps_list_append(PlutoDepList *list, Dep *dep);

struct pluto_dep_list_list {
  PlutoDepList *dep_list;
  struct pluto_dep_list_list *next;
};
typedef struct pluto_dep_list_list PlutoDepListList;

PlutoDepListList *pluto_dep_list_list_alloc();

struct pluto_constraints_list {
  PlutoConstraints *constraints;
  PlutoDepList *deps;
  struct pluto_constraints_list *next;
};

typedef struct pluto_constraints_list PlutoConstraintsList;

PlutoConstraintsList *pluto_constraints_list_alloc(PlutoConstraints *cst);

void pluto_constraints_list_free(PlutoConstraintsList *cstlist);

void pluto_constraints_list_add(PlutoConstraintsList *list,
                                const PlutoConstraints *cst, Dep *dep,
                                int copyDep);

void pluto_constraints_list_replace(PlutoConstraintsList *list,
                                    PlutoConstraints *cst);

typedef struct band {
  /* Root loop of this band */
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
int pluto_compute_dep_satisfaction_precise(PlutoProg *prog);
bool dep_is_satisfied(Dep *dep);
void pluto_detect_transformation_properties(PlutoProg *prog);
void pluto_detect_hyperplane_types_stmtwise(PlutoProg *prog);

void pluto_compute_satisfaction_vectors(PlutoProg *prog);
void pluto_compute_dep_directions(PlutoProg *prog);
void pluto_dep_satisfaction_reset(PlutoProg *prog);

void compute_pairwise_permutability(Dep *dep, PlutoProg *prog);
PlutoConstraints *get_permutability_constraints(PlutoProg *);
PlutoConstraints *get_scc_permutability_constraints(int, PlutoProg *);
PlutoConstraints *get_feautrier_schedule_constraints(PlutoProg *prog, Stmt **,
                                                     int);
PlutoConstraints **get_stmt_ortho_constraints(Stmt *stmt, const PlutoProg *prog,
                                              const PlutoConstraints *currcst,
                                              int *orthonum);
PlutoConstraints **get_array_ortho_constraints(Array *array,
                                               const PlutoProg *prog,
                                               const PlutoConstraints *currcst,
                                               int *orthonum);
PlutoConstraints *get_global_independence_cst(PlutoConstraints ***ortho_cst,
                                              int *orthonum,
                                              const PlutoProg *prog);
PlutoConstraints *get_non_trivial_sol_constraints(const PlutoProg *, bool);
PlutoConstraints *get_coeff_bounding_constraints(const PlutoProg *);

int64 *pluto_prog_constraints_lexmin(PlutoConstraints *cst, PlutoProg *prog);

int pluto_auto_transform(PlutoProg *prog);
int pluto_dist_auto_transform(PlutoProg *prog);

PlutoConstraints *pluto_dist_add_constraints_array_access(PlutoAccess *access,
                                                          Stmt *stmt,
                                                          PlutoProg *prog);

int pluto_multicore_codegen(FILE *fp, FILE *outfp, PlutoProg *prog);
int pluto_dynschedule_graph_codegen(PlutoProg *prog, FILE *sigmafp, FILE *outfp,
                                    FILE *headerfp);
int pluto_dynschedule_codegen(PlutoProg *prog, FILE *sigmafp, FILE *outfp,
                              FILE *headerfp);
int pluto_distmem_codegen(PlutoProg *prog, FILE *cloogfp, FILE *sigmafp,
                          FILE *outfp, FILE *headerfp);

int find_permutable_hyperplanes(PlutoProg *prog, bool lin_ind_mode,
                                int max_sols, int band_depth);

void detect_hyperplane_type(Stmt *stmts, int nstmts, Dep *deps, int ndeps, int,
                            int, int);
DepDir get_dep_direction(const Dep *dep, const PlutoProg *prog, int level);

void getInnermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd);
void getOutermostTilableBand(PlutoProg *prog, int *bandStart, int *bandEnd);

void pluto_gen_cloog_file(FILE *fp, const PlutoProg *prog);
void cut_lightest_edge(Stmt *stmts, int nstmts, Dep *deps, int ndeps, int);
void pluto_tile(PlutoProg *);
bool pluto_create_tile_schedule(PlutoProg *prog, Band **bands, int nbands);
int pluto_detect_mark_unrollable_loops(PlutoProg *prog);

int pluto_omp_parallelize(PlutoProg *prog);
int pluto_dynschedule_graph_parallelize(PlutoProg *prog, FILE *sigmafp,
                                        FILE *headerfp);
int pluto_dynschedule_parallelize(PlutoProg *prog, FILE *sigmafp,
                                  FILE *headerfp, FILE *pifp);
int pluto_distmem_parallelize(PlutoProg *prog, FILE *sigmafp, FILE *headerfp,
                              FILE *pifp);

void ddg_update(Graph *g, PlutoProg *prog);
void ddg_compute_scc(PlutoProg *prog);
void ddg_compute_cc(PlutoProg *prog);
Graph *ddg_create(PlutoProg *prog);
int ddg_sccs_direct_connected(Graph *g, PlutoProg *prog, int scc1, int scc2);
int cut_between_sccs(PlutoProg *prog, Graph *ddg, int scc1, int scc2);
int cut_all_sccs(PlutoProg *prog, Graph *ddg);
void cut_smart(PlutoProg *prog, Graph *ddg);

void pluto_print_dep_directions(PlutoProg *prog);
void pluto_print_depsat_vectors(PlutoProg *prog, int levels);
PlutoConstraints *pluto_stmt_get_schedule(const Stmt *stmt);

int is_access_scalar(PlutoAccess *access);

PlutoConstraints *pluto_arr_get_schedule(const Array *arr);
void pluto_update_deps(Stmt *stmt, PlutoConstraints *cst, PlutoProg *prog);
int dep_satisfaction_update(PlutoProg *prog, int level);
int deps_satisfaction_check(PlutoProg *prog);

PlutoMatrix *get_new_access_func(const Stmt *stmt, const PlutoMatrix *acc,
                                 const PlutoProg *prog);
PlutoConstraints *pluto_get_new_domain(const Stmt *stmt);

PlutoConstraints *pluto_get_new_arr_domain(const Array *arr,
                                           PlutoConstraints *domain,
                                           int copy_level);

PlutoConstraints *pluto_compute_region_data(const Stmt *stmt,
                                            const PlutoConstraints *dom,
                                            const PlutoAccess *acc,
                                            int copy_level,
                                            const PlutoProg *prog);

int generate_declarations(PlutoProg *prog, FILE *outfp);
int pluto_gen_cloog_code(const PlutoProg *prog, int cloogf, int cloogl,
                         FILE *cloogfp, FILE *outfp);
void pluto_mark_vec_stmts(FILE *cloogfp, PlutoProg *prog);
Stmt *create_helper_stmt(const Stmt *stmt, unsigned level, const char *text,
                         PlutoStmtType type, int ploop_num);
void pluto_add_given_stmt(PlutoProg *prog, Stmt *stmt);
Stmt *pluto_create_stmt(int dim, const PlutoConstraints *domain,
                        const PlutoMatrix *trans, char **iterators,
                        const char *text, PlutoStmtType type);

void pluto_add_stmt_at_given_pos(PlutoProg *prog, Stmt *stmt, int pos);
PlutoConstraints *compute_flow_in_of_dep(Dep *dep, int copy_level,
                                         PlutoProg *prog,
                                         int use_src_unique_dpolytope);
PlutoConstraints *compute_flow_in(struct stmt_access_pair *wacc_stmt,
                                  int copy_level, PlutoProg *prog);
PlutoConstraints *compute_flow_out_of_dep(Dep *dep, int src_copy_level,
                                          int *copy_level, PlutoProg *prog,
                                          int split, PlutoConstraints **dcst1,
                                          int *pi_mappings);
PlutoConstraints *compute_flow_out(struct stmt_access_pair *wacc_stmt,
                                   int src_copy_level, int *copy_level,
                                   PlutoProg *prog, int *pi_mappings);
void compute_flow_out_partitions(struct stmt_access_pair *wacc_stmt,
                                 int src_copy_level, int *copy_level,
                                 PlutoProg *prog,
                                 PlutoConstraintsList *atomic_flowouts,
                                 int *pi_mappings);
PlutoConstraints *compute_write_out(struct stmt_access_pair *wacc_stmt,
                                    int copy_level, PlutoProg *prog);
void split_deps_acc_flowout(PlutoConstraintsList *atomic_flowouts,
                            int copy_level, int access_nrows, PlutoProg *prog);

void generate_pi(FILE *pifp, FILE *headerfp, int outer_dist_loop_level,
                 int inner_dist_loop_level, const PlutoProg *prog, Ploop *loop,
                 int loop_num, int *dimensions_to_skip, int num_dimensions);
void generate_outgoing(Stmt **loop_stmts, int nstmts, int *copy_level,
                       PlutoProg *prog, int loop_num, int *pi_mappings,
                       char *tasks_loops_decl, FILE *outfp, FILE *headerfp);
void generate_incoming(Stmt **loop_stmts, int nstmts, int *copy_level,
                       PlutoProg *prog, int loop_num, int *pi_mappings,
                       FILE *outfp, FILE *headerfp);
void generate_sigma(struct stmt_access_pair **wacc_stmts, int naccs,
                    int *copy_level, PlutoProg *prog, int loop_num,
                    int *pi_mappings, FILE *outfp, FILE *headerfp);
void generate_sigma_dep_split(struct stmt_access_pair **wacc_stmts, int naccs,
                              int *copy_level, PlutoProg *prog,
                              PlutoConstraintsList *list, int loop_num,
                              int *pi_mappings, FILE *outfp, FILE *headerfp);
void generate_tau(struct stmt_access_pair **racc_stmts, int naccs,
                  int *copy_level, PlutoProg *prog, int loop_num,
                  int *pi_mappings, FILE *outfp, FILE *headerfp);
PlutoConstraints *get_receiver_tiles_of_dep(Dep *dep, int src_copy_level,
                                            int dest_copy_level,
                                            PlutoProg *prog,
                                            int use_src_unique_dpolytope);
PlutoConstraints *get_receiver_tiles(struct stmt_access_pair **wacc_stmt,
                                     int naccs, int src_copy_level,
                                     int dest_copy_level, PlutoProg *prog,
                                     int dep_loop_num, int *pi_mappings);

void print_dynsched_file(char *srcFileName, FILE *cloogfp, FILE *outfp,
                         PlutoProg *prog);
int get_outermost_parallel_loop(const PlutoProg *prog);

int is_loop_dominated(Ploop *loop1, Ploop *loop2, const PlutoProg *prog);
Ploop **pluto_get_parallel_loops(const PlutoProg *prog, unsigned *nploops);
Ploop **pluto_get_all_loops(const PlutoProg *prog, unsigned *num);
Ploop **pluto_get_dom_parallel_loops(const PlutoProg *prog, unsigned *nploops);
Band **pluto_get_dom_parallel_bands(PlutoProg *prog, unsigned *nbands,
                                    int **comm_placement_levels);
void pluto_loop_print(const Ploop *loop);
void pluto_loops_print(Ploop **loops, int num);
void pluto_loops_free(Ploop **loops, int nloops);
int pluto_loop_compar(const void *_l1, const void *_l2);
Band *pluto_band_alloc(Ploop *loop, int width);
void pluto_bands_print(Band **bands, int num);
void pluto_band_print(const Band *band);

Band **pluto_get_outermost_permutable_bands(PlutoProg *prog, unsigned *ndbands);
Ploop *pluto_loop_dup(Ploop *l);
int pluto_loop_is_parallel(const PlutoProg *prog, Ploop *loop);
int pluto_loop_is_parallel_for_stmt(const PlutoProg *prog, const Ploop *loop,
                                    const Stmt *stmt);
int pluto_loop_satisfies_inter_stmt_dep(const PlutoProg *prog,
                                        const Ploop *loop);
int pluto_loop_has_satisfied_dep_with_component(const PlutoProg *prog,
                                                const Ploop *loop);
void pluto_bands_free(Band **bands, unsigned nbands);
int pluto_is_hyperplane_loop(const Stmt *stmt, int level);
void pluto_detect_hyperplane_types(PlutoProg *prog);
void pluto_tile_band(PlutoProg *prog, Band *band, int *tile_sizes);

Ploop **pluto_get_loops_under(Stmt **stmts, unsigned nstmts, unsigned depth,
                              const PlutoProg *prog, unsigned *num);
Ploop **pluto_get_loops_immediately_inner(Ploop *ploop, PlutoProg *prog,
                                          unsigned *num);
int pluto_intra_tile_optimize(PlutoProg *prog, int is_tiled);
int pluto_intra_tile_optimize_band(Band *band, int is_tiled, PlutoProg *prog);

int pluto_is_band_innermost(const Band *band, int is_tiled);
Band **pluto_get_innermost_permutable_bands(PlutoProg *prog, unsigned *ndbands);
int pluto_loop_is_innermost(const Ploop *loop, const PlutoProg *prog);

PlutoConstraints *pluto_get_transformed_dpoly(const Dep *dep, Stmt *src,
                                              Stmt *dest);

PlutoConstraints *pluto_find_dependence(PlutoConstraints *domain1,
                                        PlutoConstraints *domain2, Dep *dep1,
                                        Dep *dep2, PlutoProg *prog,
                                        PlutoMatrix *access_matrix);

PlutoDepList *pluto_dep_list_alloc(Dep *dep);

void pluto_detect_scalar_dimensions(PlutoProg *prog);
void pluto_opencl_codegen();

Array *pluto_dist_get_corrs_array(char *arr_name, const PlutoProg *prog);

PlutoConstraints *pluto_dist_get_required_data_tiles(PlutoConstraints *domain,
                                                     int copylevel,
                                                     char *arr_name,
                                                     PlutoProg *prog);

char *pluto_dist_modify_stmt_text(char *stmt_text, int use_stride,
                                  PlutoProg *prog);
int pluto_dist_apply_mod_optimization(char *stmt_text, char *new_stmt_text,
                                      char *decls, const Stmt *stmt,
                                      PlutoProg *prog);
char *pluto_dist_malloc_stmt_text(char *arr_name, PlutoProg *prog,
                                  int usr_strides);

char *pluto_dist_ptr_init_stmt_text(char *arr_name, PlutoProg *prog);

char *pluto_dist_get_ptr_dim_stride(Array *arr, int curr_dim, PlutoProg *prog);
char *pluto_dist_get_ptr_dim_stride_malloc(Array *arr, int curr_dim,
                                           PlutoProg *prog);

void add_data_dist_parm_decl(FILE *fp, char *acc_name, PlutoProg *prog);
void pluto_dist_print_tiles_strided_access(char *new_stmt_text, Array *arr,
                                           PlutoProg *prog);

void pluto_dist_add_upper_lower_bounds(PlutoConstraints *cst, PlutoProg *prog);
void pluto_data_dist_identity_trans(PlutoProg *prog);
void pluto_dist_print_tiles_index_string(char *new_stmt_text, Array *arr,
                                         PlutoProg *prog);
void restore_data_hyperplane(Array *arr, int row);

int pluto_shared_memory_data_dist(PlutoProg *prog, FILE *headerfp, FILE *outfp);
int pluto_sharedmem_data_dist_codegen(PlutoProg *prog, Ploop **loops,
                                      int nloops, int *copy_level,
                                      FILE *outfp);

int pluto_dist_get_innermost_loop(const Stmt *stmt, const PlutoProg *prog);

void pluto_dist_generate_copy_back_code(Ploop **loops, int nloops,
                                        int *copy_level, FILE *headerfp,
                                        PlutoProg *prog);

void gen_copy_back_function(Stmt ****copy_back_stmts, int *num_write_data,
                            int num_copy_back_stmts, int nloops,
                            PlutoProg *prog, FILE *headerfp);

void print_data_dist_parm_call(char *str, Array *arr);
PlutoConstraints *compute_read_in(struct stmt_access_pair *racc_stmt,
                                  int copy_level, PlutoProg *prog);

void generate_pack_or_unpack(FILE *packfp, PlutoProg *prog,
                             PlutoConstraints *constraints, char *stmttext,
                             PlutoStmtType stmttype, char **iters,
                             int src_copy_level, int acc_nrows, char *acc_name,
                             char *send_count);

int pluto_dist_gen_intial_code(int nloops, Ploop **loops, int *copy_level,
                               PlutoProg *prog, FILE *headerfp);

int pluto_dist_gen_compute_pi(int nloops, Ploop **loops, int *copy_level,
                              PlutoProg *prog, int *pi_mappings, FILE *outfp,
                              FILE *headerfp);

void generate_sigma_common(struct stmt_access_pair **wacc_stmts, int naccs,
                           int *copy_level, PlutoProg *prog, int loop_num,
                           int *pi_mappings, int isRemoteOnly, int *broadcast,
                           int *num_dest_loops,
                           PlutoConstraints **tdpoly_list,
                           int *dep_loop_num_list,
                           PlutoMatrix **trans,
                           char ***iters,
                           char **indices);

int pluto_dist_read_tile_sizes(int *tile_sizes, int num_tile_dims);

int pluto_are_stmts_fused(Stmt **stmts, int nstmts, const PlutoProg *prog);

void pluto_iss_dep(PlutoProg *prog);
PlutoConstraints *pluto_find_iss(const PlutoConstraints **doms, int ndoms,
                                 int npar, PlutoConstraints *);
void pluto_iss(Stmt *stmt, PlutoConstraints **cuts, int num_cuts,
               PlutoProg *prog);

void populate_scaling_csr_matrices_for_pluto_program(int ***index,
                                                     double ***val, int nrows,
                                                     PlutoProg *prog);
PlutoMatrix *construct_cplex_objective(const PlutoConstraints *cst,
                                       const PlutoProg *prog);

#ifdef GLPK
Graph *build_fusion_conflict_graph(PlutoProg *prog, int *colour, int num_nodes,
                                   int current_colour);
void find_permutable_dimensions_scc_based(int *colour, PlutoProg *prog);
void introduce_skew(PlutoProg *prog);
#endif

#if defined(__cplusplus)
}
#endif

#endif
