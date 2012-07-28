#ifndef __LIBPLUTO__
#define __LIBPLUTO__
#include "isl/union_set.h"
#include "isl/union_map.h"

#if defined(__cplusplus)
extern "C" {
#endif
struct plutoOptions{

    /* To tile or not? */
    int tile;

    /* Intra-tile optimization */
    int intratileopt;

    /* parallelization */
    int parallel;

    /* prefer pure inner parallelism to pipelined parallelism */
    int innerpar;

    /* Automatic unroll/unroll-jamming of loops */
    int unroll;

    /* unroll/jam factor */
    int ufactor;

    /* Enable or disable post-transformations to make code amenable to
     * vectorization (default - enabled) */
    int prevector;

    /* consider RAR dependences */
    int rar;

    /* Decides the fusion algorithm (MAXIMAL_FUSE, NO_FUSE, or SMART_FUSE) */
    int fuse;

    /* for debugging - print default cloog-style total */
    int scancount;

    /* parameters will be more than at least this much */
    /* setting parameter context for cloog */
    int context;

    /* multiple (currently two) degrees of pipelined parallelism */
    int multipipe;

    /* Tile for L2 too */
    /* By default, only L1 tiling is done; under parallel execution, every
     * processor executes a sequence of L1 tiles (OpenMP adds another blocking
     * on the parallel loop). With L2 tiling, each processor executes a
     * sequence of L2 tiles and barrier is done after a group of L2 tiles is
     * exectuted -- causes load imbalance due to pipe startup when problem
     * sizes are not huge */
    int l2tile;


    /* NOTE: --ft and --lt are to manually force tiling depths */
    /* First depth to tile (starting from 0) */
    int ft;
    /* Last depth to tile (indexed from 0)  */
    int lt;

    /* Output for debugging */
    int debug;

    /* More debugging output */
    int moredebug;

    /* Not implemented yet: Don't output anything unless something fails */
    int quiet;

    /* Pure polyhedral unrolling (instead of postpass) */
    int polyunroll;

    /* Identity transformation */
    int identity;

    /* Generate scheduling pragmas for Bee+Cl@k */
    int bee;

    /* Force this for cloog's -f */
    int cloogf;

    /* Force this for cloog's -l */
    int cloogl;

    /* Use isl to compute dependences */
    int isldep;

    /* Compact dependences with ISL */
    int isldepcompact;

    /* Compute lastwriter for dependences */
    int lastwriter;

    /* DEV: Don't use cost function */
    int nobound;

    /* Ask candl to privatize */
    int scalpriv;

    /* No output from Pluto if everything goes right */
    int silent;

    /* Read input from a scoplib file */
    int readscoplib;

    /* Use isl as ilp solver. */
    int islsolve;
    /* Output file name supplied from -o */
    char *out_file;
};
typedef struct plutoOptions PlutoOptions;

PlutoOptions *pluto_options_alloc();
void pluto_options_free(PlutoOptions *);

__isl_give isl_union_map *pluto_schedule(isl_union_set *domains,
        isl_union_map *dependences,
        PlutoOptions *options);

#if defined(__cplusplus)
}
#endif
#endif
