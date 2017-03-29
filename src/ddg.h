#ifndef _DDG_H_
#define _DDG_H_

#include "math_support.h"

/* Vertex of a graph */
struct vertex{
    /* In PLUTO code, id here is same as Stmt ID */
    int id;
    int vn;
    int fn;

    /* Id of the SCC this vertex belongs to */
    int scc_id;
};
typedef struct vertex Vertex;

struct scc{
    /* Number of vertices in it */
    int size;

    /* Maximum dimensionality statement in this SCC */
    int max_dim;

    /* Id of this SCC */
    int id;

    /* Stmt Id's of the vertices in the SCC */
    int* vertices;
};
typedef struct scc Scc;


struct graph{
    /* List of vertices. For Pluto, this list vertices directly corresponds to
     * the list of statements in the order they appear in Stmt *stmts with a
     * vertex ID being same as the statement ID */
    Vertex *vertices;
    int nVertices;

    /* Adjacency matrix */
    PlutoMatrix *adj;

    Scc *sccs;
    int num_sccs;
};
typedef struct graph Graph;

Graph *graph_alloc (int nVertices);
void graph_free(Graph *g);
void graph_print_sccs (Graph *g);
void dfs_for_scc (Graph *g);
Graph *graph_transpose (Graph *g);
void dfs (Graph *g);
void dfs_for_scc (Graph *g);
Vertex *ddg_get_vertex_by_id(Graph *g, int id);

void transitive_closure(Graph *graph);
void compute_scc_vertices(Graph *ddg);

#endif
