#ifndef _DDG_H_
#define _DDG_H_

#include <stdbool.h>

typedef struct pluto_matrix PlutoMatrix;
typedef struct plutoContext PlutoContext;

/* Vertex of a graph */
struct vertex {
  /* In PLUTO code, id here is same as Stmt ID */
  int id;
  int vn;
  int fn;

  /* This gives the offset of the Statement in the FCG. In the
   * FCG this field will give the id of the statement */
  int fcg_stmt_offset;
  /* Id of the SCC this vertex belongs to */
  int scc_id;
  /* Id of the Connected Component this vertex belongs to */
  int cc_id;
};
typedef struct vertex Vertex;

struct scc {
  /* Number of vertices in it */
  int size;

  /* Maximum dimensionality statement in this SCC */
  int max_dim;

  /* Id of this SCC */
  int id;

  /* Stmt Id's of the vertices in the SCC */
  int *vertices;

  /* Set to 1 if the scc is parallel. */
  int is_parallel;

  /* Rational hyperplane that weakly satisfies all the dependences in the SCC.
   * This can be scaled to integeral hyperplane */
  double *sol;

  /* Points to the first vertex in the FCG corresponding to this scc,
   * when used with SCC based clustering heuristic */
  int fcg_scc_offset;

  /* Set to true if the scc is coloured with current colour else false */
  bool is_scc_coloured;

  /* Set to true if there is a parallel hyperplane has already been found for
   * this scc */
  bool has_parallel_hyperplane;
  /* Set to true if the scc has stencil dependence pattern. */
  bool is_scc_stencil;
};
typedef struct scc Scc;

struct graph {
  /* List of vertices. For Pluto, this list vertices directly corresponds to
   * the list of statements in the order they appear in Stmt *stmts with a
   * vertex ID being same as the statement ID */
  Vertex *vertices;
  int nVertices;

  /* Number of vertices that have already been coloured */
  int num_coloured_vertices;

  /* Adjacency matrix */
  PlutoMatrix *adj;

  Scc *sccs;
  int num_sccs;

  int num_ccs;

  /* Indicates whether the graph has to be rebuilt. Used by pluto-dfp in order
   * to
   * rebuild only when necessary */
  bool to_be_rebuilt;

  PlutoContext *context;
};
typedef struct graph Graph;

#if defined(__cplusplus)
extern "C" {
#endif

Graph *graph_alloc(int nVertices, PlutoContext *context);
void graph_free(Graph *g);
void graph_print_sccs(Graph *g);
void dfs_for_scc(Graph *g);
Graph *get_undirected_graph(const Graph *g);
Graph *graph_transpose(Graph *g);
void dfs(Graph *g);
void dfs_for_scc(Graph *g);
void dfs_vertex(Graph *g, Vertex *v, int *time);
void transitive_closure(Graph *g);
Vertex *ddg_get_vertex_by_id(Graph *g, int id);

bool is_adjecent(Graph *, int, int);
int *get_ssc_topological_order(Graph *ddg);
void compute_scc_vertices(Graph *ddg);
void print_scc_vertices(int scc_id, Graph *g);
void free_scc_vertices(Graph *ddg);

#if defined(__cplusplus)
}
#endif

#endif
