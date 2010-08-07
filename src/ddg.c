#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "pluto.h"
#include "math_support.h"
#include "ddg.h"

/* Allocate a graph */
Graph *graph_alloc (int nVertices)
{
    Graph *g;
    int i, j;

    g = (Graph *) malloc(sizeof(Graph));

    g->nVertices = nVertices;

    g->vertices = (Vertex *) malloc(nVertices*sizeof(Vertex));
    for (i=0; i<nVertices; i++) {
        g->vertices[i].id = i;
    }

    g->adj = pluto_matrix_alloc(nVertices, nVertices);
    g->adj->nrows = nVertices;
    g->adj->ncols = nVertices;

    for (i=0; i<nVertices; i++) {
        g->vertices[i].vn = -1;
        g->vertices[i].fn = -1;
        for (j=0; j<nVertices; j++)
            g->adj->val[i][j] = 0;
    }

    g->sccs = (Scc *) malloc(nVertices*sizeof(Scc));

    /* Not computed yet */
    g->num_sccs = -1;

    return g;
}


/* Print the strongly-connected components */
void graph_print_sccs (Graph *g)
{
    int i;

    /* SCCs should have been computed */
    assert(g->num_sccs != -1);

    for (i=0; i<g->num_sccs; i++)  {
        IF_DEBUG(printf("SCC_id %d: size: %d: max stmt dim: %d\n", 
                g->sccs[i].id, g->sccs[i].size, g->sccs[i].max_dim));
    }
}


/* Return transpose of a graph G
 * G^T has an edge b/w u and v  iff G has an edge b/w v and u */
Graph *graph_transpose (Graph *g)
{
    int i, j;
    Graph *gT;

    gT = graph_alloc(g->nVertices);

    for (i=0; i<g->nVertices; i++)  {
        for (j=0; j<g->nVertices; j++)  {
            gT->adj->val[j][i] = g->adj->val[i][j];
        }

        gT->vertices[i].fn = g->vertices[i].fn;
    }

    return gT;
}


/* Depth first search from a given vertex */
void dfs_vertex (Graph *g, Vertex *v, int *time)
{
    int j;

    *time = *time + 1;
    v->vn = *time;

    /* matrix_print(stdout, g->adj, g->nVertices, g->nVertices); */

    for (j=0; j< g->nVertices; j++) {
        if (g->adj->val[v->id][j])   {
            // Vertex *w = ddg_get_vertex_by_id(g, j);
            if (g->vertices[j].vn == 0) {
                dfs_vertex(g, &g->vertices[j], time);
            }
        }
    }

    *time = *time + 1;
    v->fn = *time;
}


/* Depth first search */
void dfs (Graph *g)
{
    int i;
    int time = 0;
    
    for (i=0; i<g->nVertices; i++)  {
        g->vertices[i].vn = 0;
    }
    
    for (i=0; i<g->nVertices; i++)  {
        if (g->vertices[i].vn == 0)  {
            // printf("DFS vertex: %d\n", i);
            dfs_vertex (g, &g->vertices[i], &time);
        }
    }
}


/* Comparison function for sorting graph vertices by their finish time */
static int compar (const void *e1, const void *e2)
{
    Vertex *v1, *v2;

    v1 = (Vertex *) e1;
    v2 = (Vertex *) e2;

    if (v1->fn < v2->fn)    {
        return -1;
    }else if (v1->fn == v2->fn) {
        return 0;
    }else return 1;
}

/* Depth first search - this version stores additional data
 * that is useful for computing SCCs */
void dfs_for_scc (Graph *g)
{
    int i, j;
    int time = 0;

    Vertex *vCopy = (Vertex *) malloc(g->nVertices*sizeof(Vertex));
    memcpy(vCopy, g->vertices, g->nVertices*sizeof(Vertex));
    
    for (i=0; i<g->nVertices; i++)  {
        g->vertices[i].vn = 0;
        g->vertices[i].scc_id = -1;
    }

    /* Sort by fn */
    qsort(vCopy, g->nVertices, sizeof(Vertex), compar);

    /* Reverse order of fn */
    int numScc = 0;
    for (i=g->nVertices-1; i>=0; i--)  {
        if (g->vertices[vCopy[i].id].vn == 0)  {
            g->sccs[numScc].id = numScc;
            dfs_vertex (g, &g->vertices[vCopy[i].id], &time);

            IF_DEBUG(printf("SCC %d: Stmts: ", numScc));
            for (j=0; j<g->nVertices; j++) {
                if (g->vertices[j].scc_id == -1 && g->vertices[j].vn > 0) {
                    g->vertices[j].scc_id = numScc;
                    IF_DEBUG(printf(" %d", g->vertices[j].id));
                }
            }
            IF_DEBUG(printf("\n"));
            numScc++;
        }
    }
    IF_DEBUG(printf("\n\n"));

    g->num_sccs = numScc;

    free(vCopy);
}

void graph_free(Graph *g)
{
    pluto_matrix_free(g->adj);
    free(g->vertices);
    free(g->sccs);
    free(g);
}

