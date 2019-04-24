/*
 * PLUTO: An automatic parallelizer and locality optimizer
 * 
 * Copyright (C) 2007-2012 Uday Bondhugula
 *
 * This file is part of Pluto.
 *
 * Pluto is free software; you can redistribute it and/or modify
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
 */
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "pluto.h"
#include "math_support.h"
#include "ddg.h"

/* Allocate a graph */
Graph *graph_alloc(int nVertices)
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

    /* Not computed yet */
    g->num_ccs = -1;
    return g;
}


/* Print the strongly-connected components */
void graph_print_sccs(Graph *g)
{
    int i;

    /* SCCs should have been computed */
    assert(g->num_sccs != -1);

    for (i=0; i<g->num_sccs; i++)  {
        IF_DEBUG(printf("\tSCC %d: size: %d: max stmt dim: %d\n",
                g->sccs[i].id, g->sccs[i].size, g->sccs[i].max_dim));
    }
}


/* Return transpose of a graph G
 * G^T has an edge b/w u and v  iff G has an edge b/w v and u */
Graph *graph_transpose(Graph *g)
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

/* Returns an undirected graph corresponding to the input directed graph */
/* This is used to find the connected components in the graph */
Graph* get_undirected_graph(const Graph *g)
{
    int i,j;
    Graph *gU;
    gU=graph_alloc(g->nVertices);
    for (i=0;i<g->nVertices;i++){
        for (j=0;j<=i;j++){
            gU->adj->val[i][j]=(g->adj->val[i][j]==0)?g->adj->val[j][i]:g->adj->val[i][j];
            gU->adj->val[j][i]=gU->adj->val[i][j];
        }
    }
    return gU;
}

/* Floyd-Warshall algorithm to find the transitive closure of a graph */
void transitive_closure (Graph *g)
{
    int i,j,k;
    for (i=0; i<g->nVertices; i++){
        for (j=0; j<g->nVertices; j++){
            for (k=0;k<g->nVertices;k++){
                g->adj->val[i][j] = ((g->adj->val[i][j]) || (g->adj->val[i][k] && g->adj->val[k][j]));
            }
        }
    }
}

/* Depth first search from a given vertex */
void dfs_vertex(Graph *g, Vertex *v, int *time)
{
    int j;

    *time = *time + 1;
    v->vn = *time;

    /* matrix_print(stdout, g->adj, g->nVertices, g->nVertices); */

    for (j=0; j< g->nVertices; j++) {
        if (g->adj->val[v->id][j])   {
            // Vertex *w = ddg_get_vertex_by_id(g, j);
            g->vertices[j].cc_id = v->cc_id;
            if (g->vertices[j].vn == 0) {
                dfs_vertex(g, &g->vertices[j], time);
            }
        }
    }

    *time = *time + 1;
    v->fn = *time;
}


/* Depth first search */
void dfs(Graph *g)
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

/*
 * Depth first search in the descending order of finish
 * times of previous dfs - used to compute SCCs
 *
 * SCC related information is stored
 *
 **/
void dfs_for_scc(Graph *g)
{
    int i, j;
    int time = 0;

    Vertex *vCopy = (Vertex *) malloc(g->nVertices*sizeof(Vertex));
    memcpy(vCopy, g->vertices, g->nVertices*sizeof(Vertex));
    
    for (i=0; i<g->nVertices; i++)  {
        g->vertices[i].vn = 0;
        g->vertices[i].scc_id = -1;
    }

    /* Sort (in the ascending order) by finish number */
    qsort(vCopy, g->nVertices, sizeof(Vertex), compar);

    /* Reverse order of fn */
    int numScc = 0;
    for (i=g->nVertices-1; i>=0; i--)  {
        if (g->vertices[vCopy[i].id].vn == 0)  {
            g->sccs[numScc].id = numScc;
            dfs_vertex(g, &g->vertices[vCopy[i].id], &time);

            IF_MORE_DEBUG(printf("[pluto] dfs_for_scc: SCC %d: Stmt ids: ", numScc));
            for (j=0; j<g->nVertices; j++) {
                if (g->vertices[j].scc_id == -1 && g->vertices[j].vn > 0) {
                    g->vertices[j].scc_id = numScc;
                    IF_MORE_DEBUG(printf(" %d", g->vertices[j].id));
                }
            }
            IF_MORE_DEBUG(printf("\n"));
            numScc++;
        }
    }
    IF_MORE_DEBUG(printf("\n\n"));

    g->num_sccs = numScc;

    free(vCopy);
}

bool is_adjecent(Graph *g, int i, int j){
    /* PlutoMatrix *adj; */
    /* adj = g->adj; */
    if(g->adj->val[i][j] != 0|| g->adj->val[j][i] != 0){
        return true;
    }
    return false;
}

void compute_scc_vertices(Graph *ddg){
    int i,j,k;
    int n_sccs;
    int *vertices;

    n_sccs = ddg->num_sccs;
    for(i=0; i<n_sccs; i++){
        vertices = (int*)malloc((ddg->sccs[i].size)*sizeof(int));
        k = 0;
        for(j=0;j<ddg->nVertices; j++){
            if((ddg->vertices[j].scc_id) == i){
                vertices[k] = ddg->vertices[j].id;
                k++;
            }
        }
        ddg->sccs[i].vertices = vertices;
    }
}

void print_scc_vertices(int scc_id, Graph *g){
    int i;
    for (i=0; i<g->sccs[scc_id].size; i++){
        printf("S%d, ",g->sccs[scc_id].vertices[i]);
    }
    printf("\n");
}

void free_scc_vertices(Graph *ddg){
    int i;
    int *vertices;
    for(i=0;i<ddg->num_sccs;i++){
        vertices = ddg->sccs[i].vertices;
        if(vertices!=NULL){
            free(vertices);
        }
    }
}
void graph_free(Graph *g)
{
    pluto_matrix_free(g->adj);
    free(g->vertices);
    free(g->sccs);
    free(g);
}

