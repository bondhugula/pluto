#ifdef __USE_DAG_SCHEDULER_AS_LIB__
#include "scheduler.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <math.h>

#include <string.h>
#include <stdint.h>

#include <omp.h>

#define WHITE 1
#define GREY 2
#define BLACK 3

void calculate_tlevel_blevel(void);
void dfs_visit(int vertexId, int *nextIndex, int *reverseTopSort);

int dag_vertices_not_complete();
void execute_task(int taskId, int threadId);

#define DATA_ACCESS_BW 1 

uint32_t hashword(const uint32_t *k, size_t length, uint32_t initVal);
int get_hash_index(int i, int j, int k, int l);

#define MAX_VALUE 1000000
#define MIN_VALUE -1000000
#define HTABLE_SIZE 256

typedef struct taskDescr /* The tuple describing the task tile */
{
	int i;
	int j;
	int k;
	int l;
}TaskDescr;

typedef struct edge
{
	int svertexId; /* Source vertex id */
	int dvertexId; /* Destination vertex id */
	double edgeWt; /* Edge wt is a measure of reuse volume between the incident vertices (tasks) */
	struct edge *pNext;	
}Edge; 

typedef struct vertex
{
	int vertexId;
	TaskDescr taskDescr;
	double vertexWt; /* Vertex wt is a measure of the expected execution time of the corresponding task */
	int nParents;
	Edge *parentList; /* Populated with the corresponding edge ids */
	int nChildren;
	Edge *childList;
	int procId;
	int DFSColor; /* Denotes color of vertex (white, grey, black) during DFS */
	double tLevel; /* top level */
	double bLevel; /* bottom level */
	double est; /* Earliest start time */
	double ast; /* Actual start time */
	double complTime; /* Completion Time (ast + data access time + exec time*/
	int complFlag;
}Vertex;

typedef struct dag
{
	int nVertices;
	int nEdges;
	Vertex *vertexList;
}DAG;

#ifdef __STATIC_SCHEDULE__

typedef struct resourceFreeList /* Data structure for insertion scheduling */
{
	int nFree;
	int *freeSet;
	double startTime;
	double endTime;
	struct resourceFreeList *pNext;
}ResourceFreeList;

#endif

typedef struct resourceInfo
{
	int nProcs;
#ifdef __STATIC_SCHEDULE__
	ResourceFreeList *pFreeList;
#endif
}ResourceInfo;

#ifdef __STATIC_SCHEDULE__
typedef struct schedInfo
{
    int **taskQueue;
    int *nTasks;
}ScheduleInfo;
#endif

#ifndef __STATIC_SCHEDULE__
typedef struct pqueue
{ 
    int size;
    int avail;
    int *pqData;
}PQueue;
#endif

typedef struct hashnode
{
   int vertexId;
   int i;
   int j;
   int k;
   int l;
   struct hashnode *pNext;
}HashNode; 

#ifndef __STATIC_SCHEDULE__

PQueue *pq_init(int nElems);
void pq_insert(PQueue *pq, int n);
int pq_remove(PQueue *pq);
void pq_delete(PQueue *pq);
double get_priority(int n);
#endif

#ifdef __STATIC_SCHEDULE__
void schedule_task(int vertexId);
void populate_task_queue(void);
#endif

#ifdef __STATIC_SCHEDULE__

int is_present(int item, int *arr, int nElems);
void create_reservation(double ast, double dur, int proc);
double find_min_starttime_for_task_on_proc(int proc, double est, double dur);

#endif

#define DEC_BLEVEL 1
#define INC_ALAP 2
#define INC_AST 3

void sort_priorities(int *sortList, int nEntries, int sortFlag);
void merge_sort(int *initList, int *sortList, int left, int right, int sortFlag);
void merge(int *list1, int *list2, int left, int mid, int right, int sortFlag);

/* C files */

DAG dag;
int nVertices = 0;
HashNode* hashTable[HTABLE_SIZE];

#ifdef __STATIC_SCHEDULE__
ResourceInfo resource;
ScheduleInfo schedule;
#endif

void calculate_tlevel_blevel(void)
{
    int *reverseTopSort, *topSort, i, *nextIndex, j;
    double maxValue = 0;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering calculate_tlevel_blevel()\n");
#endif

    reverseTopSort = (int *)malloc(sizeof(int)*dag.nVertices);    
    if (reverseTopSort == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    topSort = (int *)malloc(sizeof(int)*dag.nVertices);    
    if (topSort == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    nextIndex = (int *)malloc(sizeof(int));
    if (nextIndex == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    for (i=0; i<dag.nVertices; i++)
    {
        reverseTopSort[i] = -1;
        topSort[i] = -1;
        dag.vertexList[i].DFSColor = WHITE;    
        dag.vertexList[i].bLevel = 0;
        dag.vertexList[i].tLevel = 0;
    }
    *nextIndex = 0;
    /* DFS to sort vertices in reverse topological order */
    for (i=0; i<dag.nVertices; i++)
    {
        if (dag.vertexList[i].DFSColor == WHITE)
        {
            dfs_visit(i, nextIndex, reverseTopSort);
        }
    }
    /* Calculate top level */
    j = dag.nVertices - 1;
    for (i=0; i<dag.nVertices; i++)
    {
        topSort[j] = reverseTopSort[i];
        j--;
    }
    i = 0;
    while (i < dag.nVertices)
    {
        int vertexId;
        Edge *pEdge;

        vertexId = topSort[i];
        maxValue = 0;
        pEdge = dag.vertexList[vertexId].parentList;
        while (pEdge)
        {
            int svertexId;

            svertexId = pEdge->svertexId;
            if (dag.vertexList[svertexId].tLevel + dag.vertexList[svertexId].vertexWt + pEdge->edgeWt > maxValue)
            {
                maxValue = dag.vertexList[svertexId].tLevel + dag.vertexList[svertexId].vertexWt + pEdge->edgeWt;
            }
            pEdge = pEdge->pNext;
        }
        dag.vertexList[vertexId].tLevel = maxValue;
        i++;
    }

    /* Calculate bottom level */
    i = 0;
    while (i < dag.nVertices)
    {
        int vertexId;
        Edge *cEdge;

        vertexId = reverseTopSort[i];
        maxValue = 0;
        cEdge = dag.vertexList[vertexId].childList;
        while (cEdge)
        {
            int dvertexId;

            dvertexId = cEdge->dvertexId;
            if (cEdge->edgeWt + dag.vertexList[dvertexId].bLevel > maxValue)
            {
                maxValue = cEdge->edgeWt + dag.vertexList[dvertexId].bLevel;
            }
            cEdge = cEdge->pNext;
        }
        dag.vertexList[vertexId].bLevel = maxValue + dag.vertexList[vertexId].vertexWt;
        i++;
    }
    free(reverseTopSort);
    free(topSort);
    free(nextIndex);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting calculate_tlevel_blevel()\n");
#endif
    return;
}

void dfs_visit(int vertexId, int *nextIndex, int *reverseTopSort)
{
    Edge *cEdge;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering dfs_visit() for vertexId %d\n", vertexId);
#endif

    dag.vertexList[vertexId].DFSColor = GREY;

    cEdge = dag.vertexList[vertexId].childList;
    while (cEdge)
    {
        if (dag.vertexList[cEdge->dvertexId].DFSColor == WHITE)
        {
            dfs_visit(cEdge->dvertexId, nextIndex, reverseTopSort);
        }
        cEdge = cEdge->pNext;
    }
    dag.vertexList[vertexId].DFSColor = BLACK;
    reverseTopSort[*nextIndex] = vertexId;
    *nextIndex = *nextIndex + 1;
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting dfs_visit() for vertexId %d\n", vertexId);
#endif
    return;
}

int dag_vertices_not_complete()
{
	int i;
	for (i = 0; i < dag.nVertices; ++i)
		if (dag.vertexList[i].complFlag != 1)
			return 1;
	return 0;
}

void dag_execute()
{
    int i;
#ifndef __STATIC_SCHEDULE__
    PQueue *taskQ;
    omp_lock_t taskQLock;
#else
    int index, *syncFlag;
#endif

#ifdef __DEBUG2__
    fprintf(stderr, "Entering dag_execute()\n");
#endif

#ifdef __STATIC_SCHEDULE__
    index = 0;
    syncFlag = (int*)malloc(sizeof(int)*dag.nVertices);
    if (syncFlag == 0)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    for (i=0; i<dag.nVertices; i++)
    {
        syncFlag[i] = 0;
    }
#endif
#ifdef __STATIC_SCHEDULE__
#pragma omp parallel private(index) shared(syncFlag)
    {
        int threadId = omp_get_thread_num();
        for (index = 0; index < schedule.nTasks[threadId]; index++)
        {
            int taskId, i, syncOut;
            Edge *pEdge;

            taskId = schedule.taskQueue[threadId][index];
            if (dag.vertexList[taskId].vertexWt > 0) /* Zero vertex weight implies a dummy source or sink vertex */
            {
                while (syncFlag[taskId] == 0); /* How to prevent starvation? */
#ifdef __DEBUG1__
    fprintf(stderr, "Executing task %d (<%d, %d, %d, %d>) by proc %d\n", taskId, dag.vertexList[taskId].taskDescr.i, dag.vertexList[taskId].taskDescr.j, dag.vertexList[taskId].taskDescr.k, dag.vertexList[taskId].taskDescr.l, threadId);
#endif
/* To assert if predecessors have finished */
                pEdge = dag.vertexList[taskId].parentList;
                for (i=0; i<dag.vertexList[taskId].nParents; i++)
                {    
                    assert(pEdge != NULL);
                    assert(dag.vertexList[pEdge->svertexId].complFlag == 1);
                    pEdge = pEdge->pNext;
                }        
                task_computation(dag.vertexList[taskId].taskDescr.i, dag.vertexList[taskId].taskDescr.j, dag.vertexList[taskId].taskDescr.k, dag.vertexList[taskId].taskDescr.l);

            }
            dag.vertexList[taskId].complFlag = 1;
            pEdge = dag.vertexList[taskId].childList;
            for (i=0; i<dag.vertexList[taskId].nChildren; i++)
            {
                int childId, j;
                Edge *pPreds;

                assert(pEdge != NULL);
                childId = pEdge->dvertexId;
                pPreds = dag.vertexList[childId].parentList;
                syncOut = 1;
                for (j=0; j<dag.vertexList[childId].nParents; j++)
                {
                    assert(pPreds != NULL);
                    if (dag.vertexList[pPreds->svertexId].complFlag == 0)
                    {
                        syncOut = 0;
                        break;
                    }
                    pPreds = pPreds->pNext; 
                }
                if (syncOut == 1)
                {
                #pragma omp atomic
                    syncFlag[childId] += 1;
                }
                pEdge = pEdge->pNext;
            }
        }
    }
    free(syncFlag);
#else
    calculate_tlevel_blevel();
    taskQ = pq_init(dag.nVertices);    
    for (i=0; i<dag.nVertices; i++)
    {
        if (dag.vertexList[i].nParents == 0)
        {
            pq_insert(taskQ, i);
        }
    }
    omp_init_lock(&taskQLock);
#pragma omp parallel shared(taskQLock, taskQ)
    {
        int taskId, i, syncOut;
        Edge *pEdge;

#ifdef __DEBUG1__
        int threadId = omp_get_thread_num();
#endif
        while (dag_vertices_not_complete())
        {
            omp_set_lock(&taskQLock);
            taskId = pq_remove(taskQ);
            omp_unset_lock(&taskQLock);
            if (taskId == -1) continue;
            if (dag.vertexList[taskId].vertexWt > 0) /* Zero vertex weight implies dummy source or sink vertex */
            {
#ifdef __DEBUG1__
                fprintf(stderr, "Executing task %d (<%d, %d, %d, %d>) by proc %d\n", taskId, dag.vertexList[taskId].taskDescr.i, dag.vertexList[taskId].taskDescr.j, dag.vertexList[taskId].taskDescr.k, dag.vertexList[taskId].taskDescr.l, threadId);
#endif
                /* To assert if predecessors have finished */
                pEdge = dag.vertexList[taskId].parentList;
                for (i=0; i<dag.vertexList[taskId].nParents; i++)
                {    
                    assert(pEdge != NULL);
                    assert(dag.vertexList[pEdge->svertexId].complFlag == 1);
                    pEdge = pEdge->pNext;
                }        
                task_computation(dag.vertexList[taskId].taskDescr.i, dag.vertexList[taskId].taskDescr.j, dag.vertexList[taskId].taskDescr.k, dag.vertexList[taskId].taskDescr.l);
            }
            dag.vertexList[taskId].complFlag = 1;
            pEdge = dag.vertexList[taskId].childList;
            for (i=0; i<dag.vertexList[taskId].nChildren; i++)
            {
                int childId, j;
                Edge *pPreds;

                assert(pEdge != NULL);
                childId = pEdge->dvertexId;
                pPreds = dag.vertexList[childId].parentList;
                syncOut = 1;
                for (j=0; j<dag.vertexList[childId].nParents; j++)
                {
                    assert(pPreds != NULL);
                    if (dag.vertexList[pPreds->svertexId].complFlag == 0)
                    {
                        syncOut = 0;
                        break;
                    }
                    pPreds = pPreds->pNext; 
                }
                if (syncOut == 1)
                {
                    omp_set_lock(&taskQLock);
                    pq_insert(taskQ, childId);
                    omp_unset_lock(&taskQLock);
                }
                pEdge = pEdge->pNext;
            }
        }
    }
    omp_destroy_lock(&taskQLock);
    pq_delete(taskQ);
#endif 

#ifdef __DEBUG2__
    fprintf(stderr, "Exiting dag_execute()\n");
#endif
}

void init_dag(int nVertices)
{
    int i;
#ifdef __DEBUG2__
    fprintf(stderr, "Entering init_dag() with %d vertices\n", nVertices);
#endif
    dag.nVertices = nVertices;
    dag.nEdges = 0;
    dag.vertexList = (Vertex *)malloc(sizeof(Vertex)*nVertices);
    if (dag.vertexList == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    for (i=0; i<HTABLE_SIZE; i++)
    {
        hashTable[i] = NULL;
    }
    
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting init_dag() with %d vertices\n", nVertices);
#endif
    return;
}

void dag_add_vertex(int i, int j, int k, int l, double vertexWt, int *vertexId)
{
    Vertex *vertex;
    HashNode *pHashNode;
    int hashIndex;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering dag_add_vertex() with tuple <%d, %d, %d, %d> and wt %.2lf to DAG\n", i, j, k, l, vertexWt);
#endif
    hashIndex = get_hash_index(i, j, k, l);
    pHashNode = hashTable[hashIndex];
    while (pHashNode != NULL)
    {
        if ((pHashNode->i == i) && (pHashNode->j == j) && (pHashNode->k == k) && (pHashNode->l == l))
        {
#ifdef __DEBUG2__
            fprintf(stderr, "Duplicate vertex <%d, %d, %d, %d> \n", i, j, k, l);
#endif
            return;
        }
        pHashNode = pHashNode->pNext;
    } 
        
    *vertexId = nVertices;
    if (nVertices > (dag.nVertices - 1))
    {
#ifdef __DEBUG2__
        fprintf(stderr, "Wrong estimation of %d vertices in DAG\n", dag.nVertices);
#endif
        dag.vertexList = (Vertex *)realloc(dag.vertexList, nVertices + 1);
        if (dag.vertexList == NULL)
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
            exit(-1);
        }
        dag.nVertices = nVertices + 1;
    }
    vertex = &dag.vertexList[*vertexId];
    vertex->vertexId = *vertexId;
    vertex->taskDescr.i = i;
    vertex->taskDescr.j = j;
    vertex->taskDescr.k = k;
    vertex->taskDescr.l = l;
    vertex->vertexWt = vertexWt;
    vertex->nParents = 0;
    vertex->parentList = NULL;
    vertex->nChildren = 0;
    vertex->childList = NULL;
    vertex->tLevel = 0;
    vertex->bLevel = 0;
    vertex->est = 0;
    vertex->ast = 0;
    vertex->complTime = 0;
    vertex->complFlag = 0;
    nVertices++;

    pHashNode = (HashNode *)malloc(sizeof(HashNode));
    if (pHashNode == NULL)
    { 
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    pHashNode->vertexId = *vertexId;
    pHashNode->i = i;
    pHashNode->j = j;
    pHashNode->k = k;
    pHashNode->l = l;
    pHashNode->pNext = hashTable[hashIndex];
    hashTable[hashIndex] = pHashNode;
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting dag_add_vertex() with tuple <%d, %d, %d, %d> and wt %.2lf to DAG\n", i, j, k, l, vertexWt);
#endif
    return;
} 

void dag_add_edge(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2, double reuseVolume)
{
    Edge *edge;
    int i, svertexId, dvertexId, hashIndex;
    HashNode *pHashNode;
#ifdef __DEBUG2__
    fprintf(stderr, "Entering dag_add_edge(), adding edge from vertex <%d, %d, %d, %d> to vertex <%d, %d, %d, %d> with reuse volume %.2lf to DAG\n", i1, j1, k1, l1, i2, j2, k2, l2, reuseVolume);
#endif
    
    if ((i1 == i2) && (j1 == j2) && (k1 == k2) && (l1 == l2))
    {
#ifdef __DEBUG2__
        fprintf(stderr, "Self edge \n");
#endif
        return;
    }

    hashIndex = get_hash_index(i1, j1, k1, l1);
    pHashNode = hashTable[hashIndex];
    /*if (pHashNode == NULL) return;	*/
    assert(pHashNode != NULL);
    svertexId = -1;
    while (pHashNode)
    {
        if ((pHashNode->i == i1) && (pHashNode->j == j1) && (pHashNode->k == k1) && (pHashNode->l == l1))
        {
            svertexId = pHashNode->vertexId;
            break;
        }
        pHashNode = pHashNode->pNext;
    }
    /*if (svertexId == -1) return;*/
    assert(svertexId != -1);

    hashIndex = get_hash_index(i2, j2, k2, l2);
    pHashNode = hashTable[hashIndex];
    assert(pHashNode != NULL);
    dvertexId = -1;
    while (pHashNode)
    {
        if ((pHashNode->i == i2) && (pHashNode->j == j2) && (pHashNode->k == k2) && (pHashNode->l == l2))
        {
            dvertexId = pHashNode->vertexId;
            break;
        }
        pHashNode = pHashNode->pNext;
    }
    assert(dvertexId != -1);
    edge = dag.vertexList[svertexId].childList;
    for (i=0; i<dag.vertexList[svertexId].nChildren; i++)
    {
        assert(edge != NULL);
        if ((edge->svertexId == svertexId) && (edge->dvertexId == dvertexId))
        {
#ifdef __DEBUG2__
            fprintf(stderr, "Duplicate Edge \n");
#endif
            return;
        }
        edge = edge->pNext;
    } 

    edge = (Edge *)malloc(sizeof(Edge));
    if (edge == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    edge->svertexId = svertexId;
    edge->dvertexId = dvertexId;
    edge->edgeWt = reuseVolume/DATA_ACCESS_BW;
    edge->pNext = dag.vertexList[svertexId].childList;
    dag.vertexList[svertexId].childList = edge;
    dag.vertexList[svertexId].nChildren++;
    
    edge = (Edge *)malloc(sizeof(Edge));
    if (edge == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    edge->svertexId = svertexId;
    edge->dvertexId = dvertexId;
    edge->edgeWt = reuseVolume/DATA_ACCESS_BW;
    edge->pNext = dag.vertexList[dvertexId].parentList;
    dag.vertexList[dvertexId].parentList = edge;
    dag.vertexList[dvertexId].nParents++;
    dag.nEdges++;
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting dag_add_edge(), adding edge from vertex <%d, %d, %d, %d> to vertex <%d, %d, %d, %d> with reuse volume %.2lf to DAG\n", i1, j1, k1, l1, i2, j2, k2, l2, reuseVolume);
#endif
    return;
}

void update_number_vertices(void)
{
    dag.nVertices = nVertices;
    return;
}

void free_dag(void)
{
    int i;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering free_dag()\n");
#endif 
    for (i = 0; i < dag.nVertices; i++)
    {
        Vertex *vertex;

        vertex = &dag.vertexList[i];
        while (vertex->parentList)
        {
            Edge *tmpPointer;

            tmpPointer = vertex->parentList;
            vertex->parentList = vertex->parentList->pNext;
            free(tmpPointer);
        }
        while (vertex->childList)
        {
            Edge *tmpPointer;

            tmpPointer = vertex->childList;
            vertex->childList = vertex->childList->pNext;
            free(tmpPointer);
        }
    }
    free(dag.vertexList);
    dag.nVertices = 0;
    dag.nEdges = 0;
    for (i=0; i<HTABLE_SIZE; i++)
    {
        HashNode *pHashNode;

        pHashNode = hashTable[i];
        while (pHashNode)
        {
            HashNode *pTmpHashNode;
            
            pTmpHashNode = pHashNode->pNext;
            free(pHashNode);
            pHashNode = pTmpHashNode;
        }
    } 
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting free_dag()\n");
#endif 
    return;
}

void print_dag(void)
{
    int i, j;
    Edge *pEdge;

	fprintf(stderr, "Printing DAG details ....\n");
	fprintf(stderr, "Vertex <Vertexid> <Task Tuple> <Number of Children> <Child VertexId> <Child Tuple>\n");
    for (i=0; i<dag.nVertices; i++)
    {
        fprintf(stderr, "Vertex %d <%d, %d, %d, %d> %d ", i, dag.vertexList[i].taskDescr.i, dag.vertexList[i].taskDescr.j, dag.vertexList[i].taskDescr.k, dag.vertexList[i].taskDescr.l, dag.vertexList[i].nChildren);
        pEdge = dag.vertexList[i].childList;
        for (j = 0; j <dag.vertexList[i].nChildren; j++)
        {
            assert(pEdge != NULL);
            fprintf(stderr, "%d <%d, %d, %d, %d>", pEdge->dvertexId, dag.vertexList[pEdge->dvertexId].taskDescr.i, dag.vertexList[pEdge->dvertexId].taskDescr.j, dag.vertexList[pEdge->dvertexId].taskDescr.k, dag.vertexList[pEdge->dvertexId].taskDescr.l);
            pEdge = pEdge->pNext;
        }
        fprintf(stderr, "\n");
    }    
    return;
}

/* Robert Jenkins hashing function. http://burtleburtle.net/bob/c/lookup3.c
*/

#define hashsize(n) ((uint32_t)1<<(n))
#define hashmask(n) (hashsize(n)-1)
#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))
/*
-------------------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.

This is reversible, so any information in (a,b,c) before mix() is
still in (a,b,c) after mix().

If four pairs of (a,b,c) inputs are run through mix(), or through
mix() in reverse, there are at least 32 bits of the output that
are sometimes the same for one pair and different for another pair.
This was tested for:
* pairs that differed by one bit, by two bits, in any combination
  of top bits of (a,b,c), or in any combination of bottom bits of
  (a,b,c).
* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
  is commonly produced by subtraction) look like a single 1-bit
  difference.
* the base values were pseudorandom, all zero but one bit set, or 
  all zero plus a counter that starts at zero.

Some k values for my "a-=c; a^=rot(c,k); c+=b;" arrangement that
satisfy this are
    4  6  8 16 19  4
    9 15  3 18 27 15
   14  9  3  7 17  3
Well, "9 15 3 18 27 15" didn't quite get 32 bits diffing
for "differ" defined as + with a one-bit base and a two-bit delta.  I
used http://burtleburtle.net/bob/hash/avalanche.html to choose 
the operations, constants, and arrangements of the variables.

This does not achieve avalanche.  There are input bits of (a,b,c)
that fail to affect some output bits of (a,b,c), especially of a.  The
most thoroughly mixed value is c, but it doesn't really even achieve
avalanche in c.

This allows some parallelism.  Read-after-writes are good at doubling
the number of bits affected, so the goal of mixing pulls in the opposite
direction as the goal of parallelism.  I did what I could.  Rotates
seem to cost as much as shifts on every machine I could lay my hands
on, and rotates are much kinder to the top and bottom bits, so I used
rotates.
-------------------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= c;  a ^= rot(c, 4);  c += b; \
  b -= a;  b ^= rot(a, 6);  a += c; \
  c -= b;  c ^= rot(b, 8);  b += a; \
  a -= c;  a ^= rot(c,16);  c += b; \
  b -= a;  b ^= rot(a,19);  a += c; \
  c -= b;  c ^= rot(b, 4);  b += a; \
}

/*
-------------------------------------------------------------------------------
final -- final mixing of 3 32-bit values (a,b,c) into c

Pairs of (a,b,c) values differing in only a few bits will usually
produce values of c that look totally different.  This was tested for
* pairs that differed by one bit, by two bits, in any combination
  of top bits of (a,b,c), or in any combination of bottom bits of
  (a,b,c).
* "differ" is defined as +, -, ^, or ~^.  For + and -, I transformed
  the output delta to a Gray code (a^(a>>1)) so a string of 1's (as
  is commonly produced by subtraction) look like a single 1-bit
  difference.
* the base values were pseudorandom, all zero but one bit set, or 
  all zero plus a counter that starts at zero.

These constants passed:
 14 11 25 16 4 14 24
 12 14 25 16 4 14 24
and these came close:
  4  8 15 26 3 22 24
 10  8 15 26 3 22 24
 11  8 15 26 3 22 24
-------------------------------------------------------------------------------
*/
#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}

/*
--------------------------------------------------------------------
 This works on all machines.  To be useful, it requires
 -- that the key be an array of uint32_t's, and
 -- that the length be the number of uint32_t's in the key */

uint32_t hashword(
const uint32_t *k,                   /* the key, an array of uint32_t values */
size_t          length,               /* the length of the key, in uint32_ts */
uint32_t        initval)         /* the previous hash, or an arbitrary value */
{
  uint32_t a,b,c;

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + (((uint32_t)length)<<2) + initval;

  /*------------------------------------------------- handle most of the key */
  while (length > 3)
  {
    a += k[0];
    b += k[1];
    c += k[2];
    mix(a,b,c);
    length -= 3;
    k += 3;
  }

  /*------------------------------------------- handle the last 3 uint32_t's */
  switch(length)                     /* all the case statements fall through */
  { 
  case 3 : c+=k[2];
  case 2 : b+=k[1];
  case 1 : a+=k[0];
    final(a,b,c);
  case 0:     /* case 0: nothing left to add */
    break;
  }
  /*------------------------------------------------------ report the result */
  return c;
}

int get_hash_index(int i, int j, int k, int l)
{
    uint32_t key[4], hashVal, initVal = 13;
    
    
    key[0] = i;
    key[1] = j;
    key[2] = k;
    key[3] = l;
   
    hashVal = hashword(key, 4, initVal); 
    return (hashVal & hashmask((int)(log(HTABLE_SIZE)/log(2))));
}


#ifndef __STATIC_SCHEDULE__

PQueue *pq_init(int nElems)
{
    PQueue *pq;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering pq_init() with %d nElems\n", nElems);
#endif

    pq = (PQueue *)malloc(sizeof(PQueue));
    if (pq == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    pq->size = 0;
    pq->avail = nElems;
    pq->pqData = (int *)malloc(sizeof(int)*(nElems+1));
    if (pq->pqData == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting pq_init() with %d nElems\n", nElems);
#endif
    return pq;
}

void pq_insert(PQueue *pq, int n)
{
    int i;

#ifdef __DEBUG1__
    fprintf(stderr, "Entering pq_insert() with n=%d and priority %.2lf\n", n, get_priority(n));
#endif
    assert(pq != NULL);
    for (i=1; i<=pq->size; i++)
    {
        if (pq->pqData[i] == n) 
        {
#ifdef __DEBUG2__
            fprintf(stderr, "Exiting pq_insert() with n=%d\n", n);
#endif
            return;
        }
    } 
    for (i=pq->size+1; (i>1) && (get_priority(pq->pqData[i/2]) < get_priority(n)); i/=2) 
    {
        pq->pqData[i] = pq->pqData[i/2];
    }
    pq->pqData[i] = n;
    pq->size++;
    assert(pq->size <= pq->avail);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting pq_insert() with n=%d\n", n);
#endif
    return;
} 

int pq_remove(PQueue *pq)
{
    int i, maxElem, lastElem, child;
#ifdef __DEBUG2__
    fprintf(stderr, "Entering pq_remove()\n");
#endif    
    assert(pq != NULL);
    if (pq->size == 0)
    {
#ifdef __DEBUG2__
        fprintf(stderr, "Task Queue is empty\n");
    	fprintf(stderr, "Exiting pq_remove()\n");
#endif    
        return -1;
    }
    maxElem = pq->pqData[1];
    lastElem = pq->pqData[pq->size];
    pq->size--;
    for (i=1; (i*2)<=pq->size; i=child)
    {
        child = i*2;
        if ((child != pq->size) && (get_priority(pq->pqData[child + 1]) > get_priority(pq->pqData[child])))
        {
            child++;
        }
        if (get_priority(lastElem) < get_priority(pq->pqData[child]))
        {
            pq->pqData[i] = pq->pqData[child];
        }
        else
        {
            break;
        }
    }
    pq->pqData[i] = lastElem;
#ifdef __DEBUG1__
    fprintf(stderr, "Exiting pq_remove() returning %d with priority %.2lf\n", maxElem, get_priority(maxElem));
#endif    
    return maxElem;
}

void pq_delete(PQueue *pq)
{
#ifdef __DEBUG2__
    fprintf(stderr, "Entering pq_delete()\n");
#endif    
    free(pq->pqData);
    free(pq);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting pq_delete()\n");
#endif    
    return;
}

double get_priority(int vertexId)
{
    return (dag.vertexList[vertexId].bLevel);
}    
#endif    
 
#ifdef __STATIC_SCHEDULE__

void init_schedule()
{
    int i, nProcs;

#pragma omp parallel
#pragma omp master
    {
        nProcs = omp_get_num_threads();
    }
#ifdef __DEBUG2__
    fprintf(stderr, "Entering init_schedule() with %d nProcs\n", nProcs);
#endif
    resource.nProcs = nProcs;
    resource.pFreeList = NULL;
    for (i=0; i<dag.nVertices; i++)
    {
        dag.vertexList[i].est = 0;
        dag.vertexList[i].ast = 0;
        dag.vertexList[i].complTime = 0;
    }
#ifdef __STATIC_SCHEDULE__
    schedule.taskQueue = (int**)malloc(sizeof(int*)*nProcs);
    if (schedule.taskQueue == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
           exit(-1);
    }
    schedule.nTasks = (int*)malloc(sizeof(int)*nProcs);
    if (schedule.nTasks == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
           exit(-1);
    }
    for (i=0; i<nProcs; i++)
    {
        schedule.taskQueue[i] = NULL;
        schedule.nTasks[i] = 0;
    }
#endif
    
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting init_schedule() with %d nProcs\n", nProcs);
#endif
    return;
}

void free_schedule(void)
{
    ResourceFreeList *pTmp1, *pTmp2;
    int i;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering free_schedule()\n");
#endif

    pTmp1 = resource.pFreeList;
    while (pTmp1 != NULL)
    {
        pTmp2 = pTmp1->pNext;
        if (pTmp1->freeSet)
        {
            free(pTmp1->freeSet);
        }
        free(pTmp1);
        pTmp1 = pTmp2;
    }
    resource.pFreeList = NULL;
    free(schedule.nTasks);
    for (i = 0; i < resource.nProcs; i++)
    {
        if (schedule.taskQueue[i])
        {
            free(schedule.taskQueue[i]);
        }
    }
    free(schedule.taskQueue);
    schedule.taskQueue = NULL;
    schedule.nTasks = NULL;
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting free_schedule()\n");
#endif
    return;
}

int is_present(int item, int *arr, int nElems)
{
    int i;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering is_present() \n");
#endif
    for (i=0; i<nElems; i++)
    {
        if (arr[i] == item)
        {
            return 1;
        }
    }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting is_present() \n");
#endif
    return 0;
}

void create_reservation(double ast, double dur, int proc)
{
    ResourceFreeList *pFreeList, *pPrevPtr, *newEntry;
    int i, *tmpProcSet, index;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering create_reservation() with %.2lf ast, %.2lf dur, %d proc\n", ast, dur, proc);
#endif

    pFreeList = resource.pFreeList;
    pPrevPtr = NULL;
    assert(dur > 0);
    while ((pFreeList) && (pFreeList->endTime <= ast))
    {
        pPrevPtr = pFreeList;
        pFreeList = pFreeList->pNext;
    }
    if (!pFreeList)
    {
        if ((!pPrevPtr) && (ast > 0))
        {
            newEntry = (ResourceFreeList*)malloc(sizeof(ResourceFreeList));
            if (newEntry == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                exit(-1);
            }
            newEntry->startTime = 0;
            newEntry->endTime = ast;
            newEntry->nFree = resource.nProcs;
            newEntry->freeSet = (int*)malloc(sizeof(int)*resource.nProcs);
            if (newEntry->freeSet == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                exit(-1);
            }
            for (i=0; i<resource.nProcs; i++)
            {
                newEntry->freeSet[i] = i;
            }
            newEntry->pNext = NULL;
            pPrevPtr = newEntry;
            resource.pFreeList = newEntry;
        }
        else if ((pPrevPtr) && ((float)pPrevPtr->endTime < (float)ast))
        {
            newEntry = (ResourceFreeList*)malloc(sizeof(ResourceFreeList));
            if (newEntry == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                exit(-1);
            }
            newEntry->startTime = pPrevPtr->endTime;
            newEntry->endTime = ast;
            newEntry->nFree = resource.nProcs;
            newEntry->freeSet = (int*)malloc(sizeof(int)*resource.nProcs);
            if (newEntry->freeSet == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                exit(-1);
            }
            for (i=0; i<resource.nProcs; i++)
            {
                newEntry->freeSet[i] = i;
            }
            newEntry->pNext = NULL;
            pPrevPtr->pNext = newEntry;
            pPrevPtr = newEntry;
        }
        newEntry = (ResourceFreeList*)malloc(sizeof(ResourceFreeList));
        if (newEntry == NULL)
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
               exit(-1);
        }
        newEntry->startTime = ast;
        newEntry->endTime = ast + dur;
        newEntry->nFree = resource.nProcs - 1;
        assert(newEntry->nFree >= 0);
        if (newEntry->nFree > 0)
        {
            newEntry->freeSet = (int*)malloc(sizeof(int)*newEntry->nFree);
            if (newEntry->freeSet == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                   exit(-1);
            }
            index = 0;
            for (i=0; i<resource.nProcs; i++)
            {
                if (i != proc)
                {
                    newEntry->freeSet[index] = i;
                    index++;
                }
            }
        }
        else
        {
            newEntry->freeSet = NULL;
        }
        newEntry->pNext = NULL;
        if (pPrevPtr)
        {
            pPrevPtr->pNext = newEntry;
        }
        else
        {
            resource.pFreeList = newEntry;
        }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting create_reservation() with %.2lf ast, %.2lf dur, %d proc\n", ast, dur, proc);
#endif
        return;
    }
    if (ast > pFreeList->startTime)
    {
        newEntry = (ResourceFreeList*)malloc(sizeof(ResourceFreeList));
        if (newEntry == NULL)
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
               exit(-1);
        }
        newEntry->startTime = pFreeList->startTime;
        newEntry->endTime = ast;
        newEntry->nFree = pFreeList->nFree;
        assert((newEntry->nFree >= 0) && (newEntry->nFree <= resource.nProcs));
        newEntry->freeSet = (int*)malloc(sizeof(int)*newEntry->nFree);
        if (newEntry->freeSet == NULL)    
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
               exit(-1);
        }
        for (i=0; i<newEntry->nFree; i++)
        {
            newEntry->freeSet[i] = pFreeList->freeSet[i];
        }
        if (pPrevPtr)
        {
            pPrevPtr->pNext = newEntry;
        }
        else
        {
            resource.pFreeList = newEntry;
        }
        newEntry->pNext = pFreeList;
        pFreeList->startTime = ast;
        pPrevPtr = newEntry;
    }
    while ((pFreeList) && ((float)(ast + dur) >= (float)pFreeList->endTime))
    {
        if ((pFreeList->nFree - 1) > 0)
        {
            tmpProcSet = (int*)malloc(sizeof(int)*(pFreeList->nFree - 1));
            if (tmpProcSet == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                   exit(-1);
            }
            index = 0;
            for (i=0; i<pFreeList->nFree; i++)
            {
                if (pFreeList->freeSet[i] != proc)
                {
                    tmpProcSet[index] = pFreeList->freeSet[i];
                    index++;
                }
            }
        }
        else
        {
            tmpProcSet = NULL;
        }
        if (pFreeList->freeSet)
        {
            free(pFreeList->freeSet);
        }
        pFreeList->freeSet = tmpProcSet;
        pFreeList->nFree--;
        assert((pFreeList->nFree >= 0) && (pFreeList->nFree <= resource.nProcs));
        pPrevPtr = pFreeList;
        pFreeList = pFreeList->pNext;
    }
    if ((pPrevPtr) && ((float) pPrevPtr->endTime == (float) (ast + dur)))
    {
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting create_reservation() with %.2lf ast, %.2lf dur, %d proc\n", ast, dur, proc);
#endif
        return;
    }
    if (!pFreeList)
    {
        assert(pPrevPtr != NULL);
        newEntry = (ResourceFreeList*)malloc(sizeof(ResourceFreeList));
        if (newEntry == NULL)
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
               exit(-1);
        }
        newEntry->startTime = pPrevPtr->endTime;
        newEntry->endTime = ast + dur;
        newEntry->nFree = resource.nProcs - 1;
        if ((resource.nProcs - 1) > 0)
        {
            newEntry->freeSet = (int*)malloc(sizeof(int)*(resource.nProcs - 1));
            if (newEntry->freeSet == NULL)
            {
                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                   exit(-1);
            }
            index = 0;
            for (i=0; i<resource.nProcs; i++)
            {
                if (i != proc)
                {
                    newEntry->freeSet[index] = i;
                    index++;
                }
            }
        }
        else
        {
            newEntry->freeSet = NULL;
        }
        newEntry->pNext = NULL;
        if (pPrevPtr)
        {
            pPrevPtr->pNext = newEntry;
        }
        else
        {
            resource.pFreeList = newEntry;
        }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting create_reservation() with %.2lf ast, %.2lf dur, %d proc\n", ast, dur, proc);
#endif
        return;
    }
    newEntry = (ResourceFreeList*)malloc(sizeof(ResourceFreeList));
    if (newEntry == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
           exit(-1);
    }
    newEntry->startTime = pFreeList->startTime;
    newEntry->endTime = ast + dur;
    newEntry->nFree = pFreeList->nFree - 1;
    if ((pFreeList->nFree - -1) > 0)
    {
        newEntry->freeSet = (int*)malloc(sizeof(int)*(pFreeList->nFree - 1));
        if (newEntry->freeSet == NULL)
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
               exit(-1);
        }
        index = 0;
        for (i=0; i<pFreeList->nFree; i++)
        {
            if (pFreeList->freeSet[i] != proc)
            {
                newEntry->freeSet[index] = pFreeList->freeSet[i];
                index++;
            }
        }
    }
    else
    {
        newEntry->freeSet = NULL;
    }
    if (pPrevPtr)
    {
        pPrevPtr->pNext = newEntry;
    }
    else
    {
        resource.pFreeList = newEntry;
    }
    newEntry->pNext = pFreeList;
    pFreeList->startTime = ast + dur;
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting create_reservation() with %.2lf ast, %.2lf dur, %d proc\n", ast, dur, proc);
#endif
    return;
}

double find_min_starttime_for_task_on_proc(int proc, double est, double dur)
{
    double timeNeeded, ast;
    ResourceFreeList *pFreeList;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering find_min_starttime_for_task_on_proc() with %d proc, %.2lf est, %.2lf dur\n", proc, est, dur);
#endif

    timeNeeded = dur;
    pFreeList = resource.pFreeList;
    ast = est;
    while ((pFreeList) && ((float)pFreeList->endTime <= (float)est))    
    {
        pFreeList = pFreeList->pNext;
    }
    if (!pFreeList)
    {
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting find_min_starttime_for_task_on_proc() with %d proc, %.2lf est, %.2lf dur\n", proc, est, dur);
#endif
        return ast;
    }
    if (pFreeList->nFree > 0)
    {
        if (is_present(proc, pFreeList->freeSet, pFreeList->nFree))
        {
            timeNeeded -= (pFreeList->endTime - est);
        }
        else
        {
            ast = pFreeList->endTime;
        }
    }
    else
    {
        ast = pFreeList->endTime;
    }
    while ((pFreeList->pNext) && (timeNeeded > 0))
    {
        pFreeList = pFreeList->pNext;
        assert(pFreeList->nFree >= 0);
        if (pFreeList->nFree == 0)
        {
            ast = pFreeList->endTime;
            timeNeeded = dur;
            continue;
        }
        else
        {
            if (is_present(proc, pFreeList->freeSet, pFreeList->nFree))
            {
                timeNeeded -= (pFreeList->endTime - pFreeList->startTime);
            }
            else
            {
                ast = pFreeList->endTime;
                timeNeeded = dur;
            }
        }
    }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting find_min_starttime_for_task_on_proc() with %d proc, %.2lf est, %.2lf dur\n", proc, est, dur);
#endif
    return ast;
}
#endif

void sort_priorities(int *sortList, int nEntries, int sortFlag)
{
    int i, *initList;
    

#ifdef __DEBUG2__
    fprintf(stderr, "Entering sort_priorities() with %d nEntries and %d sortFlag\n", nEntries, sortFlag);
#endif

    initList = (int *)malloc(sizeof(int)*nEntries);
    if (initList == NULL)
    {    
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    for (i=0; i<nEntries; i++)
    {
        initList[i] = sortList[i];
    }
    merge_sort(initList, sortList, 0, nEntries - 1, sortFlag);
    free(initList);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting sort_priorities() with %d nEntries and %d sortFlag\n", nEntries, sortFlag);
#endif
    return;
}

void merge_sort(int *initList, int *sortList, int left, int right, int sortFlag)
{
    int mid;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering merge_sort()\n");
#endif

    if (right > left)
    {
        mid = (right + left)/2;
        merge_sort(initList, sortList, left, mid, sortFlag);
        merge_sort(initList, sortList, mid+1, right, sortFlag);
        merge(initList, sortList, left, mid+1, right, sortFlag);
    }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting merge_sort()\n");
#endif
    return;
}

void merge(int *list1, int *list2, int left, int mid, int right, int sortFlag)
{
    int i, leftEnd, nElements, tmpPos;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering merge()\n");
#endif
    leftEnd = mid - 1;
    tmpPos = left;
    nElements = right - left + 1;
    while ((left <= leftEnd) && (mid <= right))
    {
        switch (sortFlag)
        {
        case INC_AST:
                        if ((dag.vertexList[list1[left]].vertexWt != 0) && (dag.vertexList[list1[mid]].vertexWt != 0))
                        {
                            assert((float)dag.vertexList[list1[left]].ast != (float)dag.vertexList[list1[mid]].ast);
                        }
                        if (dag.vertexList[list1[left]].ast < dag.vertexList[list1[mid]].ast)
                        {
                            list2[tmpPos] = list1[left];
                            tmpPos++;
                            left++;
                        }
                        else if (dag.vertexList[list1[left]].ast > dag.vertexList[list1[mid]].ast)
                        {
                            list2[tmpPos] = list1[mid];
                            tmpPos++;
                            mid++;
                        }
			else if (dag.vertexList[list1[left]].vertexWt == 0)
                        {
                            list2[tmpPos] = list1[left];
                            tmpPos++;
                            left++;
			}
			else
                        {
                            list2[tmpPos] = list1[mid];
                            tmpPos++;
                            mid++;
                        }
                        break;
        case DEC_BLEVEL:
                        if (dag.vertexList[list1[left]].bLevel >= dag.vertexList[list1[mid]].bLevel)
                        {
                            list2[tmpPos] = list1[left];
                            tmpPos++;
                            left++;
                        }
                        else
                        {
                            list2[tmpPos] = list1[mid];
                            tmpPos++;
                            mid++;
                        }
                        break;
        case INC_ALAP:
                        if (dag.vertexList[list1[left]].bLevel > dag.vertexList[list1[mid]].bLevel)
                        {
                            list2[tmpPos] = list1[left];
                            tmpPos++;
                            left++;
                        }
                        else if (dag.vertexList[list1[left]].bLevel == dag.vertexList[list1[mid]].bLevel)
                        {
                            int i, *cList1, *cList2, cIndex1, cIndex2;
                            Edge *pEdge;

                            cList1 = (int *)malloc(dag.vertexList[list1[left]].nChildren * sizeof(int));
                            if (cList1 == NULL)
                            {
                                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                                exit(-1);
                            }                        
                            cList2 = (int *)malloc(dag.vertexList[list1[mid]].nChildren * sizeof(int));
                            if (cList2 == NULL)
                            {
                                fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
                                exit(-1);
                            }                        
                            pEdge = dag.vertexList[list1[left]].childList;
                            for (i=0; i<dag.vertexList[list1[left]].nChildren; i++)
                            {
                                assert(pEdge != NULL);
                                cList1[i] = pEdge->dvertexId;
                                pEdge = pEdge->pNext;
                            } 
                            pEdge = dag.vertexList[list1[mid]].childList;
                            for (i=0; i<dag.vertexList[list1[mid]].nChildren; i++)
                            {
                                assert(pEdge != NULL);
                                cList2[i] = pEdge->dvertexId;
                                pEdge = pEdge->pNext;
                            } 
                            sort_priorities(cList1, dag.vertexList[list1[left]].nChildren, DEC_BLEVEL);
                            sort_priorities(cList2, dag.vertexList[list1[mid]].nChildren, DEC_BLEVEL);
                            cIndex1 = 0; cIndex2 = 0;
                            while ((cIndex1 < dag.vertexList[list1[left]].nChildren) && (cIndex2 < dag.vertexList[list1[mid]].nChildren))
                            {
                                if (dag.vertexList[cList1[cIndex1]].bLevel > dag.vertexList[cList2[cIndex2]].bLevel)
                                {
                                    list2[tmpPos] = list1[left];
                                    tmpPos++;
                                    left++;
                                    break;
                                }        
                                if (dag.vertexList[cList1[cIndex1]].bLevel < dag.vertexList[cList2[cIndex2]].bLevel)
                                {
                                    list2[tmpPos] = list1[mid];
                                    tmpPos++;
                                    mid++;
                                    break;
                                }
                                else
                                {
                                    cIndex1++; cIndex2++;
                                }
                            }
                            if (cIndex1 == dag.vertexList[list1[left]].nChildren)
                            {
                                list2[tmpPos] = list1[mid];
                                tmpPos++;
                                mid++;
                            }
                            else if (cIndex2 == dag.vertexList[list1[mid]].nChildren)
                            {
                                list2[tmpPos] = list1[left];
                                tmpPos++;
                                left++;
                            }        
                            free(cList1); free(cList2);
                        }
                        else
                        {
                            list2[tmpPos] = list1[mid];
                            tmpPos++;
                            mid++;
                        }
                        break;
        }
    }
    while (left <= leftEnd)
    {
        list2[tmpPos] = list1[left];
        left++;
        tmpPos++;
    }
    while (mid <= right)
    {
        list2[tmpPos] = list1[mid];
        mid++;
        tmpPos++;
    }
    for (i=0; i<nElements; i++)
    {
        list1[right] = list2[right];
        right--;
    }
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting merge()\n");
#endif
    return;
}

#ifdef __STATIC_SCHEDULE__

void schedule_dag_using_mcp(void)
{
    int *sortedVertexList, i, j;
    Edge *pEdge;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering schedule_dag_using_mcp() \n");
#endif
    calculate_tlevel_blevel();
#ifdef __DEBUG1__
    fprintf(stderr, "Scheduling the DAG using MCP... \n");
    fprintf(stderr, "Calculating top and bottom levels... \n");
#endif
    sortedVertexList = (int *)malloc(sizeof(int)*dag.nVertices);
    if (sortedVertexList == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    for (i=0; i<dag.nVertices; i++)
    {
        sortedVertexList[i] = i;
    }
    sort_priorities(sortedVertexList, dag.nVertices, INC_ALAP);
#ifdef __DEBUG1__ 
    fprintf(stderr, "Sorting tasks in increasing order of ALAP time... \n");
#endif
    for (i=0; i<dag.nVertices; i++)
    {
        int vertexId;

        vertexId = sortedVertexList[i];
        assert((vertexId >= 0) && (vertexId < dag.nVertices));
        schedule_task(vertexId);
#ifdef __DEBUG1__ 
    fprintf(stderr, "Scheduling task %d on proc %d with %.2lf est, %.2lf ast, %.2lf compl time... \n", vertexId, dag.vertexList[vertexId].procId, dag.vertexList[vertexId].est, dag.vertexList[vertexId].ast, dag.vertexList[vertexId].complTime);
#endif
        pEdge = dag.vertexList[vertexId].childList;
        for (j = 0; j < dag.vertexList[vertexId].nChildren; j++)
        {
            assert(pEdge != NULL);
            dag.vertexList[pEdge->dvertexId].est = (dag.vertexList[pEdge->dvertexId].est < dag.vertexList[vertexId].complTime)?dag.vertexList[vertexId].complTime: dag.vertexList[pEdge->dvertexId].est;
            pEdge = pEdge->pNext;
        } 
    }
    populate_task_queue();
    free(sortedVertexList);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting schedule_dag_using_mcp() \n");
#endif 
    return;
}

void schedule_task(int vertexId)
{
    int i, *reuse;
    double ast, maxDataAccessTime, minStartTime;
    Edge *pEdge;

#ifdef __DEBUG1__ 
    fprintf(stderr, "Entering schedule_task() for vertexId %d\n", vertexId);
#endif
    if ((float)dag.vertexList[vertexId].vertexWt == 0)
    {
        assert((dag.vertexList[vertexId].nParents == 0) || (dag.vertexList[vertexId].nChildren == 0));
        dag.vertexList[vertexId].ast = dag.vertexList[vertexId].est;
        dag.vertexList[vertexId].procId = 0;
        dag.vertexList[vertexId].complTime = dag.vertexList[vertexId].ast + dag.vertexList[vertexId].vertexWt;
        schedule.nTasks[dag.vertexList[vertexId].procId]++;
#ifdef __DEBUG2__
        fprintf(stderr, "Exiting schedule_task() for vertexId %d\n", vertexId);
#endif 
        return;
    }    
    reuse = (int*)calloc(resource.nProcs, sizeof(int));
    if (reuse == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    maxDataAccessTime = 0;
    pEdge = dag.vertexList[vertexId].parentList;
    for (i=0; i<dag.vertexList[vertexId].nParents; i++)
    {
        int proc;

        assert(pEdge != NULL);
        proc = dag.vertexList[pEdge->svertexId].procId;
        reuse[proc] += pEdge->edgeWt;
        maxDataAccessTime += pEdge->edgeWt;
        pEdge = pEdge->pNext;
    }
    minStartTime = MAX_VALUE;
    for (i=0; i<resource.nProcs; i++)
    {
        ast = find_min_starttime_for_task_on_proc(i, dag.vertexList[vertexId].est, (maxDataAccessTime - reuse[i]) + dag.vertexList[vertexId].vertexWt);
        if (minStartTime > (ast + (maxDataAccessTime - reuse[i])))
        {
            minStartTime = ast + (maxDataAccessTime - reuse[i]);
            dag.vertexList[vertexId].ast = ast;
            dag.vertexList[vertexId].procId = i;
            dag.vertexList[vertexId].complTime = minStartTime + dag.vertexList[vertexId].vertexWt;
        }
    }
#ifdef __DEBUG1__
    fprintf(stderr, "Scheduling vertexId %d on proc %d at time %.2lf\n", vertexId, dag.vertexList[vertexId].procId, dag.vertexList[vertexId].ast);
#endif
    create_reservation(dag.vertexList[vertexId].ast, (maxDataAccessTime - reuse[dag.vertexList[vertexId].procId]) + dag.vertexList[vertexId].vertexWt, dag.vertexList[vertexId].procId);
    schedule.nTasks[dag.vertexList[vertexId].procId]++;
    free(reuse);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting schedule_task() for vertexId %d\n", vertexId);
#endif 
    return;
}

void populate_task_queue(void)
{
     int i, *index;

#ifdef __DEBUG2__
    fprintf(stderr, "Entering populate_task_queue()\n");
#endif
    index = (int*)malloc(sizeof(int)*resource.nProcs);
    if (index == NULL)
    {
        fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
        exit(-1);
    }
    for (i=0; i<resource.nProcs; i++)
    {
        schedule.taskQueue[i] = (int*)malloc(sizeof(int)*schedule.nTasks[i]);
        if (schedule.taskQueue[i] == NULL)
        {
            fprintf(stderr, "%s:%d: Run out of memory. Exiting program...\n", __FILE__, __LINE__);
            exit(-1);
        }
        index[i] = 0;
    }
    for (i=0; i<dag.nVertices; i++)
    {
        int proc;

        proc = dag.vertexList[i].procId;
        assert((proc >= 0) && (proc < resource.nProcs));
        schedule.taskQueue[proc][index[proc]] = i;
        index[proc]++;
    }
    for (i=0; i<resource.nProcs; i++)
    {
        sort_priorities(schedule.taskQueue[i], schedule.nTasks[i], INC_AST);
    }
    free(index);
#ifdef __DEBUG2__
    fprintf(stderr, "Exiting populate_task_queue()\n");
#endif
}
#endif
