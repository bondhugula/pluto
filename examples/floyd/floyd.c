#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif


double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

double t_start, t_end;

#define NUM_NODES 2048
#pragma declarations
double pathDistanceMatrix[NUM_NODES][NUM_NODES];
#pragma enddeclarations

void init_array()
{
    int i, j;

    srand(123);
    for(i= 0; i < NUM_NODES; i++) {
        for(j= 0; j < NUM_NODES; j++) {
            double ra = 1.0 + ((double)NUM_NODES * rand() / (RAND_MAX + 1.0));
            pathDistanceMatrix[i][j] = ra;
        }
    }
    for(i= 0; i < NUM_NODES; i++) {
        pathDistanceMatrix[i][i] = 0;
    }
}


void print_array()
{
    int i, j;
    for(i= 0; i < NUM_NODES; i++) {
        for(j= 0; j < NUM_NODES; j++) {
            fprintf(stderr, "%lf ", pathDistanceMatrix[i][j]); 
        }
        fprintf(stderr, "\n");
    }
}



int main()
{
    int i, j, k, x, y;
    unsigned int distanceYtoX, distanceYtoK, distanceKtoX;
    /*
     * pathDistanceMatrix is the adjacency matrix (square) with
     * dimension length equal to number of nodes in the graph.
     */
    unsigned int width = NUM_NODES;
    unsigned int yXwidth;

    init_array();

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

#pragma scop
    for(k=0; k < NUM_NODES; k++)
    {
        for(y=0; y < NUM_NODES; y++)
        {
            for(x=0; x < NUM_NODES; x++)
            {
                pathDistanceMatrix[y][x] = ((pathDistanceMatrix[y][k] + pathDistanceMatrix[k][x]) < pathDistanceMatrix[y][x]) ? (pathDistanceMatrix[y][k] + pathDistanceMatrix[k][x]):pathDistanceMatrix[y][x];
            }
        }
    }
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "time = %0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

    if (fopen(".test", "r")) {
#ifdef MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_array();
#endif
    }

    return 0;
}
