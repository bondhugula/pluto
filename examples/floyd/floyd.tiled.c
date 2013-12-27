#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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



  int t1, t2, t3, t4, t5;

 register int lbv, ubv;

/* Start of CLooG code */
if (NUM_NODES >= 1) {
  for (t1=0;t1<=NUM_NODES-1;t1++) {
    for (t2=0;t2<=floord(NUM_NODES-1,32);t2++) {
      for (t3=0;t3<=floord(NUM_NODES-1,32);t3++) {
        for (t4=32*t2;t4<=min(NUM_NODES-1,32*t2+31);t4++) {
          for (t5=32*t3;t5<=min(NUM_NODES-1,32*t3+31);t5++) {
            pathDistanceMatrix[t4][t5]=((pathDistanceMatrix[t4][t1]+pathDistanceMatrix[t1][t5])<pathDistanceMatrix[t4][t5])?(pathDistanceMatrix[t4][t1]+pathDistanceMatrix[t1][t5]):pathDistanceMatrix[t4][t5];;
          }
        }
      }
    }
  }
}
/* End of CLooG code */

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
