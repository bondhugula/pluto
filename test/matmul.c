#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define M 2048
#define N 2048
#define K 2048
#define alpha 1
#define beta 1
double A[M][K+13];
double B[K][N+13];
double C[M][N+13];

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{
    int i, j;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = (i + j);
            B[i][j] = (double)(i*j);
            C[i][j] = 0.0;
        }
    }
}


void print_array()
{
    int i, j;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fprintf(stderr, "%lf ", C[i][j]);
            if (j%80 == 79) fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    }
}

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

int main()
{
    int i, j, k;
    register double s;


    init_array();

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

#pragma scop
    /* pluto start (M,N,K) */
    for(i=0; i<M; i++)
        for(j=0; j<N; j++)  
            for(k=0; k<K; k++)
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
    /* pluto end */
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(printf("%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

#ifdef TEST
    print_array();
#endif
    return 0;
}
