#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#define N 1000
#define T 500

#pragma declarations
double a[N][N];
#pragma enddeclarations

#include "util.h"

int main()
{
    int i, j, k, t;

    double t_start, t_end;

    init_array() ;

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

#pragma scop
    for (t=0; t<=T-1; t++)  {
        for (i=1; i<=N-2; i++)  {
            for (j=1; j<=N-2; j++)  {
                a[i][j] = (a[i-1][j-1] + a[i-1][j] + a[i-1][j+1] 
                        + a[i][j-1] + a[i][j] + a[i][j+1]
                        + a[i+1][j-1] + a[i+1][j] + a[i+1][j+1])/9.0;
            }
        }
    }
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

    if (fopen(".test", "r")) {
        print_array();
    }
    return 0;
}
