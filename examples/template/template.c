#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#define N 1000
double f[N][N+13];

#include "util.h"

double t_start, t_end;

int main()
{
    int i, j, k, t;

    init_array() ;

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

    /* pluto start (N) */
#pragma scop
    for (i=1; i<=N-2; i++)  {
        for (j=1; j<=N-2; j++)  {
            f[i][j] = f[j][i] + f[i][j-1];
        }
    }
#pragma endscop
    /* pluto end */

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
        print_array();
    }

    return 0;
}
