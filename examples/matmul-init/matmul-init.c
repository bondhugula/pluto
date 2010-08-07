#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "decls.h"

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "util.h"

int main()
{
	int i, j, k;
    register double s;
    double t_start, t_end;

	init_array();

#ifdef PERFCTR
	PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

#pragma scop
    for(i=0; i<N; i++)  {
        for(j=0; j<N; j++)  {
            C[i][j] = 0;
            for(k=0; k<N; k++)
                C[i][j] = C[i][j] + A[i][k] * B[k][j];
        }
    }
#pragma endscop


    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

#ifdef TEST
    print_array();
#endif
    return 0;
}
