#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <assert.h>

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "decls.h"
#include "util.h"

double t_start, t_end;

int main()
{
	int i, j, k, l, m, n, o, t;

	init_array();

#ifdef PERFCTR
	PERF_INIT; 
#endif

	IF_TIME(t_start = rtclock());


    /* pluto start (M) */
    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++) {
            temp2d[i][j] = 0.0;
            for (k = 0; k < M; k++)
                temp2d[i][j] = temp2d[i][j] + block[i][k] * cos1[j][k];
        }

    for (i = 0; i < M; i++)
        for (j = 0; j < M; j++) {
            sum2[i][j] = 0.0;
            for (k = 0; k < M; k++)
                sum2[i][j] = sum2[i][j] + cos1[i][k] * temp2d[k][j];
            block[i][j] = (sum2[i][j]);
        }
    /* pluto end */

    /* End of CLooG code */
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
