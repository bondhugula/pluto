#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "decls.h"

#include "util.h"

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

double t_start, t_end;

int main()
{
	assert(nx >= 100);
	assert(ny >= 100);
	assert(nz >= 100);
	int i, j, k, l, m, n, t;

	init_array() ;

#ifdef PERFCTR
	PERF_INIT; 
#endif

	IF_TIME(t_start = rtclock());
#define reps 10

    /* pluto start (nx,ny,nz) */
    for (j = 4; j <= ny+9-2; j++)
        for (i = 4; i <= nx+9-3; i++)
            for (k = 4; k <= nz+9-3; k++)
                ab[j][i][k] = (f60 * (a[j-1][i][k] + a[j][i][k]) + f61 
                        * (a[j-2][i][k] + a[j+1][i][k]) + f62 * (a[j-3][i][k] 
                            + a[j+2][i][k])) * thirddtbydy * uyb[j][i][k];

    for (j = 4; j <= ny+9-3; j++)
        for (i = 4; i <= nx+9-2; i++)
            for (k = 4; k <= nz+9-3; k++)
                al[j][i][k] = (f60 * (a[j][i-1][k] + a[j][i][k]) + f61 
                        * (a[j][i-2][k] + a[j][i+1][k]) + f62 * (a[j][i-3][k] + a[j][i+2][k]))
                    * thirddtbydx * uxl[j][i][k];

    for (j = 4; j <= ny+9-3; j++)
        for (i = 4; i <= nx+9-3; i++)
            for (k = 4; k <= nz+9-2; k++)
                af[j][i][k] = (f60 * (a[j][i][k-1] + a[j][i][k]) + f61 
                        * (a[j][i][k-2] + a[j][i][k+1]) + f62 * (a[j][i][k-3] + a[j][i][k+2]))
                    * thirddtbydz * uzf[j][i][k];


    for (j = 4; j <= ny+9-3; j++)
        for (i = 4; i <= nx+9-3; i++)
            for (k = 4; k <= nz+9-3; k++)
                athird[j][i][k] = a[j][i][k] + (al[j][i+1][k] - al[j][i][k])
                    + (ab[j+1][i][k] - ab[j][i][k]) + (af[j][i][k+1] - af[j][i][k]);

    /* pluto end */
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
