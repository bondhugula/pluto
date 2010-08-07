#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <omp.h>

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "decls.h"

#include "util.h"

int main()
{
    int i, j, k, l, m, n, t;

    int i1, i2;

    double t_start, t_end;

    init_array();

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());


    /* pluto start (T,N) */
    for (t = 0; t < T; t++) {

        for (i1=0; i1<N; i1++) {
            for (i2 = 1; i2 < N; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
            }
        }

        for (i1=0; i1<N; i1++) {
                X[i1][N-1] = X[i1][N-1] / B[i1][N-1];
        }

        for (i1=0; i1<N; i1++) {
            for (i2 = 0; i2 <= N-2; i2++) {
                X[i1][N-i2-2] = (X[i1][N-2-i2] - X[i1][N-2-i2-1] * A[i1][N-i2-3]) / B[i1][N-3-i2];
            }
        }

        for (i1=1; i1<N; i1++) {
            for (i2 = 0; i2 < N; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
            }
        }

        for (i2=0; i2<N; i2++) {
                X[N-1][i2] = X[N-1][i2] / B[N-1][i2];
        }

        for (i1 = 0; i1 <= N-2; i1++) {
            for (i2=0; i2<N; i2++) {
                X[N-2-i1][i2] = (X[N-2-i1][i2] - X[N-i1-3][i2] * A[N-3-i1][i2]) / B[N-2-i1][i2];
            }
        }

    }
    /* pluto end */
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
