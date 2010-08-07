#include <stdio.h>

#include "decls.h"
#include "util.h"

double t_start, t_end;

main()
{
    int t, p, q, r, s;

    init_array();

    IF_TIME(t_start = rtclock());

#ifdef TIME
    for (t=0; t<1000; t++)  {
#endif

#pragma scop
        for( r = 0; r< N; r++)  {
            for( q = 0; q< N; q++)  {
                for( p = 0; p< N; p++)  {
                    sum[r][q][p] = 0;
                    for( s = 0; s< N; s++)  {
                        sum[r][q][p] = sum[r][q][p] + A[r][q][s]*C4[s][p];
                    }
                }
                for( p = 0; p< N; p++)  {
                    A[r][q][p] = sum[r][q][p];
                }
            }
        }
#pragma endscop

#ifdef TIME
    }
#endif

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef TEST
    print_array();
#endif

    return 0;

}
