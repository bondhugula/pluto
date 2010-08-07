#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#include "decls.h"
#include "util.h"

int main()
{
    int t, i, j, k, l;

    double t_start, t_end;

    init_array();

	IF_TIME(t_start = rtclock());

    /* pluto start (N,T) */
    for (t=1; t<=T; t++){
	    for (i=1; i<=N-1; i++)
		    e[i] = e[i] - coeff1*(h[i]-h[i-1]);
	    for (i=0; i<=N-1; i++)
		    h[i] = h[i] - coeff2*(e[i+1]-e[i]);
    }
    /* pluto end */

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
        print_array();
    }

    return 0;
}
