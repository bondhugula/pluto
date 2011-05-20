#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "decls.h"
#include "util.h"

#define TIME 1

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void strsm(long N) {
	int i,j,k;

#pragma scop
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			for (k=i+1; k<N; k++) {
				if (k == i+1) b[j][i] /= a[i][i];
                b[j][k] -= a[i][k] * b[j][i];
            }
        }
    }
#pragma endscop
}


int main()
{
    double t_start, t_end;
    long N=NMAX;
    int i,j;

    IF_TIME(t_start = rtclock());
    strsm(N);
    IF_TIME(t_end = rtclock());

    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                fprintf(stdout, "%lf ",  b[i][j]);
            }
            fprintf(stdout, "\n");
        }
    }

    return 0;
}
