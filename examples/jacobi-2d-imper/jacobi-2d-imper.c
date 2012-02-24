#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#define N 2000
#define T 2000

#pragma declarations
double a[N][N];
double b[N][N];
#pragma enddeclarations

#include "util.h"

int main()
{
    int t, i, j;
    double t_start, t_end;

    init_array();

    IF_TIME(t_start = rtclock());

#pragma scop
    for (t=0; t<T; t++) {
        for (i=2; i<N-1; i++) {
            for (j=2; j<N-1; j++) {
                b[i][j]= 0.2*(a[i][j]+a[i][j-1]+a[i][1+j]+a[1+i][j]+a[i-1][j]);
            }
        }
        for (i=2; i<N-1; i++) {
            for (j=2; j<N-1; j++)   {
                a[i][j]=b[i][j];
            }
        }
    }
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r")) {
        print_array();
    }
    return 0;
}
