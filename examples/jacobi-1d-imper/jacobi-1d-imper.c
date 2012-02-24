#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"

#define N 1000
#define T 50

#pragma declarations
double a[N];
double b[N];
#pragma enddeclarations

double t_start, t_end;

void init_array()
{
    int j;

    for (j=0; j<N; j++) {
        a[j] = ((double)j)/N;
    }
}


void print_array()
{
    int j;

    for (j=0; j<N; j++) {
        fprintf(stdout, "%lf ", a[j]);
        if (j%80 == 20) fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}


int main()
{
    int t, i, j;

    init_array();

    IF_TIME(t_start = rtclock());

#pragma scop
    for (t = 0; t < T; t++) {
        for (i = 2; i < N - 1; i++) {
            b[i] = 0.33333 * (a[i-1] + a[i] + a[i + 1]);
        }
        for (j = 2; j < N - 1; j++) {
            a[j] = b[j];
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
