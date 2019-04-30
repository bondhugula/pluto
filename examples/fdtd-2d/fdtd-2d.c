#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

#include <assert.h>

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#define tmax 128
#define nx 2048
#define ny 2048

#pragma declarations
double ex[nx][ny+1];
double ey[nx+1][ny];
double hz[nx][ny];
#pragma enddeclarations


#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

void init_array()
{
    int i, j;

    for (i=0; i<nx+1; i++)  {
        for (j=0; j<ny; j++)  {
            ey[i][j] = 0;
        }
    }

    for (i=0; i<nx; i++)  {
        for (j=0; j<ny+1; j++)  {
            ex[i][j] = 0;
        }
    }

    for (j=0; j<ny; j++)  {
        ey[0][j] = ((double)j)/ny;
    }

    for (i=0; i<nx; i++)    {
        for (j=0; j<ny; j++)  {
            hz[i][j] = 0;
        }
    }
}


void print_array()
{
    int i, j;

    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            fprintf(stderr, "%lf ", hz[i][j]);
            if (j%80 == 20) fprintf(stderr, "\n");
        }
    }
    fprintf(stdout, "\n");
}

double t_start, t_end;

int main()
{
	int t, i, j, k, l, m, n;

	init_array() ;

#ifdef PERFCTR
    PERF_INIT;
#endif

	IF_TIME(t_start = rtclock());

#pragma scop
    for(t=0; t<tmax; t++)  {
        for (j=0; j<ny; j++)
            ey[0][j] = t;
        for (i=1; i<nx; i++)
            for (j=0; j<ny; j++)
                ey[i][j] = ey[i][j] - 0.5*(hz[i][j]-hz[i-1][j]);
        for (i=0; i<nx; i++)
            for (j=1; j<ny; j++)
                ex[i][j] = ex[i][j] - 0.5*(hz[i][j]-hz[i][j-1]);
        for (i=0; i<nx; i++)
            for (j=0; j<ny; j++)
                hz[i][j]=hz[i][j]-0.7*(ex[i][j+1]-ex[i][j]+ey[i+1][j]-ey[i][j]);
    }
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT;
#endif

    if (fopen(".test", "r")) {
#ifdef MPI
        if (my_rank == 0) {
            print_array();
        }
#else
        print_array();
#endif
    }

    return 0;
}
