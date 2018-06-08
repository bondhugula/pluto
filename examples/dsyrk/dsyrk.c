#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <unistd.h>
#include <sys/time.h>

#define NMAX 2000
#pragma declarations
double a[NMAX][NMAX], c[NMAX][NMAX];
#pragma enddeclarations

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

#ifdef TIME
double rtclock() {
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
#endif


int main() {
    double t_start, t_end;

    int i,j, k;

    for (i = 0; i < NMAX; i++) {
        for (j = 0; j < NMAX; j++) {
            c[i][j] = 0.0;
            a[i][j] = i*j*3.2345 / NMAX;
        }
    }

    IF_TIME(t_start = rtclock());

#pragma scop
    for (i=0; i<NMAX; i++) {
        for (j=0; j<NMAX; j++) {
            for (k=j; k<NMAX; k++) {
                c[j][k] += a[i][j] * a[i][k];
            }
        }
    }
#pragma endscop


    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));

    if (fopen(".test", "r"))  {
#ifdef MPI
        if (my_rank == 0) {
#endif
            for (i = 0; i < NMAX; i++) {
                for (j = 0; j < NMAX; j++) {
                    fprintf(stderr, "%lf ", c[i][j]);
                }
            }
        }
#ifdef MPI
    }
#endif

    return 0;

}
