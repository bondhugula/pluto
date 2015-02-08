/*
 * Discretized 1D heat equation stencil with non periodic boundary conditions
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

/*
 * N is the number of points
 * T is the number of timesteps
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 1600000L
#define T 1000L
#endif

#define NUM_FP_OPS 4

/* Define our arrays */
double A[2][N];
double total=0; double sum_err_sqr=0;
long int chtotal=0;
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {
    if (x->tv_usec < y->tv_usec) {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }

    if (x->tv_usec - y->tv_usec > 1000000) {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;

        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    return x->tv_sec < y->tv_sec;
}

int main(int argc, char * argv[]) {
    long int t, i, j, k;
    const int BASE = 1024;

    // for timekeeping
    int ts_return = -1;
    struct timeval start, end, result;
    double tdiff = 0.0;
    long count=0;

    fprintf(stderr, "Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

    /* Initialization */
    srand(42); // seed with a constant value to verify results

    for (i = 0; i < N; i++) {
        A[0][i] = 1.0 * (rand() % BASE);
    }

#ifdef TIME
    gettimeofday(&start, 0);
#endif

    int _N=N-1;
    int _T=T;
#pragma scop
    for (t = 0; t <_T; t++) {
        for (i =0; i <=_N; i++) {
            A[(t+1)%2][i]     =   0.125 * (A[t%2][i==_N?0:i+1]                    - 2.0 * A[t%2][i]     + A[t%2][((i==0)?(_N):(i-1))]);
            /* A[(t+1)%2][_N-1-i] =   0.125 * (A[t%2][((i==0)?(0):(_N-i-1+1))] - 2.0 * A[t%2][_N-1-i] + A[t%2][_N-1-i-1]); */
        }
    }
#pragma endscop

#ifdef TIME
    gettimeofday(&end, 0);

    ts_return = timeval_subtract(&result, &end, &start);
    tdiff = (double)(result.tv_sec + result.tv_usec * 1.0e-6);

    printf("|Time taken =  %7.5lfms\t", tdiff * 1.0e3);
    printf("|MFLOPS =  %f\n", ((((double)NUM_FP_OPS * N *  T) / tdiff) / 1000000L));
#endif

    if (fopen(".test", "r")) {
        total=0;
        for (i = 0; i < N; i++) {
            total+= A[T%2][i] ;
        }
        fprintf(stderr, "|sum: %e\t", total);
        for (i = 0; i < N; i++) {
            sum_err_sqr += (A[T%2][i] - (total/N))*(A[T%2][i] - (total/N));
        }
        fprintf(stderr, "|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
        for (i = 0; i < N; i++) {
            chtotal += ((char *)A[T%2])[i];
        }
        fprintf(stderr, "|sum(rep(A)) = %ld\n", chtotal);
    }
    return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1; boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t3=8; T1t4=8; ) @*/
// /* @ end @*/
// ,t2,t3,t4,t5,t6)
