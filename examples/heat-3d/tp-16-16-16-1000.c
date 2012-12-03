#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/*
 * Discretized 3D heat equation stencil with non periodic boundary conditions
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
#define N 200L
#define T 150L
#endif

#define NUM_FP_OPS 15

/* Define our arrays */

double A[2][N+2][N+2][N+2];
double total=0; double sum_err_sqr=0;
int chtotal=0;
/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {
    /* Perform the carry for the later subtraction by updating y. */
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

    /* Compute the time remaining to wait.
     * tv_usec is certainly positive.
     */
    result->tv_sec = x->tv_sec - y->tv_sec;
    result->tv_usec = x->tv_usec - y->tv_usec;

    /* Return 1 if result is negative. */
    return x->tv_sec < y->tv_sec;
}


int main(int argc, char * argv[]) {
    long int t, i, j, k;
    const int BASE = 1024;
    long count=0;
    // for timekeeping
    int ts_return = -1;
    struct timeval start, end, result;
    double tdiff = 0.0;

    printf("Number of points = %ld\t|Number of timesteps = %ld\t", N, T);

    /* Initialization */
    srand(42); // seed with a constant value to verify results

    for (i = 1; i < N+1; i++) {
        for (j = 1; j < N+1; j++) {
            for (k = 1; k < N+1; k++) {
                A[0][i][j][k] = 1.0 * (rand() % BASE);
            }
        }
    }

#ifdef TIME
    gettimeofday(&start, 0);
#endif




  int t1, t2, t3, t4, t5, t6, t7, t8;

 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;

/* Start of CLooG code */
if ((N >= 1) && (T >= 2)) {
  for (t1=0;t1<=floord(2*T+N-4,16);t1++) {
    lbp=max(ceild(t1,2),ceild(16*t1-T+2,16));
    ubp=min(min(floord(T+N-2,16),floord(16*t1+N+15,32)),t1);
#pragma omp parallel for private(lbv,ubv)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(ceild(16*t2-N-14,16),t1-t2);t3<=min(min(floord(T+N-2,16),floord(16*t2+N+14,16)),floord(16*t1-16*t2+N+15,16));t3++) {
        for (t4=max(max(ceild(2*t1-2*t2-124,125),ceild(16*t2-N-998,1000)),ceild(16*t3-N-998,1000));t4<=min(min(min(floord(T+N-2,1000),floord(16*t2+N+14,1000)),floord(16*t3+N+14,1000)),floord(16*t1-16*t2+N+15,1000));t4++) {
          for (t5=max(max(max(16*t1-16*t2,16*t2-N),16*t3-N),1000*t4-N);t5<=min(min(min(min(T-2,16*t2+14),16*t3+14),1000*t4+998),16*t1-16*t2+15);t5++) {
if(t5%2==0){            for (t6=max(16*t2,t5+1);t6<=min(16*t2+15,t5+N);t6++) {
              for (t7=max(16*t3,t5+1);t7<=min(16*t3+15,t5+N);t7++) {
                lbv=max(1000*t4,t5+1);
                ubv=min(1000*t4+999,t5+N);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[1][-t5+t6][-t5+t7][-t5+t8]=0.125*(A[0][-t5+t6+1][-t5+t7][-t5+t8]-2.0*A[0][-t5+t6][-t5+t7][-t5+t8]+A[0][-t5+t6-1][-t5+t7][-t5+t8])+0.125*(A[0][-t5+t6][-t5+t7+1][-t5+t8]-2.0*A[0][-t5+t6][-t5+t7][-t5+t8]+A[0][-t5+t6][-t5+t7-1][-t5+t8])+0.125*(A[0][-t5+t6][-t5+t7][-t5+t8-1]-2.0*A[0][-t5+t6][-t5+t7][-t5+t8]+A[0][-t5+t6][-t5+t7][-t5+t8+1])+A[0][-t5+t6][-t5+t7][-t5+t8];;
                }
              }
            }
}else{
            for (t6=max(16*t2,t5+1);t6<=min(16*t2+15,t5+N);t6++) {
              for (t7=max(16*t3,t5+1);t7<=min(16*t3+15,t5+N);t7++) {
                lbv=max(1000*t4,t5+1);
                ubv=min(1000*t4+999,t5+N);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[0][-t5+t6][-t5+t7][-t5+t8]=0.125*(A[1][-t5+t6+1][-t5+t7][-t5+t8]-2.0*A[1][-t5+t6][-t5+t7][-t5+t8]+A[1][-t5+t6-1][-t5+t7][-t5+t8])+0.125*(A[1][-t5+t6][-t5+t7+1][-t5+t8]-2.0*A[1][-t5+t6][-t5+t7][-t5+t8]+A[1][-t5+t6][-t5+t7-1][-t5+t8])+0.125*(A[1][-t5+t6][-t5+t7][-t5+t8-1]-2.0*A[1][-t5+t6][-t5+t7][-t5+t8]+A[1][-t5+t6][-t5+t7][-t5+t8+1])+A[1][-t5+t6][-t5+t7][-t5+t8];;
                }
              }
            }
}          }
        }
      }
    }
  }
}
/* End of CLooG code */

#ifdef TIME
    gettimeofday(&end, 0);

    ts_return = timeval_subtract(&result, &end, &start);
    tdiff = (double)(result.tv_sec + result.tv_usec * 1.0e-6);

    printf("|Time taken: %7.5lfms\t", tdiff * 1.0e3);
    printf("|MFLOPS: %f\t", ((((double)NUM_FP_OPS * N *N * N * (T-1)) / tdiff) / 1000000L));
#endif

#ifdef VERIFY
    for (i = 1; i < N+1; i++) {
        for (j = 1; j < N+1; j++) {
            for (k = 1; k < N+1; k++) {
                total+= A[T%2][i][j][k] ;
            }
        }
    }
    printf("|sum: %e\t", total);
    for (i = 1; i < N+1; i++) {
        for (j = 1; j < N+1; j++) {
            for (k = 1; k < N+1; k++) {
                sum_err_sqr += (A[T%2][i][j][k] - (total/N))*(A[T%2][i][j][k] - (total/N));
            }
        }
    }
    printf("|rms(A) = %7.2f\t", sqrt(sum_err_sqr));
    for (i = 1; i < N+1; i++) {
        for (j = 1; j < N+1; j++) {
            for (k = 1; k < N+1; k++) {
                chtotal += ((char *)A[T%2][i][j])[k];
            }
        }
    }
    printf("|sum(rep(A)) = %d\n", chtotal);
#endif

    return 0;
}

// icc -O3 -fp-model precise heat_1d_np.c -o op-heat-1d-np -lm
// /* @ begin PrimeTile (num_tiling_levels=1; first_depth=1; last_depth=-1; boundary_tiling_level=-1;) @*/
// /* @ begin PrimeRegTile (scalar_replacement=0; T1t5=4; T1t6=4; T1t7=4; T1t8=4; ) @*/
// /* @ end @*/
