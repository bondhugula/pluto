#include <stdio.h>
#include <stdlib.h>

#include "mkl.h"

#include "decls.h"
#include "util.h"

double t_start, t_end;

main()
{
    int t, p, q, r, s;
    double *_C4, *_sum, *_A;
    int LDA, LDB, LDC;
    int i, j, k;

    LDA=N;
    LDB=N;
    LDC=N;

    init_array();

    _A = (double *) malloc(sizeof(double)*N*N*N);
    _C4 = (double *) malloc(sizeof(double)*N*N);
    _sum = (double *) malloc(sizeof(double)*N*N*N);

    for( i=0; i<N; i++){
        for( j=0; j<N; j++){
            for( k=0; k<N; k++){
                _A[i*N*N+j*N+k]= A[i][j][k];
            }
            _C4[i*N+j]= C4[i][j]; 
        }
    }

    IF_TIME(t_start = rtclock());

#ifndef TEST
    for (t=0; t<1000; t++)  {
#endif

    /* pluto start (N) */
    for( r = 0; r < N; r++)  {
        for( q = 0; q< N; q++)  {
            for( p = 0; p< N; p++)  {
                _sum[N*N*r + N*q + p] = 0.0;
            }
        }
    }

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,
            N*N,N,N,1,_A,LDA,_C4,LDB,1,_sum,LDC);

    for( r = 0; r<N; r++)  {
        for( q = 0; q<N; q++)  {
            for( p = 0; p< N; p++)  {
                A[r][q][p] = _sum[N*N*r+N*q+p];
            }
        }
    }
    /* pluto end */

#ifndef TEST
    }
#endif

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef TEST
    print_array();
#endif

    return 0;

}
