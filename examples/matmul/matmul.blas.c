#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// #include "mkl.h"

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "mkl.h"
#include "decls.h"

#include "util.h"

double t_start, t_end;

int main()
{
    int i, j, k;
    register double s;
    int LDA, LDB, LDC;
    double *a, *b, *c;
    // char transa='n', transb='n';
    LDA=M;
    LDB=N;
    LDC=M;

    init_array();

    a = (double *) malloc(sizeof(double)*M*M);
    b = (double *) malloc(sizeof(double)*M*M);
    c = (double *) malloc(sizeof(double)*M*M);

    for( i=0; i<M; i++){
        for( j=0; j<N; j++){
            a[i*M+j]= A[i][j]; 
            b[i*M+j]= B[i][j]; 
            c[i*M+j]= C[i][j]; 
        }
    }

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,M,N,K,alpha,a,LDA,b,LDB,beta,c,LDC);

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

#ifdef TEST
    for( i=0; i<M; i++){
        for( j=0; j<N; j++){
            C[i][j] = c[i*M+j];
        }
    }

    print_array();
#endif
    return 0;
}
