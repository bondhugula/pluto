#include <stdio.h>
#include <stdlib.h>
#include <acml.h>

#include "decls.h"
#include "util.h"

main()
{
    double t_start, t_end;
    int i, j;
    int M=N;
    double *a, *b;
    int n;

    init_array();

    a = (double *) malloc(sizeof(double)*M*M);
    b = (double *) malloc(sizeof(double)*M*M);

    for( i=0; i<M; i++){
        for( j=0; j<N; j++){
            a[i*M+j]= A[j][i]; 
        }
    }

    n = M*N;

    IF_TIME(t_start = rtclock());

    dcopy(M*N, a, 1, b, 1);
    dger(M, N, 1.0, u1, 1, v1, 1, b, N);
    dger(M, N, 1.0, u2, 1, v2, 1, b, N);
    dcopy(N, z, 1, x, 1);
    dgemv('T',M,N,beta,b,M,y,1,1.0,x,1);
    dgemv('N',M,N,alpha,b,M,x,1,1.0,w,1);

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef TEST
    print_array();
#endif

    return 0;
}

