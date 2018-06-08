
/*@ begin PolySyn(
 parallel = True;
 tiles = [64,64,64,1,1,1];
 unroll_factors = [2,2,2];
 scalar_replace = True;
 vectorize = True;

 skeleton_code = 'polysyn_profiling.c';
 compile_cmd = 'gcc';
 compile_opts = '-lm -openmp';
) @*/

int i,j,k;

/* pluto start (N) */
for (k=0; k<=N-1; k++) {
    for (j=k+1; j<=N-1; j++)
        A[k][j] = A[k][j] / A[k][k];
    for(i=k+1; i<=N-1; i++)
        for (j=k+1; j<=N-1; j++)
            A[i][j] = A[i][j] - A[i][k]*A[k][j];
}
/* pluto end */

/*@ end @*/


