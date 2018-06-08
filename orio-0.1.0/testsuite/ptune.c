/*
 * This is a simple example of how to use the annotations tool to automatically tune
 * a matrix-matrix multiplication kernel. Here only the unroll-and-jam optimization is used.
 */


/*@ begin PerfTuning (
  import spec unrolljam_mm_mul;
) @*/

int i, j, k;

/*@ begin Loop (
transform UnrollJam(ufactor=Ui)
for (i = 0; i <= M-1; i++)
  transform UnrollJam(ufactor=Uj)
  for (j = 0; j <= N-1; j++)
    transform UnrollJam(ufactor=Uk)
    for (k = 0; k <= O-1; k++)
      {
        X[i][j] = X[i][j] + A[i][k] * B[k][j];
      }
) @*/
for (i = 0; i <= M-1; i++)
    for (j = 0; j <= N-1; j++)
        for (k = 0; k <= O-1; k++) {
            X[i][j] = X[i][j] + A[i][k] * B[k][j];
            //X[i*N+j] = X[i*N+j] + A[i*O+k] * B[k*N+j];
        }
/*@ end @*/

/*@ end @*/


