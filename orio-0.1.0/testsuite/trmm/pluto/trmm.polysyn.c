/*@ begin PerfTuning (
  def build
  {
  arg command = 'icc';
  arg options = '-fast -openmp -I/usr/local/icc/include -lm';
  }

  def performance_counter
  {
  arg method = 'basic timer';
  arg repetitions = 1;
  }

  def performance_params
  {
#    param T1_1[] = [1,32,64,128,256,512];
#    param T1_2[] = [1,32,64,128,256,512];
#    param T1_3[] = [1,32,64,128,256,512];
#    param T2_1[] = [1,4,8,16,32];
#    param T2_2[] = [1,4,8,16,32];
#    param T2_3[] = [1,4,8,16,32];

    param T1_1[] = [128];
    param T1_2[] = [128];
    param T1_3[] = [256];
    param T2_1[] = [8];
    param T2_2[] = [4];
    param T2_3[] = [1];

    param U1[] = [1];
    param U2[] = [11];
    param U3[] = [8];

#    param U1[] = [1,2,4,8];
#    param U2[] = [1,2,4,8];
#    param U3[] = [1,2,4,8];

    param PERM[] = [
#      [0,1,2],
#      [1,0,2],
      [1,2,0],

#      [0,2,1],
#      [2,0,1],
#      [2,1,0],
    ];

    param PAR[] = [True];
    param SCREP[] = [True];
    param IVEC[] = [True];
  }

  def search
  {
  arg algorithm = 'Exhaustive';
# arg algorithm = 'Simplex';
# arg algorithm = 'Random';
#  arg time_limit = 10;
 arg total_runs = 1;
  }

  def input_params
  {
  param N[] = [4096];
  param alpha = 1;
  decl in static double A[N][N+20] = random;
  decl out static double B[N][N+20] = random;
  }
) @*/

register int i,j,k,t;
register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t;
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6,
         newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6,
         newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;

/*@ begin PolySyn(
  parallel = PAR;
  tiles = [T1_1,T1_2,T1_3,T2_1,T2_2,T2_3];
  permut = PERM;
  unroll_factors = [U1,U2,U3];
  scalar_replace = SCREP;
  vectorize = IVEC;

  profiling_code = 'trmm_profiling.c';
  compile_cmd = 'gcc';
  compile_opts = '-lm';
  ) @*/

/* pluto start (N,alpha) */
for (i=0; i<=N-1; i++)
    for (j=0; j<=N-1; j++)
        for (k=i+1; k<=N-1; k++)
            B[i][j] = B[i][j] + alpha*A[i][k]*B[k][j];
/* pluto end */

/*@ end @*/
/*@ end @*/
