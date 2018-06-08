/*@ begin PerfTuning (
  def build
  {
  arg command = 'icc';
  arg options = '-fast -openmp -I/usr/local/icc/include -lm';
  }

  def performance_counter
  {
  arg method = 'basic timer';
  arg repetitions = 4;
  }

  def performance_params
  {
    param T1_1[] = [512];
    param T1_2[] = [1];
    param T2_1[] = [8];
    param T2_2[] = [4];

    param U1[] = [1];
    param U2[] = [4];
    param U3[] = [4];
    param U4[] = [1];

    param PERM1[] = [
      [0,1],
#      [1,0],
    ];
    param PERM2[] = [
     ['i','j'],
#     ['j','i'],
    ];

    param PAR[] = [True];
    param SCREP1[] = [True];
    param SCREP2[] = [True];
    param IVEC1[] = [True];
    param IVEC2[] = [True];
  }

  def search
  {
  arg algorithm = 'Exhaustive';
# arg algorithm = 'Simplex';
# arg algorithm = 'Random';
# arg time_limit = 10;
  arg total_runs = 1;
  }

  def input_params
  {
  param N[] = [10000];
  param alpha = 1;
  param beta = 1;

  decl in static double A[N][N+20] = 0;
  decl in static double B[N][N+20] = 0;
  decl in static double x[N] = 0;
  decl in static double u1[N] = 0;
  decl in static double u2[N] = 0;
  decl in static double v1[N] = 0;
  decl in static double v2[N] = 0;
  decl out static double w[N] = 0;
  decl in static double y[N] = 0;
  decl in static double z[N] = 0;
  }
) @*/

register int i,j,k,t;
register int it,jt,kt,tt;
register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t;
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6,
         newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6,
         newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;

/*@ begin PolySyn(
  parallel = PAR;
  tiles = [T1_1,T1_2,T2_1,T2_2];
  permut = PERM1;
  unroll_factors = [U1,U2];
  scalar_replace = SCREP1;
  vectorize = IVEC1;

  profiling_code = 'gemver_profiling.c';
  compile_cmd = 'gcc';
  compile_opts = '-lm';
  ) @*/


/* pluto start (N,alpha,beta) */
for (i=0; i<=N-1; i++)
    for (j=0; j<=N-1; j++)
        B[i][j] = A[i][j] + u1[i]*v1[j] + u2[i]*v2[j];
for (i=0; i<=N-1; i++)
    for (j=0; j<=N-1; j++)
        x[i] = x[i] + beta* B[j][i]*y[j];
/* pluto end */

/*@ end @*/

for (i=0; i<=N-1; i++)
    x[i] = x[i] + z[i];

/*@ begin Loop(
 transform Composite(
  permut = [PERM2],
  regtile = [['i','j'],[U3,U4]],
  scalarreplace = (SCREP2, 'double', 'itv_'),
  vector = (IVEC2, ['ivdep','vector always'])
 )
 for (i=0; i<=N-1; i++)
   for (j=0; j<=N-1; j++)
     w[i] = w[i] + alpha* B[i][j]*x[j];
     ) @*/
/*@ end @*/

/*@ end @*/
