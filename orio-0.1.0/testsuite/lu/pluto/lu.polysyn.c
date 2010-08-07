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
    param PERM_A[] = [
      ('i','j'),
      #('j','i'),
      ];

#    param T1_1[] = [8,16,32,64,128];
#    param T1_2[] = [8,16,32,64,128];
#    param T1_3[] = [8,16,32,64,128];
#    param T2_1[] = [1,2,4,8];
#    param T2_2[] = [1,2,4,8];
#    param T2_3[] = [1,2,4,8];

    param T1_1[] = [32];
    param T1_2[] = [16];
    param T1_3[] = [32];
    param T2_1[] = [2];
    param T2_2[] = [8];
    param T2_3[] = [2];

    constraint LU_tiles = ((T1_1 == T1_3) and (T2_1 == T2_3));
#    constraint LU_tiles = ((T2_1 == T2_3));

    param PERM_B[] = [
      [0,1,2],
      #[0,2,1],
      #[2,0,1],
      ];

    param U1[] = [32];
    param U2[] = [8];
    param U3[] = [1];

    param PAR[] = [True];
    param SCREP[] = [True];
    param IVEC[] = [True];
  }

  def search 
  { 
    arg algorithm = 'Exhaustive'; 
#    arg algorithm = 'Simplex'; 
#    arg algorithm = 'Random'; 
    arg time_limit = 5;
#    arg total_runs = 2;
  } 
   
  def input_params 
  {
    param N[] = [256];
    decl in static double L[N][N] = 0; 
    decl in static double U[N][N] = 0; 
    decl out static double A[N][N+13] = random; 
  }
) @*/ 

register int i,j,k;
register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t;
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6,
  newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6,
  newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;


/*@ begin PolySyn(
  l1_tiles = [T1_1,T1_2,T1_3];
  l2_tiles = [T2_1,T2_2,T2_3];
  hotspot_permut = PERM_B;
  unroll_factors = [U1,U2,U3];
  parallelize = PAR;
  scalar_replace = SCREP;
  icc_vectorize = IVEC;
) @*/

/* pluto start (N) */
for (k=0; k<=N-1; k++)
  {
    for (j=k+1; j<=N-1; j++)
      A[k][j] = A[k][j]/A[k][k];
/*@ begin Loop(
    transform Composite(permut = [PERM_A])
    for(i=k+1; i<=N-1; i++)
      for (j=k+1; j<=N-1; j++)
        A[i][j] = A[i][j] - A[i][k]*A[k][j];
) @*/
/*@ end @*/
  }
/* pluto end */
/*@ end @*/
/*@ end @*/

