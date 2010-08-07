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

    param T1_1[] = [32];
    param T1_2[] = [64];
    param T1_3[] = [8];
    param T2_1[] = [1];
    param T2_2[] = [1];
    param T2_3[] = [1];

    param U1[] = [3];
    param U2[] = [1];
    param U3[] = [1];

#    param U1[] = [1,2,4,8];
#    param U2[] = [1,2,4,8];
#    param U3[] = [1,2,4,8];

    param PERM[] = [
#      [0,1,2],
      [0,2,1],
#      [1,0,2],
#      [1,2,0],
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
# arg time_limit = 10;
  arg total_runs = 1;
  }

  def input_params
  {
  param tmax[] = [500];
  param nx[] = [2000];
  param ny[] = [2000]; 
  decl in static double ex[nx][ny+1] = random;
  decl in static double ey[nx+1][ny] = random;
  decl out static double hz[nx][ny] = random;
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
 
  profiling_code = 'fdtd-2d_profiling.c'; 
  compile_cmd = 'gcc'; 
  compile_opts = '-lm'; 
  ) @*/ 



/* pluto start (tmax,nx,ny) */
for(t=0; t<tmax; t++) 
  {
    for (j=0; j<ny; j++)
      ey[0][j] = t;
    for (i=1; i<nx; i++)
      for (j=0; j<ny; j++)
	ey[i][j] = ey[i][j] - 0.5*(hz[i][j]-hz[i-1][j]);
    for (i=0; i<nx; i++)
      for (j=1; j<ny; j++)
	ex[i][j] = ex[i][j] - 0.5*(hz[i][j]-hz[i][j-1]);
    for (i=0; i<nx; i++)
      for (j=0; j<ny; j++)
	hz[i][j]=hz[i][j]-0.7*(ex[i][j+1]-ex[i][j]+ey[i+1][j]-ey[i][j]);
  }
/* pluto end */



/*@ end @*/
/*@ end @*/
