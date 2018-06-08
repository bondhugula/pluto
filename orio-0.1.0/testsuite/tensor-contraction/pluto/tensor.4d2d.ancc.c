
/*@ begin PerfTuning (
  def build
  {
    arg command = 'icc';
    arg options = '-fast -openmp -I/usr/local/icc/include -lm';
  }

  def performance_counter
  {
    arg method = 'basic timer';
    arg repetitions = 5;
  }

  let VR = 1500;
  let OR = 10;

  def performance_params
  {
    param VRANGE = VR;
    param ORANGE = OR;

    param T_O1[] = [1,4,8,16];
    param T_O2[] = [1,4,8,16];
    param T_OX[] = [1,4,8,16];
    param T_V1[] = [1,32,64,128,512,1024];
    param T_V2[] = [1,32,64,128,512,1024];

    param ACOPY_A2[] = [True,False];
    param ACOPY_R[]  = [True,False];
    param ACOPY_T[]  = [True,False];

    param U_O1[] = [1,2,3,4];
    param U_O2[] = [1,2,3,4];
    param U_OX[] = [1,2,3,4];
    param U_V1[] = [1,2,3,4];
    param U_V2[] = [1,2,3,4];

    param SCREP[] = [False,True];
    param VEC[]   = [False,True];
    param OMP[]   = [False,True];

    param PERMUTS[] = [
      [['tv1'],['tv2'],['to1'],['to2'],['tox'],'v1','v2','o1','o2','ox'],
      [['tv1'],['tv2'],['to1'],['to2'],['tox'],'v1','ox','o1','o2','v2'],
    ];

    constraint reg_capacity = (U_V1*U_V2*U_O1*U_O2 +  U_V1*U_OX*U_O1*U_O2 + U_V2*U_OX <= 130);
    constraint unroll_limit = ((U_V1==1) or (U_V2==1) or (U_O1==1) or (U_O2==1) or (U_OX==1));

    constraint copyR = ((not ACOPY_R) or (ACOPY_R and
                         (T_V1 if T_V1>1 else VRANGE)*(T_V2 if T_V2>1 else VRANGE)*
                         (T_O1 if T_O1>1 else ORANGE)*(T_O2 if T_O2>1 else ORANGE) <= 512*512));
    constraint copyT = ((not ACOPY_T) or (ACOPY_T and
                         (T_V1 if T_V1>1 else VRANGE)*(T_OX if T_OX>1 else ORANGE)*
                         (T_O1 if T_O1>1 else ORANGE)*(T_O2 if T_O2>1 else ORANGE) <= 512*512));
    constraint copyA2 = ((not ACOPY_A2) or (ACOPY_A2 and
                         (T_V2 if T_V2>1 else VRANGE)*(T_OX if T_OX>1 else ORANGE) <= 512*512));
    constraint copy_limit =  (not ACOPY_R or not ACOPY_T or not ACOPY_A2);
  }

  def search
  {
    arg algorithm = 'Simplex';
    arg time_limit = 10;
    arg total_runs = 1;
  }

  def input_params
  {
    param VSIZE = VR;
    param OSIZE = OR;
    decl int V = VSIZE;
    decl int O = OSIZE;
    decl double A2[V][O] = random;
    decl double T[V][O][O][O] = random;
    decl double R[V][V][O][O] = 0;
  }
) @*/

int v1,v2,o1,o2,ox;
int tv1,tv2,to1,to2,tox;


/*@ begin Loop(
  transform Composite(
    tile = [('v1',T_V1,'tv1'),('v2',T_V2,'tv2'),('o1',T_O1,'to1'),('o2',T_O2,'to2'),
            ('ox',T_OX,'tox')],
    permut = [PERMUTS],
    arrcopy = [(ACOPY_R,'R[v1][v2][o1][o2]',
            [(T_V1 if T_V1>1 else VRANGE),(T_V2 if T_V2>1 else VRANGE),(T_O1 if T_O1>1 else ORANGE),
             (T_O2 if T_O2>1 else ORANGE)],'_copy'),
           (ACOPY_T,'T[v1][ox][o1][o2]',
            [(T_V1 if T_V1>1 else VRANGE),(T_OX if T_OX>1 else ORANGE),(T_O1 if T_O1>1 else ORANGE),
             (T_O2 if T_O2>1 else ORANGE)],'_copy'),
           (ACOPY_A2,'A2[v2][ox]',
            [(T_V2 if T_V2>1 else VRANGE),(T_OX if T_OX>1 else ORANGE)],'_copy')],
    unrolljam = [('v1',U_V1),('v2',U_V2),('o1',U_O1),('o2',U_O2),('ox',U_OX)],
    scalarreplace = (SCREP, 'double', 'scv_'),
    vector = (VEC, ['ivdep','vector always']),
    openmp = (OMP, 'omp parallel for private(tv1,tv2,to1,to2,tox,v1,v2,o1,o2,ox,R_copy,A2_copy,T_copy)')
  )

  for(v1=0; v1<=V-1; v1++)
    for(v2=0; v2<=V-1; v2++)
      for(o1=0; o1<=O-1; o1++)
        for(o2=0; o2<=O-1; o2++)
	  for(ox=0; ox<=O-1; ox++)
	    R[v1][v2][o1][o2] = R[v1][v2][o1][o2] + T[v1][ox][o1][o2] * A2[v2][ox];

) @*/

for(v1=0; v1<=V-1; v1++)
    for(v2=0; v2<=V-1; v2++)
        for(o1=0; o1<=O-1; o1++)
            for(o2=0; o2<=O-1; o2++)
                for(ox=0; ox<=O-1; ox++)
                    R[v1][v2][o1][o2] = R[v1][v2][o1][o2] + T[v1][ox][o1][o2] * A2[v2][ox];

/*@ end @*/
/*@ end @*/

