
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

  let RANGE = 5000;
  let BSIZE = 512*32;

  def performance_params
  {
    param R = RANGE;

    param T_I[] = [1,256];
    param T_J[] = [1,256];

    param ACOPY_Y[]  = [False,True];
    param ACOPY_X0[] = [False,True];
    param ACOPY_X1[] = [False,True];
    param ACOPY_X2[] = [False,True];

    param U_I[] = [1];
    param U_J[] = [1];

    param SCREP[] = [True];
    param VEC[] = [True];
    param OMP[] = [True];

    param PERMUTS[] = [
		       (['ii'],['jj'],'i','j'),
		       (['jj'],['ii'],'i','j'),
		       (['ii'],['jj'],'j','i'),
		       (['jj'],['ii'],'j','i'),
		       ];

    constraint reg_capacity = (4*U_I*U_J + 3*U_I + 3*U_J <= 128);

    constraint copy_limitY = ((not ACOPY_Y) or (ACOPY_Y and
                              (T_I if T_I>1 else R)*(T_J if T_J>1 else R) <= BSIZE));
    constraint copy_limitX0 = ((not ACOPY_X0) or (ACOPY_X0 and
                               (T_I if T_I>1 else R)*(T_J if T_J>1 else R) <= BSIZE));
    constraint copy_limitX1 = ((not ACOPY_X1) or (ACOPY_X1 and
                               (T_I if T_I>1 else R)*(T_J if T_J>1 else R) <= BSIZE));
    constraint copy_limitX2 = ((not ACOPY_X2) or (ACOPY_X2 and
                               (T_I if T_I>1 else R)*(T_J if T_J>1 else R) <= BSIZE));
    constraint copy_limit = (not ACOPY_Y or not ACOPY_X0 or not ACOPY_X1 or not ACOPY_X2);
  }

  def search
  {
    arg algorithm = 'Simplex';
    arg time_limit = 5;
    arg total_runs = 2;
  }

  def input_params
  {
    param SIZE = RANGE;
    decl int N = SIZE;
    decl double X0[N][N] = random;
    decl double X1[N][N] = random;
    decl double X2[N][N] = random;
    decl double Y[N][N] = 0;
    decl double u0[N] = random;
    decl double u1[N] = random;
    decl double u2[N] = random;
    decl double a0 = 32.12;
    decl double a1 = 3322.12;
    decl double a2 = 1.123;
    decl double b00 = 1321.9;
    decl double b01 = 21.55;
    decl double b02 = 10.3;
    decl double b11 = 1210.313;
    decl double b12 = 9.373;
    decl double b22 = 1992.31221;
  }
) @*/

int i,j,ii,jj,iii,jjj;

/*@ begin Loop(
  transform Composite(
    tile = [('i',T_I,'ii'),('j',T_J,'jj')],
    permut = [PERMUTS],
    arrcopy = [(ACOPY_Y,'Y[i][j]',[(T_I if T_I>1 else R),(T_J if T_J>1 else R)],'_copy'),
               (ACOPY_X0,'X0[i][j]',[(T_I if T_I>1 else R),(T_J if T_J>1 else R)],'_copy'),
               (ACOPY_X1,'X1[i][j]',[(T_I if T_I>1 else R),(T_J if T_J>1 else R)],'_copy'),
               (ACOPY_X2,'X2[i][j]',[(T_I if T_I>1 else R),(T_J if T_J>1 else R)],'_copy')],
    unrolljam = [('j',U_J),('i',U_I)],
    scalarreplace = (SCREP, 'double', 'scv_'),
    vector = (VEC, ['ivdep','vector always']),
    openmp = (OMP, 'omp parallel for private(i,j,ii,jj,iii,jjj,Y_copy,X0_copy,X1_copy,X2_copy)')
  )

for (i=0; i<=N-1; i++)
  for (j=0; j<=N-1; j++)
    {
      Y[i][j]=a0*X0[i][j] + a1*X1[i][j] + a2*X2[i][j]
	+ 2.0*b00*u0[i]*u0[j]
	+ 2.0*b11*u1[i]*u1[j]
	+ 2.0*b22*u2[i]*u2[j]
	+ b01*u0[i]*u1[j] + b01*u1[i]*u0[j]
	+ b02*u0[i]*u2[j] + b02*u2[i]*u0[j]
	+ b12*u1[i]*u2[j] + b12*u2[i]*u1[j];
    }

) @*/

for (i=0; i<=N-1; i++)
    for (j=0; j<=N-1; j++) {
        Y[i][j]=a0*X0[i][j] + a1*X1[i][j] + a2*X2[i][j]
                + 2.0*b00*u0[i]*u0[j]
                + 2.0*b11*u1[i]*u1[j]
                + 2.0*b22*u2[i]*u2[j]
                + b01*(u0[i]*u1[j] + u1[i]*u0[j])
                + b02*(u0[i]*u2[j] + u2[i]*u0[j])
                + b12*(u1[i]*u2[j] + u2[i]*u1[j]);
    }

/*@ end @*/
/*@ end @*/
