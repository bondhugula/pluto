/*@ begin PerfTuning (
  def build 
  {
    arg command = 'icc';  
    arg options = '-fast  -I/usr/local/icc/include -lm';   
  }
  
  def performance_counter 
  {
    arg method = 'basic timer';
    arg repetitions = 1;  
  }
   
  def performance_params  
  {
    param T1[] = [1];
    param T2[] = [1];
    
    param PERMUTS[] = [   
      #('c5','c6'),
      ('c6','c5'),
      ];  
  
      param U1[] = [1];
      param U2[] = [4];
      
      param VEC[] = [False]; 
  }
  
  def search
  {
    #arg algorithm = 'Simplex';
    arg algorithm = 'Exhaustive';  
    #arg time_limit = 5;  
    #arg total_runs = 1;  
  }

  def input_params   
  {
    param TVAL = 500;
    param NXVAL = 500;   
    param NYVAL = 500;   
    decl int tmax = TVAL; 
    decl int nx = NXVAL;  
    decl int ny = NYVAL;  
    decl double ex[nx][ny+1] = random; 
    decl double ey[nx+1][ny] = random; 
    decl double hz[nx][ny] = random;   
  }
  ) @*/


int t, i, j, k, l, m, n,ii,jj;

	#define S1(zT0,zT1,t,j)	{ey[0][j]=t;}
	#define S2(zT0,zT1,zT2,t,i,j)	{ey[i][j]=ey[i][j]-((double)(1))/2*(hz[i][j]-hz[i-1][j]);}
	#define S3(zT0,zT1,zT2,t,i,j)	{ex[i][j]=ex[i][j]-((double)(1))/2*(hz[i][j]-hz[i][j-1]);}
	#define S4(zT0,zT1,zT2,t,i,j)	{hz[i][j]=hz[i][j]-((double)(7))/10*(ey[1+i][j]+ex[i][1+j]-ex[i][j]-ey[i][j]);}

	int c1, c2, c3, c4, c5, c6, c7;

	register int lbv, ubv;

for (c1=0;c1<=floord(tmax-1,32);c1++) {
  for (c2=max(ceild(32*c1-31,32),0);c2<=min(floord(tmax+ny-1,32),floord(32*c1+ny+31,32));c2++) {
for (c3=max(max(max(max(ceild(32*c2-ny-30,32),0),ceild(64*c1-32*c2-61,32)),ceild(32*c1-31,32)),ceild(32*c1-992*c2-1891,992));c3<=min(min(floord(32*c2+nx+30,32),floord(tmax+nx-1,32)),floord(32*c1+nx+31,32));c3++) {
      if ((c1 <= floord(32*c3-nx,32)) && (c2 <= floord(32*c3-nx+ny,32)) && (c3 >= ceild(nx,32))) {
        for (c5=max(32*c3-nx+1,32*c2);c5<=min(32*c2+31,32*c3-nx+ny);c5++) {
          S4(c1,-c1+c3,-c1+c2,32*c3-nx,nx-1,-32*c3+c5+nx-1) ;
        }
      }
      if ((c1 <= floord(32*c2-ny,32)) && (c2 >= max(ceild(32*c3-nx+ny+1,32),ceild(ny,32)))) {
        for (c6=max(32*c3,32*c2-ny+1);c6<=min(32*c2+nx-ny,32*c3+31);c6++) {
          S4(c1,-c1+c3,-c1+c2,32*c2-ny,-32*c2+c6+ny-1,ny-1) ;
        }
      }
      if (c1 == c3) {
        for (c4=max(max(32*c2-ny+1,0),32*c3);c4<=min(min(32*c3+30,32*c2-ny+31),tmax-1);c4++) {
          for (c5=32*c2;c5<=c4+ny-1;c5++) {
            S1(c1,-c1+c2,c4,-c4+c5) ;
            S3(c1,0,-c1+c2,c4,0,-c4+c5) ;
            for (c6=c4+1;c6<=32*c3+31;c6++) {
              S2(c1,0,-c1+c2,c4,-c4+c6,-c4+c5) ;
              S3(c1,0,-c1+c2,c4,-c4+c6,-c4+c5) ;
              S4(c1,0,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
            }
          }
          for (c6=c4+1;c6<=32*c3+31;c6++) {
            S4(c1,0,-c1+c2,c4,-c4+c6-1,ny-1) ;
          }
        }
      }
      if (c1 == c3) {
        for (c4=max(max(0,32*c3),32*c2-ny+32);c4<=min(min(tmax-1,32*c3+30),32*c2-1);c4++) {
          for (c5=32*c2;c5<=32*c2+31;c5++) {
            S1(c1,-c1+c2,c4,-c4+c5) ;
            S3(c1,0,-c1+c2,c4,0,-c4+c5) ;
            for (c6=c4+1;c6<=32*c3+31;c6++) {
              S2(c1,0,-c1+c2,c4,-c4+c6,-c4+c5) ;
              S3(c1,0,-c1+c2,c4,-c4+c6,-c4+c5) ;
              S4(c1,0,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
            }
          }
        }
      }
      if (c1 == c3) {
        for (c4=max(max(32*c2,0),32*c3);c4<=min(min(tmax-1,32*c3+30),32*c2+30);c4++) {
          S1(c1,-c1+c2,c4,0) ;
          for (c6=c4+1;c6<=32*c3+31;c6++) {
            S2(c1,0,-c1+c2,c4,-c4+c6,0) ;
          }
          for (c5=c4+1;c5<=32*c2+31;c5++) {
            S1(c1,-c1+c2,c4,-c4+c5) ;
            S3(c1,0,-c1+c2,c4,0,-c4+c5) ;
            for (c6=c4+1;c6<=32*c3+31;c6++) {
              S2(c1,0,-c1+c2,c4,-c4+c6,-c4+c5) ;
              S3(c1,0,-c1+c2,c4,-c4+c6,-c4+c5) ;
              S4(c1,0,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
            }
          }
        }
      }
      for (c4=max(max(max(32*c1,0),32*c2-ny+1),32*c3-nx+1);c4<=min(min(min(32*c3-nx+31,32*c2-ny+31),32*c1+31),tmax-1);c4++) {
        for (c5=32*c2;c5<=c4+ny-1;c5++) {
          for (c6=32*c3;c6<=c4+nx-1;c6++) {
            S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S3(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
          }
          S4(c1,-c1+c3,-c1+c2,c4,nx-1,-c4+c5-1) ;
        }
        for (c6=32*c3;c6<=c4+nx;c6++) {
          S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,ny-1) ;
        }
      }
      for (c4=max(max(max(32*c1,0),32*c3-nx+1),32*c2-ny+32);c4<=min(min(min(tmax-1,32*c1+31),32*c2-1),32*c3-nx+31);c4++) {
        for (c5=32*c2;c5<=32*c2+31;c5++) {
          for (c6=32*c3;c6<=c4+nx-1;c6++) {
            S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S3(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
          }
          S4(c1,-c1+c3,-c1+c2,c4,nx-1,-c4+c5-1) ;
        }
      }
      for (c4=max(max(max(32*c3-nx+32,32*c1),0),32*c2-ny+1);c4<=min(min(min(32*c2-ny+31,32*c1+31),tmax-1),32*c3-1);c4++) {
        for (c5=32*c2;c5<=c4+ny-1;c5++) {
          for (c6=32*c3;c6<=32*c3+31;c6++) {
            S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S3(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
          }
        }
        for (c6=32*c3;c6<=32*c3+31;c6++) {
          S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,ny-1) ;
        }
      }
      for (c4=max(max(max(32*c2,32*c1),0),32*c3-nx+1);c4<=min(min(min(32*c2+30,tmax-1),32*c1+31),32*c3-nx+31);c4++) {
        for (c6=32*c3;c6<=c4+nx-1;c6++) {
          S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,0) ;
        }
        for (c5=c4+1;c5<=32*c2+31;c5++) {
          for (c6=32*c3;c6<=c4+nx-1;c6++) {
            S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S3(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
          }
          S4(c1,-c1+c3,-c1+c2,c4,nx-1,-c4+c5-1) ;
        }
      }
      for (c4=max(max(max(32*c1,0),32*c3-nx+32),32*c2-ny+32);c4<=min(min(min(32*c3-1,tmax-1),32*c1+31),32*c2-1);c4++) {
/*@ begin Loop(
 transform Composite(                                                                        
  tile = [('c5',T1,'ii'),('c6',T2,'jj')],
  permut = [PERMUTS],
  unrolljam = [('c5',U1),('c6',U2)],
  vector = (VEC, ['ivdep','vector always'])
 )                        
        for (c5=32*c2;c5<=32*c2+31;c5++) 
          for (c6=32*c3;c6<=32*c3+31;c6++) 
{
            S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S3(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
}
) @*/

/*@ end @*/
      }
      for (c4=max(max(max(32*c2,32*c3-nx+32),32*c1),0);c4<=min(min(min(32*c3-1,32*c2+30),tmax-1),32*c1+31);c4++) {
        for (c6=32*c3;c6<=32*c3+31;c6++) {
          S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,0) ;
        }
        for (c5=c4+1;c5<=32*c2+31;c5++) {
          for (c6=32*c3;c6<=32*c3+31;c6++) {
            S2(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S3(c1,-c1+c3,-c1+c2,c4,-c4+c6,-c4+c5) ;
            S4(c1,-c1+c3,-c1+c2,c4,-c4+c6-1,-c4+c5-1) ;
          }
        }
      }
      if ((c1 == c3) && (c2 <= min(floord(32*c3-1,32),floord(tmax-32,32)))) {
        S1(c1,-c1+c2,32*c2+31,0) ;
        for (c6=32*c2+32;c6<=32*c3+31;c6++) {
          S2(c1,0,-c1+c2,32*c2+31,-32*c2+c6-31,0) ;
        }
      }
      if ((-c1 == -c3) && (c1 >= ceild(32*c2-31,32)) && (c1 <= min(floord(tmax-32,32),floord(32*c2-1,32)))) {
        S1(c1,-c1+c2,32*c1+31,0) ;
        for (c5=32*c1+32;c5<=32*c2+31;c5++) {
          S1(c1,-c1+c2,32*c1+31,-32*c1+c5-31) ;
          S3(c1,0,-c1+c2,32*c1+31,0,-32*c1+c5-31) ;
        }
      }
      if ((-c1 == -c3) && (c1 <= min(floord(tmax-32,32),c2-1))) {
        for (c5=32*c2;c5<=min(32*c2+31,32*c1+ny+30);c5++) {
          S1(c1,-c1+c2,32*c1+31,-32*c1+c5-31) ;
          S3(c1,0,-c1+c2,32*c1+31,0,-32*c1+c5-31) ;
        }
      }
      if ((-c1 == -c2) && (-c1 == -c3) && (c1 <= floord(tmax-32,32))) {
        S1(c1,0,32*c1+31,0) ;
      }
      if ((c1 >= c2) && (c2 <= min(c3-1,floord(tmax-32,32)))) {
        for (c6=32*c3;c6<=min(32*c2+nx+30,32*c3+31);c6++) {
          S2(c1,-c1+c3,-c1+c2,32*c2+31,-32*c2+c6-31,0) ;
        }
      }
    }
  }
}

/*@ end @*/
