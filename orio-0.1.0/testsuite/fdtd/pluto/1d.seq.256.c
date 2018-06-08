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
    param T1[] = [1];

    param U1[] = [8];
    param U2[] = [1];

    param VEC[] = [False];
  }

  def search
  {
#    arg algorithm = 'Simplex';
    arg algorithm = 'Exhaustive';
#    arg time_limit = 5;
#    arg total_runs = 1;
  }

  def input_params
  {
    param TVAL = 10000;
    param NVAL = 500000;
    decl int T = TVAL;
    decl int N = NVAL;
    decl double coeff1 = 0.5;
    decl double coeff2 = 0.7;
    decl double h[N] = random;
    decl double e[N+1] = random;
  }
) @*/



int t, i, j, k, l,ii;

#define S1(zT0,zT1,t,i)	{e[i]=e[i]-coeff1*(h[i]-h[i-1]);}
#define S2(zT0,zT1,t,i)	{h[i]=h[i]-coeff2*(e[1+i]-e[i]);}

int c1, c2, c3, c4, c5;

register int lbv, ubv;

for (c1=0; c1<=floord(T,256); c1++) {
    for (c2=max(ceild(128*c1-127,128),0); c2<=min(floord(256*c1+N+255,256),floord(N+T,256)); c2++) {
        if ((c1 <= floord(256*c2-N,256)) && (c2 >= ceild(N+1,256))) {
            S2(c1,-c1+c2,256*c2-N,N-1) ;
        }
        for (c3=max(max(256*c2-N+1,256*c1),1); c3<=min(min(T,256*c1+255),256*c2-N+255); c3++) {
            for (c4=max(c3+1,256*c2); c4<=c3+N-1; c4++) {
                S1(c1,-c1+c2,c3,-c3+c4) ;
                S2(c1,-c1+c2,c3,-c3+c4-1) ;
            }
            S2(c1,-c1+c2,c3,N-1) ;
        }
        /*@ begin Loop(
          transform Composite(
            tile = [('c3',T1,'ii')],
            unrolljam = [('c3',U1),('c4',U2)],
            vector = (VEC, ['ivdep','vector always'])
          )
            for (c3=max(max(1,256*c1),256*c2-N+256);c3<=min(min(256*c1+255,T),256*c2+254);c3++)
              for (c4=max(c3+1,256*c2);c4<=256*c2+255;c4++)
        {
                S1(c1,-c1+c2,c3,-c3+c4) ;
                S2(c1,-c1+c2,c3,-c3+c4-1) ;
        }
        ) @*/

        /*@ end @*/
    }
}
/*@ end @*/

