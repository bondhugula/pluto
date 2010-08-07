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
    param T1[] = [1,8,16,32,64];
   
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
    param NVAL = 1000000;     
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

	register int lb, ub, lb1, ub1, lb2, ub2;
	register int lbv, ubv;

for (c1=-1;c1<=floord(N+2*T,512);c1++) {
	lb1=max(max(0,ceild(256*c1-255,512)),ceild(512*c1-T,512));
	ub1=min(min(floord(256*c1+255,256),floord(512*c1+N+511,1024)),floord(N+T,512));
#pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5)
	for (c2=lb1; c2<=ub1; c2++) {
    if ((c1 <= floord(1024*c2-N,512)) && (c2 >= ceild(N+1,512))) {
      S2(c1-c2,-c1+2*c2,512*c2-N,N-1) ;
    }
    for (c3=max(max(512*c2-N+1,512*c1-512*c2),1);c3<=min(min(512*c1-512*c2+511,512*c2-N+511),T);c3++) {
      for (c4=max(512*c2,c3+1);c4<=c3+N-1;c4++) {
        S1(c1-c2,-c1+2*c2,c3,-c3+c4) ;
        S2(c1-c2,-c1+2*c2,c3,-c3+c4-1) ;
      }
      S2(c1-c2,-c1+2*c2,c3,N-1) ;
    }
/*@ begin Loop(
  transform Composite(
   tile = [('c3',T1,'ii')],
   unrolljam = [('c3',U1),('c4',U2)],
   vector = (VEC, ['ivdep','vector always'])
  )              
    for (c3=max(max(1,512*c1-512*c2),512*c2-N+512);c3<=min(min(512*c1-512*c2+511,T),512*c2+510);c3++) 
      for (c4=max(512*c2,c3+1);c4<=512*c2+511;c4++) 
{
        S1(c1-c2,-c1+2*c2,c3,-c3+c4) ;
        S2(c1-c2,-c1+2*c2,c3,-c3+c4-1) ;
}
) @*/

/*@ end @*/
  }
}
/*@ end @*/
