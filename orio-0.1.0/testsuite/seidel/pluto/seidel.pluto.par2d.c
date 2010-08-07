
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

double A[N][N+17];

void init_arrays()
{
  int i, j;
  for (i=0; i<N; i++) 
    for (j=0; j<N; j++) 
      A[i][j] = i*i+j*j;
}

double rtclock()
{
  struct timezone tzp;
  struct timeval tp;
  int stat;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main()
{
  init_arrays();

  double annot_t_start=0, annot_t_end=0, annot_t_total=0;
  int annot_i;

  for (annot_i=0; annot_i<REPS; annot_i++)
  {
    annot_t_start = rtclock();
    



#include <math.h>
#include <assert.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

	#define S1(zT0,zT1,zT2,zT3,zT4,zT5,t,i,j)	{A[i][j]=(A[1+i][1+j]+A[1+i][j]+A[1+i][j-1]+A[i][1+j]+A[i][j]+A[i][j-1]+A[i-1][1+j]+A[i-1][j]+A[i-1][j-1])/9;}

	int c1, c2, c3, c4, c5, c6, c7, c8, c9;

	register int lb, ub, lb1, ub1, lb2, ub2;
	register int lbv, ubv;

	omp_set_nested(1);
	omp_set_num_threads(2);
/* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 5.45s. */
for (c1=-2;c1<=floord(4*T+3*N-10,256);c1++) {
	lb1=max(max(max(0,ceild(256*c1-2*T-N-251,512)),ceild(256*c1-3*T-2*N+7,256)),ceild(256*c1-N-761,1024));
	ub1=min(min(min(floord(256*c1+2*N+505,1024),floord(256*c1+509,512)),floord(64*c1+127,64)),floord(T+N-3,256));
#pragma omp parallel for shared(c1,lb1,ub1) private(lb2,ub2,c2,c3,c4,c5,c6,c7,c8,c9)
	for (c2=lb1; c2<=ub1; c2++) {
	lb2=max(max(max(max(max(max(ceild(256*c1-256*c2-T+1,256),ceild(512*c1-512*c2-253,768)),0),ceild(512*c2-N-252,256)),ceild(128*c1-256*c2-127,128)),ceild(128*c2-127,128)),ceild(128*c1-127,256));
	ub2=min(min(min(min(min(min(floord(256*c1-256*c2+255,256),floord(256*c1-512*c2+N+253,256)),floord(256*c2+T+N+252,256)),floord(T+N-3,128)),floord(256*c1+N+508,512)),floord(256*c1-256*c2+N+253,384)),floord(512*c2+N+507,256));
#pragma omp parallel for shared(c1,c2,lb1,ub1,lb2,ub2) private(c3,c4,c5,c6,c7,c8,c9)
	for (c3=lb2; c3<=ub2; c3++) {
      for (c4=max(max(max(max(0,ceild(-256*c2+256*c3-N-284,32)),8*c1-8*c2-8*c3),ceild(256*c2-N-29,32)),ceild(128*c3-N-29,32));c4<=min(min(min(min(8*c1-8*c2-8*c3+7,floord(256*c3+253,64)),floord(T-1,32)),floord(128*c2+127,16)),floord(-128*c2+128*c3+127,16));c4++) {
        for (c5=max(max(max(max(max(8*c2,ceild(16*c4-15,16)),ceild(256*c3-T-N-28,32)),0),ceild(256*c3-32*c4-N-60,32)),ceild(256*c3-N-59,64));c5<=min(min(min(min(min(floord(32*c4+N+29,32),floord(128*c3+127,16)),8*c2+7),floord(128*c3-16*c4+127,16)),floord(T+N-3,32)),floord(256*c3+N+252,64));c5++) {
          for (c6=max(max(max(max(max(ceild(64*c4-29,32),8*c3),ceild(16*c5-15,16)),ceild(16*c4+16*c5-15,16)),0),ceild(64*c5-N-28,32));c6<=min(min(min(min(min(8*c3+7,floord(T+N-3,16)),floord(32*c4+32*c5+N+60,32)),floord(32*c4+N+29,16)),floord(64*c5+N+59,32)),floord(32*c5+T+N+28,32));c6++) {
            for (c7=max(max(max(max(0,32*c4),32*c5-N+2),16*c6-N+2),-32*c5+32*c6-N-29);c7<=min(min(min(min(-32*c5+32*c6+30,floord(32*c6+29,2)),T-1),32*c5+30),32*c4+31);c7++) {
/*@ begin Loop(
transform UnrollJam(ufactor=8)
              for (c8=max(max(32*c5,c7+1),32*c6-c7-N+2);c8<=min(min(32*c6-c7+30,32*c5+31),c7+N-2);c8++) 
transform Unroll(ufactor=8)
                for (c9=max(c7+c8+1,32*c6);c9<=min(32*c6+31,c7+c8+N-2);c9++) 
{
                  S1(c1-c2-c3,-c1+2*c2+c3,-c1+2*c3,c4,-c4+c5,-c4-c5+c6,c7,-c7+c8,-c7-c8+c9) ;
}
) @*/
              for (c8=max(max(32*c5,c7+1),32*c6-c7-N+2);c8<=min(min(32*c6-c7+30,32*c5+31),c7+N-2);c8++) {
                for (c9=max(c7+c8+1,32*c6);c9<=min(32*c6+31,c7+c8+N-2);c9++) {
                  S1(c1-c2-c3,-c1+2*c2+c3,-c1+2*c3,c4,-c4+c5,-c4-c5+c6,c7,-c7+c8,-c7-c8+c9) ;
                }
              }
/*@ end @*/
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */


    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return ((int) A[0][0]); 

}
                                    
