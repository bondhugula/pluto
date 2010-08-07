#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#include "decls.h"
#include "util.h"

int main()
{
    int i, j, k, t;

    double t_start, t_end;

    init_array() ;

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

#include <math.h>
#include <assert.h>

#include <omp.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

	#define S1(zT0,zT1,zT2,t,i,j)	{a[i][j]=(a[1+i][1+j]+a[1+i][j]+a[1+i][j-1]+a[i][1+j]+a[i][j]+a[i][j-1]+a[i-1][1+j]+a[i-1][j]+a[i-1][j-1])/9;}

	int c1, c2, c3, c4, c5, c6;

	register int lb, ub, lb1, ub1, lb2, ub2;
	register int lbv, ubv;

/* Generated from seidel.sched.cloog by CLooG v0.14.1 64 bits in 0.03s. */
for (c1=0;c1<=floord(3*N+4*T-10,50);c1++) {
  for (c2=max(0,ceild(50*c1-3*N-190,200));c2<=min(floord(T-1,50),floord(25*c1+23,100));c2++) {
    for (c3=max(max(max(max(ceild(25*c2-24,25),ceild(50*c1-N-2*T-94,100)),0),ceild(50*c1-100*c2-N-194,100)),ceild(50*c1-N-192,200));c3<=min(min(min(min(floord(50*c2+N+47,50),floord(N+T-3,50)),floord(25*c1+N+22,100)),floord(25*c1+24,50)),floord(25*c1-50*c2+24,50));c3++) {
      for (c4=max(max(max(max(max(100*c2+100*c3+1,200*c3-2*N+5),50*c1),100*c3+1),200*c2+3),3);c4<=min(min(min(min(min(50*c1+49,200*c2+3*N+190),200*c3+N+192),100*c3+N+2*T+94),3*N+4*T-10),100*c2+100*c3+N+194);c4++) {
        for (c5=max(max(max(max(0,50*c3-N+2),ceild(c4-3*N+6,4)),ceild(-100*c3+c4-N-96,2)),50*c2);c5<=min(min(min(min(50*c2+49,T-1),floord(-100*c3+c4-1,2)),50*c3+48),floord(c4-3,4));c5++) {
          for (c6=max(max(ceild(c4-2*c5-N+2,2),50*c3),c5+1);c6<=min(min(c5+N-2,floord(c4-2*c5-1,2)),50*c3+49);c6++) {
            S1(c2,-c2+c3,c1-2*c2-2*c3,c5,-c5+c6,c4-2*c5-2*c6) ;
          }
        }
      }
    }
  }
}

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

#ifdef TEST
    print_array();
#endif
    return 0;
}
