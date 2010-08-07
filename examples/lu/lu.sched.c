#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "decls.h"

#include "omp.h"

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "util.h"

int main()
{
	int i, j, k;
    double t_start, t_end;

	init_array() ;

	IF_TIME(t_start = rtclock());

#include <math.h>
#include <assert.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

	#define S1(zT0,zT1,k,j)	{a[k][j]=a[k][j]/a[k][k];}
	#define S2(zT0,zT1,zT2,k,i,j)	{a[i][j]=a[i][j]-a[i][k]*a[k][j];}

	int c1, c2, c3, c4, c5, c6;

	register int lbv, ubv;

/* Generated from lu.sched.cloog by CLooG v0.14.1 64 bits in 0.05s. */
for (c1=-1;c1<=floord(N-2,8);c1++) {
  for (c2=max(ceild(4*c1-149,150),0);c2<=floord(N-1,300);c2++) {
    for (c3=max(ceild(4*c1-105,120),ceild(4*c1-7,8));c3<=floord(N-1,16);c3++) {
      if ((c1 == 2*c3) && (c3 <= 0)) {
        for (c5=max(1,300*c2);c5<=min(300*c2+299,N-1);c5++) {
          if (c1%2 == 0) {
            S1(c1/2,c2,0,c5) ;
          }
        }
      }
      if ((-c1 == -2*c3) && (c1 >= 1)) {
        for (c5=max(8*c1+1,300*c2);c5<=min(300*c2+299,N-1);c5++) {
          if (c1%2 == 0) {
            S1(c1/2,c2,8*c1,c5) ;
          }
        }
      }
      if (c1 <= 2*c3-1) {
        for (c4=max(16*c1+1,1);c4<=min(min(600*c2+597,2*N-3),16*c1+31);c4++) {
          for (c5=max(ceild(c4+1,2),300*c2);c5<=min(N-1,300*c2+299);c5++) {
            for (c6=max(ceild(c4+1,2),16*c3);c6<=min(N-1,16*c3+15);c6++) {
              if (c1%2 == 0) {
                if ((c4-1)%2 == 0) {
                  S2(c1/2,c3,c2,(c4-1)/2,c6,c5) ;
                }
              }
            }
          }
        }
      }
      if (c1 >= 2*c3+1) {
        for (c4=max(1,16*c1+1);c4<=min(min(32*c3+29,600*c2+597),2*N-3);c4++) {
          for (c5=max(300*c2,ceild(c4+1,2));c5<=min(N-1,300*c2+299);c5++) {
            for (c6=ceild(c4+1,2);c6<=min(16*c3+15,N-1);c6++) {
              if (c1%2 == 0) {
                if ((c4-1)%2 == 0) {
                  S2(c1/2,c3,c2,(c4-1)/2,c6,c5) ;
                }
              }
            }
          }
        }
      }
      if (c1 == 2*c3) {
        for (c4=max(1,32*c3+1);c4<=min(32*c3+29,600*c2-2);c4++) {
          for (c5=300*c2;c5<=min(N-1,300*c2+299);c5++) {
            if (c4%2 == 0) {
              if (c1%2 == 0) {
                if (c4%2 == 0) {
                  S1(c1/2,c2,c4/2,c5) ;
                }
              }
            }
            for (c6=ceild(c4+1,2);c6<=min(16*c3+15,N-1);c6++) {
              if (c1%2 == 0) {
                if (c1%2 == 0) {
                  if ((c4-1)%2 == 0) {
                    S2(c1/2,c1/2,c2,(c4-1)/2,c6,c5) ;
                  }
                }
              }
            }
          }
        }
      }
      if (c1 == 2*c3) {
          if ((c4+1)%2 == 0) {
            for (c6=ceild(c4+1,2);c6<=min(16*c3+15,N-1);c6++) {
        for (c4=max(max(600*c2-1,1),32*c3+1);c4<=min(min(32*c3+29,600*c2+596),2*N-4);c4++) {
              if (c1%2 == 0) {
                if (c1%2 == 0) {
                  if ((c4-1)%2 == 0) {
                    if ((c4+1)%2 == 0) {
                      S2(c1/2,c1/2,c2,(c4-1)/2,c6,(c4+1)/2) ;
                    }
                  }
                }
              }
            }
          }
          for (c5=ceild(c4+2,2);c5<=min(N-1,300*c2+299);c5++) {
            if (c4%2 == 0) {
              if (c1%2 == 0) {
                if (c4%2 == 0) {
                  S1(c1/2,c2,c4/2,c5) ;
                }
              }
            }
            for (c6=ceild(c4+1,2);c6<=min(16*c3+15,N-1);c6++) {
              if (c1%2 == 0) {
                if (c1%2 == 0) {
                  if ((c4-1)%2 == 0) {
                    S2(c1/2,c1/2,c2,(c4-1)/2,c6,c5) ;
                  }
                }
              }
            }
          }
        }
      }
      if ((-c1 == -2*c3) && (c1 >= -1) && (c1 <= min(floord(N-17,8),floord(300*c2+283,8)))) {
        for (c5=max(8*c1+16,300*c2);c5<=min(300*c2+299,N-1);c5++) {
          if (c1%2 == 0) {
            S1(c1/2,c2,8*c1+15,c5) ;
          }
        }
      }
      if ((c1 == 2*c3) && (c2 >= ceild(N-300,300)) && (c3 >= ceild(N-16,16))) {
        if (c1%2 == 0) {
          if (c1%2 == 0) {
            S2(c1/2,c1/2,c2,N-2,N-1,N-1) ;
          }
        }
      }
      if ((c1 == 2*c3) && (c2 <= min(floord(4*c3-71,75),floord(2*N-601,600)))) {
        for (c6=300*c2+299;c6<=min(16*c3+15,N-1);c6++) {
          if (c1%2 == 0) {
            if (c1%2 == 0) {
              S2(c1/2,c1/2,c2,300*c2+298,c6,300*c2+299) ;
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */

	IF_TIME(t_end = rtclock());
	IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));


#ifdef TEST
    print_array();
#endif
    return 0;
}
