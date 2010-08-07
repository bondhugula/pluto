#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "decls.h"

#include "util.h"

double t_start, t_end;

int main()
{
	int t, i, j, k, l, m, n;


	init_array() ;

#ifdef PERFCTR
    PERF_INIT;
#endif

	IF_TIME(t_start = rtclock());






  int t1, t2, t3, t4, t5, t6;

 register int lbv, ubv;

/* Generated from PLUTO-produced CLooG file by CLooG 0.14.0-201-g08025b1 gmp bits in 2.74s. */
if ((nx >= 2) && (ny >= 2) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(32*t1+ny+31,32),floord(tmax+ny-1,32));t2++) {
      for (t3=max(ceild(32*t2-ny-30,32),t1);t3<=min(min(floord(32*t1+nx+31,32),floord(32*t2+nx+30,32)),floord(tmax+nx-1,32));t3++) {
        if ((t1 == t3) && (t1 <= floord(32*t2-ny,32))) {
          for (t6=32*t2-ny+1;t6<=min(32*t1+31,32*t2-ny+nx);t6++) {
            hz[-32*t2+t6+ny-1][ny-1]=hz[-32*t2+t6+ny-1][ny-1]-0.7*(ex[-32*t2+t6+ny-1][ny-1 +1]-ex[-32*t2+t6+ny-1][ny-1]+ey[-32*t2+t6+ny-1 +1][ny-1]-ey[-32*t2+t6+ny-1][ny-1]);;
          }
        }
        if ((t1 <= min(floord(32*t2-ny,32),t3-1)) && (t2 >= ceild(32*t3+ny-nx+1,32))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t2-ny+nx);t6++) {
            hz[-32*t2+t6+ny-1][ny-1]=hz[-32*t2+t6+ny-1][ny-1]-0.7*(ex[-32*t2+t6+ny-1][ny-1 +1]-ex[-32*t2+t6+ny-1][ny-1]+ey[-32*t2+t6+ny-1 +1][ny-1]-ey[-32*t2+t6+ny-1][ny-1]);;
          }
        }
        if ((t1 <= min(floord(32*t2-ny,32),floord(32*t2-ny+nx-32,32))) && (32*t2 == 32*t3+ny-nx)) {
          if ((31*ny+nx)%32 == 0) {
            hz[nx-1][ny-1]=hz[nx-1][ny-1]-0.7*(ex[nx-1][ny-1 +1]-ex[nx-1][ny-1]+ey[nx-1 +1][ny-1]-ey[nx-1][ny-1]);;
          }
        }
        if ((t1 <= floord(32*t3-nx,32)) && (t2 <= floord(32*t3+ny-nx-1,32))) {
          for (t5=max(32*t2,32*t3-nx+1);t5<=min(32*t2+31,32*t3+ny-nx);t5++) {
            hz[nx-1][-32*t3+t5+nx-1]=hz[nx-1][-32*t3+t5+nx-1]-0.7*(ex[nx-1][-32*t3+t5+nx-1 +1]-ex[nx-1][-32*t3+t5+nx-1]+ey[nx-1 +1][-32*t3+t5+nx-1]-ey[nx-1][-32*t3+t5+nx-1]);;
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=32*t1;t4<=min(min(tmax-1,32*t1-nx+31),32*t1-ny+31);t4++) {
            ey[0][0]=t4;;
            for (t6=t4+1;t6<=t4+nx-1;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
            for (t6=t4+1;t6<=t4+nx;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(32*t1,32*t1-ny+32);t4<=min(tmax-1,32*t1-nx+31);t4++) {
            ey[0][0]=t4;;
            for (t6=t4+1;t6<=t4+nx-1;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(32*t1,32*t1-nx+32);t4<=min(tmax-1,32*t1-ny+31);t4++) {
            ey[0][0]=t4;;
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(max(32*t1,32*t1-nx+32),32*t1-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            ey[0][0]=t4;;
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+1);t4<=min(min(tmax-1,32*t1-nx+31),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
            for (t6=t4+1;t6<=t4+nx;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+32);t4<=min(tmax-1,32*t1-nx+31);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(max(32*t1,32*t1-nx+32),32*t2-ny+1);t4<=min(min(32*t1+30,tmax-1),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
            for (t6=t4+1;t6<=32*t1+31;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(max(32*t1,32*t1-nx+32),32*t2-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              for (t6=t4+1;t6<=32*t1+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
          }
        }
        if ((t1 == t3) && (t1 <= min(floord(tmax-32,32),t2-1))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1+ny+30);t5++) {
            ey[0][-32*t1+t5-31]=32*t1+31;;
            ex[0][-32*t1+t5-31]=ex[0][-32*t1+t5-31]-0.5*(hz[0][-32*t1+t5-31]-hz[0][-32*t1+t5-31 -1]);;
          }
        }
        if ((t1 == t2) && (t1 == t3) && (t1 <= floord(tmax-32,32))) {
          ey[0][0]=32*t1+31;;
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(32*t1,32*t3-nx+1);t4<=min(min(tmax-1,32*t1-ny+31),32*t3-nx+31);t4++) {
            for (t6=32*t3;t6<=t4+nx-1;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
            for (t6=32*t3;t6<=t4+nx;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(max(32*t1,32*t1-ny+32),32*t3-nx+1);t4<=min(min(32*t1+30,tmax-1),32*t3-nx+31);t4++) {
            for (t6=32*t3;t6<=t4+nx-1;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(32*t1,32*t3-nx+32);t4<=min(tmax-1,32*t1-ny+31);t4++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(max(32*t1,32*t1-ny+32),32*t3-nx+32);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+1),32*t3-nx+1);t4<=min(min(min(32*t1+31,tmax-1),32*t2-ny+31),32*t3-nx+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
            for (t6=32*t3;t6<=t4+nx;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+32),32*t3-nx+1);t4<=min(min(32*t1+31,tmax-1),32*t3-nx+31);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              for (t6=32*t3;t6<=t4+nx-1;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
              hz[nx-1][-t4+t5-1]=hz[nx-1][-t4+t5-1]-0.7*(ex[nx-1][-t4+t5-1 +1]-ex[nx-1][-t4+t5-1]+ey[nx-1 +1][-t4+t5-1]-ey[nx-1][-t4+t5-1]);;
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+1),32*t3-nx+32);t4<=min(min(32*t1+31,tmax-1),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              hz[-t4+t6-1][ny-1]=hz[-t4+t6-1][ny-1]-0.7*(ex[-t4+t6-1][ny-1 +1]-ex[-t4+t6-1][ny-1]+ey[-t4+t6-1 +1][ny-1]-ey[-t4+t6-1][ny-1]);;
            }
          }
        }
        if (t1 <= min(t2-1,t3-1)) {
          for (t4=max(max(32*t1,32*t2-ny+32),32*t3-nx+32);t4<=min(32*t1+31,tmax-1);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                ey[-t4+t6][-t4+t5]=ey[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6-1][-t4+t5]);;
                ex[-t4+t6][-t4+t5]=ex[-t4+t6][-t4+t5]-0.5*(hz[-t4+t6][-t4+t5]-hz[-t4+t6][-t4+t5-1]);;
                hz[-t4+t6-1][-t4+t5-1]=hz[-t4+t6-1][-t4+t5-1]-0.7*(ex[-t4+t6-1][-t4+t5-1 +1]-ex[-t4+t6-1][-t4+t5-1]+ey[-t4+t6-1 +1][-t4+t5-1]-ey[-t4+t6-1][-t4+t5-1]);;
              }
            }
          }
        }
        if ((t1 == t2) && (t1 <= min(floord(tmax-32,32),t3-1))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t1+nx+30);t6++) {
            ey[-32*t1+t6-31][0]=ey[-32*t1+t6-31][0]-0.5*(hz[-32*t1+t6-31][0]-hz[-32*t1+t6-31 -1][0]);;
          }
        }
      }
    }
  }
}
if ((nx >= 2) && (ny == 1) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(tmax,32),t1+1);t2++) {
      for (t3=t2;t3<=min(min(floord(32*t1+nx+31,32),floord(32*t2+nx+30,32)),floord(tmax+nx-1,32));t3++) {
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=32*t1;t4<=min(32*t1+30,tmax-1);t4++) {
            ey[0][0]=t4;;
            for (t6=t4+1;t6<=min(32*t1+31,t4+nx-1);t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t6=t4+1;t6<=min(32*t1+31,t4+nx);t6++) {
              hz[-t4+t6-1][0]=hz[-t4+t6-1][0]-0.7*(ex[-t4+t6-1][0 +1]-ex[-t4+t6-1][0]+ey[-t4+t6-1 +1][0]-ey[-t4+t6-1][0]);;
            }
          }
        }
        if ((t1 == t2) && (t1 == t3) && (t1 <= floord(tmax-32,32))) {
          ey[0][0]=32*t1+31;;
        }
        if ((t1 == t2) && (t1 <= floord(32*t3-nx,32))) {
          hz[nx-1][0]=hz[nx-1][0]-0.7*(ex[nx-1][0 +1]-ex[nx-1][0]+ey[nx-1 +1][0]-ey[nx-1][0]);;
        }
        if ((t1 == t2-1) && (32*t1 == 32*t3-nx-31)) {
          if ((nx+31)%32 == 0) {
            hz[nx-1][0]=hz[nx-1][0]-0.7*(ex[nx-1][0 +1]-ex[nx-1][0]+ey[nx-1 +1][0]-ey[nx-1][0]);;
          }
        }
        if ((t1 == t2-1) && (t1 >= ceild(32*t3-nx-30,32))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t1+nx+31);t6++) {
            hz[-32*t1+t6-32][0]=hz[-32*t1+t6-32][0]-0.7*(ex[-32*t1+t6-32][0 +1]-ex[-32*t1+t6-32][0]+ey[-32*t1+t6-32 +1][0]-ey[-32*t1+t6-32][0]);;
          }
        }
        if ((t1 == t2) && (t1 <= t3-1)) {
          for (t4=max(32*t1,32*t3-nx+1);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t6=32*t3;t6<=min(32*t3+31,t4+nx-1);t6++) {
              ey[-t4+t6][0]=ey[-t4+t6][0]-0.5*(hz[-t4+t6][0]-hz[-t4+t6-1][0]);;
            }
            for (t6=32*t3;t6<=min(32*t3+31,t4+nx);t6++) {
              hz[-t4+t6-1][0]=hz[-t4+t6-1][0]-0.7*(ex[-t4+t6-1][0 +1]-ex[-t4+t6-1][0]+ey[-t4+t6-1 +1][0]-ey[-t4+t6-1][0]);;
            }
          }
        }
        if ((t1 == t2) && (t1 <= min(floord(tmax-32,32),t3-1))) {
          for (t6=32*t3;t6<=min(32*t3+31,32*t1+nx+30);t6++) {
            ey[-32*t1+t6-31][0]=ey[-32*t1+t6-31][0]-0.5*(hz[-32*t1+t6-31][0]-hz[-32*t1+t6-31 -1][0]);;
          }
        }
      }
    }
  }
}
if ((nx == 1) && (ny >= 2) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(32*t1+ny+31,32),floord(tmax+ny-1,32));t2++) {
      for (t3=max(ceild(32*t2-ny-30,32),t1);t3<=min(min(floord(tmax,32),t2),t1+1);t3++) {
        if ((t1 == t3) && (t1 <= floord(32*t2-ny,32))) {
          hz[0][ny-1]=hz[0][ny-1]-0.7*(ex[0][ny-1 +1]-ex[0][ny-1]+ey[0 +1][ny-1]-ey[0][ny-1]);;
        }
        if ((t1 == t3-1) && (32*t1 == 32*t2-ny-31)) {
          if ((ny+31)%32 == 0) {
            hz[0][ny-1]=hz[0][ny-1]-0.7*(ex[0][ny-1 +1]-ex[0][ny-1]+ey[0 +1][ny-1]-ey[0][ny-1]);;
          }
        }
        if ((t1 == t3-1) && (t1 >= ceild(32*t2-ny-30,32))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1+ny+31);t5++) {
            hz[0][-32*t1+t5-32]=hz[0][-32*t1+t5-32]-0.7*(ex[0][-32*t1+t5-32 +1]-ex[0][-32*t1+t5-32]+ey[0 +1][-32*t1+t5-32]-ey[0][-32*t1+t5-32]);;
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=32*t1;t4<=min(tmax-1,32*t1-ny+31);t4++) {
            ey[0][0]=t4;;
            for (t5=t4+1;t5<=t4+ny-1;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              hz[0][-t4+t5-1]=hz[0][-t4+t5-1]-0.7*(ex[0][-t4+t5-1 +1]-ex[0][-t4+t5-1]+ey[0 +1][-t4+t5-1]-ey[0][-t4+t5-1]);;
            }
            hz[0][ny-1]=hz[0][ny-1]-0.7*(ex[0][ny-1 +1]-ex[0][ny-1]+ey[0 +1][ny-1]-ey[0][ny-1]);;
          }
        }
        if ((t1 == t2) && (t1 == t3)) {
          for (t4=max(32*t1,32*t1-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            ey[0][0]=t4;;
            for (t5=t4+1;t5<=32*t1+31;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              hz[0][-t4+t5-1]=hz[0][-t4+t5-1]-0.7*(ex[0][-t4+t5-1 +1]-ex[0][-t4+t5-1]+ey[0 +1][-t4+t5-1]-ey[0][-t4+t5-1]);;
            }
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+1);t4<=min(min(32*t1+30,tmax-1),32*t2-ny+31);t4++) {
            for (t5=32*t2;t5<=t4+ny-1;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              hz[0][-t4+t5-1]=hz[0][-t4+t5-1]-0.7*(ex[0][-t4+t5-1 +1]-ex[0][-t4+t5-1]+ey[0 +1][-t4+t5-1]-ey[0][-t4+t5-1]);;
            }
            hz[0][ny-1]=hz[0][ny-1]-0.7*(ex[0][ny-1 +1]-ex[0][ny-1]+ey[0 +1][ny-1]-ey[0][ny-1]);;
          }
        }
        if ((t1 == t3) && (t1 <= t2-1)) {
          for (t4=max(32*t1,32*t2-ny+32);t4<=min(32*t1+30,tmax-1);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              ey[0][-t4+t5]=t4;;
              ex[0][-t4+t5]=ex[0][-t4+t5]-0.5*(hz[0][-t4+t5]-hz[0][-t4+t5-1]);;
              hz[0][-t4+t5-1]=hz[0][-t4+t5-1]-0.7*(ex[0][-t4+t5-1 +1]-ex[0][-t4+t5-1]+ey[0 +1][-t4+t5-1]-ey[0][-t4+t5-1]);;
            }
          }
        }
        if ((t1 == t3) && (t1 <= min(floord(tmax-32,32),t2-1))) {
          for (t5=32*t2;t5<=min(32*t2+31,32*t1+ny+30);t5++) {
            ey[0][-32*t1+t5-31]=32*t1+31;;
            ex[0][-32*t1+t5-31]=ex[0][-32*t1+t5-31]-0.5*(hz[0][-32*t1+t5-31]-hz[0][-32*t1+t5-31 -1]);;
          }
        }
        if ((t1 == t2) && (t1 == t3) && (t1 <= floord(tmax-32,32))) {
          ey[0][0]=32*t1+31;;
        }
      }
    }
  }
}
if ((nx <= 0) && (ny >= 2) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(32*t1+ny+30,32),floord(tmax+ny-2,32));t2++) {
      for (t4=max(32*t1,32*t2-ny+1);t4<=min(32*t1+31,tmax-1);t4++) {
        for (t5=max(32*t2,t4);t5<=min(32*t2+31,t4+ny-1);t5++) {
          ey[0][-t4+t5]=t4;;
        }
      }
    }
  }
}
if ((nx == 1) && (ny == 1) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t2=t1;t2<=min(floord(tmax,32),t1+1);t2++) {
      if (t1 == t2) {
        for (t4=32*t1;t4<=min(32*t1+30,tmax-1);t4++) {
          ey[0][0]=t4;;
          hz[0][0]=hz[0][0]-0.7*(ex[0][0 +1]-ex[0][0]+ey[0 +1][0]-ey[0][0]);;
        }
      }
      if ((t1 == t2) && (t1 <= floord(tmax-32,32))) {
        ey[0][0]=32*t1+31;;
      }
      if (t1 == t2-1) {
        hz[0][0]=hz[0][0]-0.7*(ex[0][0 +1]-ex[0][0]+ey[0 +1][0]-ey[0][0]);;
      }
    }
  }
}
if ((nx <= 0) && (ny == 1) && (tmax >= 1)) {
  for (t1=0;t1<=floord(tmax-1,32);t1++) {
    for (t4=32*t1;t4<=min(32*t1+31,tmax-1);t4++) {
      ey[0][0]=t4;;
    }
  }
}
/* End of CLooG code */

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
