#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#define N 2000000
#define T 1000

#pragma declarations
double a[N];
double b[N];
#pragma enddeclarations

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

#include "util.h"

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S1(zT0,zT1,t,i)	{b[i]=((double)(333))/1000*(a[1+i]+a[i]+a[i-1]);}
#define S2(zT0,zT1,t,j)	{a[j]=b[j];}

int main()
{
	assert(N >= 1000);
	assert(T >= 1000);
	int c1, c2, c3, c4, c5;
	int i, j, k, l, t;

    register int lb, ub;

#ifdef TEST
	init_array() ;
#endif

#ifdef PERFCTR
    PERF_INIT;
#endif
    

    /* Generated from jacobi-imper.par.cloog by CLooG v0.14.1 64 bits in 0.11s. */
    for (c1=-1;c1<=floord(2*N+5*T-8,2048);c1++) {
        lb = max(max(max(ceild(6144*c1-N-4092,10240),ceild(2048*c1-2047,4096)),0),ceild(2048*c1-N-2*T+4,2048));
        ub =  min(min(min(floord(2048*c1+2045,2048),floord(2048*c1+T+2047,4096)),floord(N+3*T-4,2048)),floord(6144*c1+6141,10240));
#pragma omp parallel for shared(c1,lb,ub,a,b) private(i,j,k,l,c2,c3,c4,c5) default(none) schedule(static)
        for (c2=lb;c2<=ub;c2++) {
            if ((c1 <= floord(4096*c2-T,2048)) && (c2 >= ceild(3*T,2048))) {
                c3 = 2048*c2-T ;
                c4 = 2048*c2 ;
                c5 = 1 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = T-1 ;
                l = 2048*c2-3*T+2 ;
                S2(-c1+2*c2,3*c1-5*c2,T-1,2048*c2-3*T+2) ;
            }
            for (c3=max(max(ceild(4096*c2,3),2),2048*c1-2048*c2);c3<=min(min(2048*c1-2048*c2+2047,floord(4096*c2+1,3)),2*T-1);c3++) {
                for (c4=max(2048*c2,c3+1);c4<=floord(3*c3,2);c4++) {
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(2048*c2,2),2048*c1-2048*c2);c3<=min(min(3,2048*c1-2048*c2+2047),2048*c2+2046);c3++) {
                c4 = c3 ;
                c5 = 0 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = 0 ;
                l = c3 ;
                S1(-c1+2*c2,3*c1-5*c2,0,c3) ;
                for (c4=c3+1;c4<=min(floord(3*c3,2),2048*c2+2047);c4++) {
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(ceild(4096*c2+2,3),2048*c2-T+1),2048*c1-2048*c2);c3<=min(min(min(min(2048*c2-1,floord(4096*c2+N-4,3)),2048*c1-2048*c2+2047),floord(4096*c2+4095,3)),2*T+1);c3++) {
                for (c4=2048*c2;c4<=min(c3+T-1,floord(3*c3-2,2));c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
                for (c4=ceild(3*c3-1,2);c4<=min(min(c3+T,floord(3*c3,2)),2048*c2+2047);c4++) {
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(2048*c2,4),2048*c1-2048*c2);c3<=min(min(min(N-2,2048*c1-2048*c2+2047),floord(4096*c2+4095,3)),2*T+1);c3++) {
                c4 = c3 ;
                c5 = 0 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = 0 ;
                l = c3 ;
                S1(-c1+2*c2,3*c1-5*c2,0,c3) ;
                for (c4=c3+1;c4<=min(c3+T-1,floord(3*c3-2,2));c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
                for (c4=ceild(3*c3-1,2);c4<=min(min(c3+T,floord(3*c3,2)),2048*c2+2047);c4++) {
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(2*T+2,2048*c2-T+1),2048*c1-2048*c2);c3<=min(min(min(floord(4096*c2+N-4,3),2048*c2-1),2048*c1-2048*c2+2047),2048*c2-T+2047);c3++) {
                for (c4=2048*c2;c4<=c3+T-1;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
                c4 = c3+T ;
                c5 = 1 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = T-1 ;
                l = c3-2*T+2 ;
                S2(-c1+2*c2,3*c1-5*c2,T-1,c3-2*T+2) ;
            }
            for (c3=max(2048*c1-2048*c2,ceild(4096*c2+4096,3));c3<=min(min(min(floord(4096*c2+N-4,3),2048*c2-1),2048*c1-2048*c2+2047),2*T+1);c3++) {
                for (c4=2048*c2;c4<=min(2048*c2+2047,c3+T-1);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(2048*c2,2*T+2),2048*c1-2048*c2);c3<=min(min(N-2,2048*c1-2048*c2+2047),2048*c2-T+2047);c3++) {
                c4 = c3 ;
                c5 = 0 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = 0 ;
                l = c3 ;
                S1(-c1+2*c2,3*c1-5*c2,0,c3) ;
                for (c4=c3+1;c4<=c3+T-1;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
                c4 = c3+T ;
                c5 = 1 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = T-1 ;
                l = c3-2*T+2 ;
                S2(-c1+2*c2,3*c1-5*c2,T-1,c3-2*T+2) ;
            }
            for (c3=max(max(2048*c2,2048*c1-2048*c2),ceild(4096*c2+4096,3));c3<=min(min(min(N-2,2048*c1-2048*c2+2047),2048*c2+2046),2*T+1);c3++) {
                c4 = c3 ;
                c5 = 0 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = 0 ;
                l = c3 ;
                S1(-c1+2*c2,3*c1-5*c2,0,c3) ;
                for (c4=c3+1;c4<=min(2048*c2+2047,c3+T-1);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(ceild(4096*c2+N-3,3),N-1),2048*c1-2048*c2);c3<=min(min(2048*c1-2048*c2+2047,floord(4096*c2+4095,3)),2*T+1);c3++) {
                for (c4=max(2048*c2,ceild(3*c3-N+2,2));c4<=floord(3*c3-N+3,2);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                }
                for (c4=ceild(3*c3-N+4,2);c4<=min(c3+T-1,floord(3*c3-2,2));c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
                for (c4=ceild(3*c3-1,2);c4<=min(min(c3+T,floord(3*c3,2)),2048*c2+2047);c4++) {
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(2048*c1-2048*c2,2*T+2),2048*c2-T+2048);c3<=min(min(floord(4096*c2+N-4,3),2048*c2-1),2048*c1-2048*c2+2047);c3++) {
                for (c4=2048*c2;c4<=2048*c2+2047;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(max(2048*c2,2048*c1-2048*c2),2*T+2),2048*c2-T+2048);c3<=min(min(N-2,2048*c1-2048*c2+2047),2048*c2+2046);c3++) {
                c4 = c3 ;
                c5 = 0 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = 0 ;
                l = c3 ;
                S1(-c1+2*c2,3*c1-5*c2,0,c3) ;
                for (c4=c3+1;c4<=2048*c2+2047;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(max(ceild(4096*c2+N-3,3),N-1),2*T+2),2048*c1-2048*c2);c3<=min(min(N+2*T-6,2048*c1-2048*c2+2047),2048*c2-T+2047);c3++) {
                for (c4=max(2048*c2,ceild(3*c3-N+2,2));c4<=floord(3*c3-N+3,2);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                }
                for (c4=ceild(3*c3-N+4,2);c4<=c3+T-1;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
                c4 = c3+T ;
                c5 = 1 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = T-1 ;
                l = c3-2*T+2 ;
                S2(-c1+2*c2,3*c1-5*c2,T-1,c3-2*T+2) ;
            }
            for (c3=max(max(max(ceild(4096*c2+N-3,3),2048*c1-2048*c2),N-1),ceild(4096*c2+4096,3));c3<=min(min(2048*c1-2048*c2+2047,floord(4096*c2+N+4090,3)),2*T+1);c3++) {
                for (c4=max(2048*c2,ceild(3*c3-N+2,2));c4<=floord(3*c3-N+3,2);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                }
                for (c4=ceild(3*c3-N+4,2);c4<=min(2048*c2+2047,c3+T-1);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(max(max(ceild(4096*c2+N-3,3),2048*c1-2048*c2),N-1),2*T+2),2048*c2-T+2048);c3<=min(2048*c1-2048*c2+2047,floord(4096*c2+N+4090,3));c3++) {
                for (c4=max(2048*c2,ceild(3*c3-N+2,2));c4<=floord(3*c3-N+3,2);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                }
                for (c4=ceild(3*c3-N+4,2);c4<=2048*c2+2047;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                    c5 = 1 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4-1 ;
                    l = 3*c3-2*c4+2 ;
                    S2(-c1+2*c2,3*c1-5*c2,-c3+c4-1,3*c3-2*c4+2) ;
                }
            }
            for (c3=max(max(N+2*T-5,2048*c2-T+1),2048*c1-2048*c2);c3<=min(min(2048*c1-2048*c2+2047,2048*c2-T+2047),N+2*T-4);c3++) {
                for (c4=max(2048*c2,ceild(3*c3-N+2,2));c4<=c3+T-1;c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                }
                c4 = c3+T ;
                c5 = 1 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = T-1 ;
                l = c3-2*T+2 ;
                S2(-c1+2*c2,3*c1-5*c2,T-1,c3-2*T+2) ;
            }
            if ((c1 >= 2*c2) && (c2 <= floord(N-2049,2048))) {
                c3 = 2048*c2+2047 ;
                c4 = 2048*c2+2047 ;
                c5 = 0 ;
                i = -c1+2*c2 ;
                j = 3*c1-5*c2 ;
                k = 0 ;
                l = 2048*c2+2047 ;
                S1(-c1+2*c2,3*c1-5*c2,0,2048*c2+2047) ;
            }
            for (c3=max(max(ceild(4096*c2+N+4091,3),2048*c1-2048*c2),N-1);c3<=min(min(2048*c1-2048*c2+2047,floord(4096*c2+N+4092,3)),N+2*T-4);c3++) {
                for (c4=ceild(3*c3-N+2,2);c4<=min(2048*c2+2047,c3+T-1);c4++) {
                    c5 = 0 ;
                    i = -c1+2*c2 ;
                    j = 3*c1-5*c2 ;
                    k = -c3+c4 ;
                    l = 3*c3-2*c4 ;
                    S1(-c1+2*c2,3*c1-5*c2,-c3+c4,3*c3-2*c4) ;
                }
            }
        }
    }

#ifdef PERFCTR
    PERF_EXIT;
#endif

#ifdef TEST
    print_array();
#endif
    return 0;
}
