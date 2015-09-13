#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#ifdef PERFCTR
#include <papi.h>
#endif

#define N 2000000
#define T 1000

#pragma declarations
double a[N];
double b[N];
#pragma enddeclarations

double t_start, t_end;

void init_array()
{
    int j;

    for (j=0; j<N; j++) {
        a[j] = ((double)j)/N;
    }
}


void print_array()
{
    int j;

    for (j=0; j<N; j++) {
        fprintf(stdout, "%lf ", a[j]);
        if (j%80 == 20) fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}


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
    PERFCTR_INIT;
#endif

/* Generated from jacobi-imper.sched.cloog by CLooG v0.14.1 64 bits in 0.01s. */
for (c1=-1;c1<=floord(N+3*T-4,2048);c1++) {
    lb = max(max(ceild(2048*c1-T+1,2048),ceild(4096*c1-2045,6144)),0);
    ub = min(min(floord(2048*c1+2047,2048),floord(4096*c1+N+4093,6144)),floord(N+2*T-3,2048));
#pragma omp parallel for shared(c1,lb,ub,a,b) private(c2,c3,c4,c5,i,j,k,l)
  for (c2=lb;c2<=ub;c2++) {
    if (c1 >= max(c2,ceild(6144*c2-N+2,4096))) {
      c3 = 4096*c1-4096*c2 ;
      for (c4=max(4096*c1-4096*c2+2,2048*c2);c4<=min(4096*c1-4096*c2+N-2,2048*c2+2047);c4++) {
        c5 = 0 ;
        if ((c1-c2)%2 == 0) {
          i = (c1-c2)/2 ;
          j = -c1+2*c2 ;
          k = 2048*c1-2048*c2 ;
          l = -4096*c1+4096*c2+c4 ;
          S1((c1-c2)/2,-c1+2*c2,2048*c1-2048*c2,-4096*c1+4096*c2+c4) ;
        }
      }
    }
    if ((c1 <= floord(4096*c2-1,4096)) && (c2 <= floord(N-2,2048))) {
      c3 = 0 ;
      for (c4=max(2048*c2,2);c4<=min(2048*c2+2047,N-2);c4++) {
        c5 = 0 ;
        if ((c1-c2)%2 == 0) {
          i = (c1-c2)/2 ;
          j = -c1+2*c2 ;
          k = 0 ;
          l = c4 ;
          S1((c1-c2)/2,-c1+2*c2,0,c4) ;
        }
      }
    }
    for (c3=max(max(1,2048*c2-N+2),4096*c1-4096*c2+1);c3<=min(min(4096*c1-4096*c2+4094,2048*c2+2045),2*T-2);c3++) {
      for (c4=max(2048*c2,c3+2);c4<=min(c3+N-2,2048*c2+2047);c4++) {
        c5 = 0 ;
        if ((c1-c2)%2 == 0) {
          i = (c1-c2)/2 ;
          j = -c1+2*c2 ;
          if (c3%2 == 0) {
            k = c3/2 ;
            l = -c3+c4 ;
            S1((c1-c2)/2,-c1+2*c2,c3/2,-c3+c4) ;
          }
        }
        c5 = 1 ;
        if ((c1-c2)%2 == 0) {
          i = (c1-c2)/2 ;
          j = -c1+2*c2 ;
          if ((c3-1)%2 == 0) {
            k = (c3-1)/2 ;
            l = -c3+c4 ;
            S2((c1-c2)/2,-c1+2*c2,(c3-1)/2,-c3+c4) ;
          }
        }
      }
    }
    if (c1 <= min(floord(3072*c2-1025,2048),floord(2048*c2+T-2048,2048))) {
      c3 = 4096*c1-4096*c2+4095 ;
      for (c4=max(4096*c1-4096*c2+4097,2048*c2);c4<=min(4096*c1-4096*c2+N+4093,2048*c2+2047);c4++) {
        c5 = 1 ;
        if ((c1-c2)%2 == 0) {
          i = (c1-c2)/2 ;
          j = -c1+2*c2 ;
          k = 2048*c1-2048*c2+2047 ;
          l = -4096*c1+4096*c2+c4-4095 ;
          S2((c1-c2)/2,-c1+2*c2,2048*c1-2048*c2+2047,-4096*c1+4096*c2+c4-4095) ;
        }
      }
    }
    if ((c1 >= ceild(4096*c2+2*T-4095,4096)) && (c2 >= ceild(T-1023,1024))) {
      c3 = 2*T-1 ;
      for (c4=max(2*T+1,2048*c2);c4<=min(2048*c2+2047,N+2*T-3);c4++) {
        c5 = 1 ;
        if ((c1-c2)%2 == 0) {
          i = (c1-c2)/2 ;
          j = -c1+2*c2 ;
          k = T-1 ;
          l = c4-2*T+1 ;
          S2((c1-c2)/2,-c1+2*c2,T-1,c4-2*T+1) ;
        }
      }
    }
  }
}

#ifdef PERFCTR
    PERFCTR_EXIT;
#endif

#ifdef TEST
	print_array();
#endif
	return 0;
}
