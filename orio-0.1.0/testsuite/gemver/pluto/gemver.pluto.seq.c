

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define alpha 1
#define beta 1
double A[N][N +20];
double B[N][N +20];
double x[N];
double u1[N];
double u2[N];
double v1[N];
double v2[N];
double w[N];
double y[N];
double z[N];

void init_arrays()
{
    int i, j;
    for (i=0; i<N; i++) {
        u1[i] = i;
        u2[i] = (i+1)/N/2.0;
        v1[i] = (i+1)/N/4.0;
        v2[i] = (i+1)/N/6.0;
        y[i] = (i+1)/N/8.0;
        z[i] = (i+1)/N/9.0;
        x[i] = 0.0;
        w[i] = 0.0;
        for (j=0; j<N; j++) {
            A[i][j] = ((double) i*j)/N;
        }
    }
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




#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

	#define S1(zT0,zT1,zT2,zT3,i,j)	{B[i][j]=u2[i]*v2[j]+u1[i]*v1[j]+A[i][j];}
	#define S2(zT0,zT1,zT2,zT3,i,j)	{x[i]=beta*B[j][i]*y[j]+x[i];}
	#define S3(i)	{x[i]=z[i]+x[i];}
	#define S4(i,j)	{w[i]=alpha*B[i][j]*x[j]+w[i];}

	int c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;

	register int lbv, ubv;

/* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 0.05s. */
for (c2=0;c2<=floord(N-1,256);c2++) {
  for (c3=0;c3<=floord(N-1,256);c3++) {
    for (c4=max(0,8*c2);c4<=min(8*c2+7,floord(N-1,32));c4++) {
      for (c5=max(8*c3,0);c5<=min(floord(N-1,32),8*c3+7);c5++) {
/*@ begin Loop(
transform UnrollJam(ufactor=32)
        for (c6=max(32*c5,0);c6<=min(N-1,32*c5+31);c6++) 
{
{
	lbv=max(32*c4,0);
	ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
	for (c7=lbv; c7<=ubv; c7++) {
            S1(c3,c2,c5,c4,c6,c7) ;
            S2(c2,c3,c4,c5,c7,c6) ;
          }
}
}
) @*/{ 

  for (c6 = max(32 * c5, 0); c6 <= min(N - 1, 32 * c5 + 31) - 31; c6 = c6 + 32) 
{
	lbv=max(32*c4,0);
	ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
	for (c7=lbv; c7<=ubv; c7++) {
        S1(c3, c2, c5, c4, c6, c7); 
        S2(c2, c3, c4, c5, c7, c6); 
        S1(c3, c2, c5, c4, (c6 + 1), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 1)); 
        S1(c3, c2, c5, c4, (c6 + 2), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 2)); 
        S1(c3, c2, c5, c4, (c6 + 3), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 3)); 
        S1(c3, c2, c5, c4, (c6 + 4), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 4)); 
        S1(c3, c2, c5, c4, (c6 + 5), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 5)); 
        S1(c3, c2, c5, c4, (c6 + 6), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 6)); 
        S1(c3, c2, c5, c4, (c6 + 7), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 7)); 
        S1(c3, c2, c5, c4, (c6 + 8), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 8)); 
        S1(c3, c2, c5, c4, (c6 + 9), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 9)); 
        S1(c3, c2, c5, c4, (c6 + 10), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 10)); 
        S1(c3, c2, c5, c4, (c6 + 11), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 11)); 
        S1(c3, c2, c5, c4, (c6 + 12), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 12)); 
        S1(c3, c2, c5, c4, (c6 + 13), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 13)); 
        S1(c3, c2, c5, c4, (c6 + 14), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 14)); 
        S1(c3, c2, c5, c4, (c6 + 15), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 15)); 
        S1(c3, c2, c5, c4, (c6 + 16), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 16)); 
        S1(c3, c2, c5, c4, (c6 + 17), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 17)); 
        S1(c3, c2, c5, c4, (c6 + 18), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 18)); 
        S1(c3, c2, c5, c4, (c6 + 19), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 19)); 
        S1(c3, c2, c5, c4, (c6 + 20), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 20)); 
        S1(c3, c2, c5, c4, (c6 + 21), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 21)); 
        S1(c3, c2, c5, c4, (c6 + 22), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 22)); 
        S1(c3, c2, c5, c4, (c6 + 23), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 23)); 
        S1(c3, c2, c5, c4, (c6 + 24), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 24)); 
        S1(c3, c2, c5, c4, (c6 + 25), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 25)); 
        S1(c3, c2, c5, c4, (c6 + 26), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 26)); 
        S1(c3, c2, c5, c4, (c6 + 27), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 27)); 
        S1(c3, c2, c5, c4, (c6 + 28), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 28)); 
        S1(c3, c2, c5, c4, (c6 + 29), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 29)); 
        S1(c3, c2, c5, c4, (c6 + 30), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 30)); 
        S1(c3, c2, c5, c4, (c6 + 31), c7); 
        S2(c2, c3, c4, c5, c7, (c6 + 31)); 
      } 
}

  for (; c6 <= min(N - 1, 32 * c5 + 31); c6 = c6 + 1) 
{
	lbv=max(32*c4,0);
	ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
	for (c7=lbv; c7<=ubv; c7++) {
        S1(c3, c2, c5, c4, c6, c7); 
        S2(c2, c3, c4, c5, c7, c6); 
      } 
}
} 
/*@ end @*/
      }
    }
  }
}
for (c2=0;c2<=N-1;c2++) {
  S3(c2) ;
}
for (c2=0;c2<=N-1;c2++) {
  for (c3=0;c3<=N-1;c3++) {
    S4(c2,c3) ;
  }
}
/* End of CLooG code */

    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);

  return ((int) w[0]); 

}
                                    


