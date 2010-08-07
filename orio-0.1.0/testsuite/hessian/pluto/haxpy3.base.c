

#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define N SIZE
double X0[N][N];
double X1[N][N];
double X2[N][N];
double Y[N][N];
double u0[N];
double u1[N];
double u2[N];
#define a0 32.12
#define a1 3322.12
#define a2 1.123
#define b00 1321.9
#define b01 21.55
#define b02 10.3
#define b11 1210.313
#define b12 9.373
#define b22 1992.31221

void init_arrays()
{
  int i1, i2;
  for (i1=0; i1<N; i1++)
   for (i2=0; i2<N; i2++)
    X0[i1][i2] = (i1+i2) % 5 + 1;
  for (i1=0; i1<N; i1++)
   for (i2=0; i2<N; i2++)
    X1[i1][i2] = (i1+i2) % 5 + 1;
  for (i1=0; i1<N; i1++)
   for (i2=0; i2<N; i2++)
    X2[i1][i2] = (i1+i2) % 5 + 1;
  for (i1=0; i1<N; i1++)
   for (i2=0; i2<N; i2++)
    Y[i1][i2] = 0;
  for (i1=0; i1<N; i1++)
   u0[i1] = (i1) % 5 + 1;
  for (i1=0; i1<N; i1++)
   u1[i1] = (i1) % 5 + 1;
  for (i1=0; i1<N; i1++)
   u2[i1] = (i1) % 5 + 1;
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
    
    int i,j,ii,jj,iii,jjj;
    
    {
      double X2_copy;
      double X1_copy;
      double X0_copy;
      double Y_copy;
      for (i=0; i<=N-1; i=i+1) 
	for (j=0; j<=N-1; j=j+1) {
	  Y[i][j]=
	    a0*X0[i][j]+a1*X1[i][j]+a2*X2[i][j]
	    +2.0*b00*u0[i]*u0[j]
	    +2.0*b11*u1[i]*u1[j]
	    +2.0*b22*u2[i]*u2[j]
	    +b01*(u0[i]*u1[j]+u1[i]*u0[j])
	    +b02*(u0[i]*u2[j]+u2[i]*u0[j])
	    +b12*(u1[i]*u2[j]+u2[i]*u1[j]);
	}
    }
    
    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return Y[0][0];
}
                                                                   
