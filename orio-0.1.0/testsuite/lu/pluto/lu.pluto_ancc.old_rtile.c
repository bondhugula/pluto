

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

double L[N][N];
double U[N][N];
double A[N][N +13];

void init_arrays()
{
    int i, j, k;

    /* have to initialize this matrix properly to prevent
    * division by zero
    */
    for (i=0; i<N; i++) {
     for (j=0; j<N; j++) {
      L[i][j] = 0.0;
      U[i][j] = 0.0;
     }
    }
    
    for (i=0; i<N; i++) {
     for (j=0; j<=i; j++) {
      L[i][j] = i+j+1;
      U[j][i] = i+j+1;
     }
    }
    
    for (i=0; i<N; i++) {
     for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
       A[i][j] += L[i][k]*U[k][j]; 
      }
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
    
   

register int i,j,k;
register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t;
register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6,
  newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;
register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6,
  newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;


/*@ begin PolySyn(
  l1_tiles = [T1_1,T1_2,T1_3];
  l2_tiles = [T2_1,T2_2,T2_3];
  hotspot_permut = PERM_B;
  unroll_factors = [U1,U2,U3];
  parallelize = PAR;
  scalar_replace = SCREP;
  icc_vectorize = IVEC;
) @*/


	int c1, c2, c3, c4, c5, c6, c7, c8, c9;

	register int lb, ub, lb1, ub1, lb2, ub2;
/* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 2.36s. */
for (c1=-1;c1<=floord(3*N-5,128);c1++) {
	lb1=max(max(ceild(64*c1-N+2,64),ceild(32*c1-63,96)),0);
	ub1=min(floord(N-1,128),floord(64*c1+63,64));
#pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5,c6,c7,c8,c9)
	for (c2=lb1; c2<=ub1; c2++) {
    for (c3=max(ceild(32*c1-32*c2-1953,2016),ceild(32*c1-32*c2-31,32));c3<=floord(N-1,64);c3++) {
      for (c4=max(max(2*c1-2*c2-64*c3-62,2*c1-2*c2),0);c4<=min(min(min(min(floord(64*c2+63,16),2*c1-2*c2+1),floord(992*c3+961,16)),floord(32*c3+31,16)),floord(N-2,32));c4++) {
        for (c5=max(max(ceild(16*c4-15,16),0),4*c2);c5<=min(4*c2+3,floord(N-1,32));c5++) {
          for (c6=max(max(max(max(ceild(16*c4-465,496),ceild(2*c1-2*c2-2*c3-c4-31,31)),ceild(-2*c1+2*c2+2*c3+c4-31,33)),2*c3),ceild(16*c4-15,16));c6<=min(2*c3+1,floord(N-1,32));c6++) {
            if ((c1 == c2+c3) && (c4 == c6)) {
              for (c7=max(0,32*c6);c7<=min(min(32*c6+30,N-2),32*c5+30);c7++) {
                for (c8=max(32*c5,c7+1);c8<=min(N-1,32*c5+31);c8++) {
                  A[c7][c8]=A[c7][c8]/A[c7][c7] ;
                  for (c9=c7+1;c9<=min(N-1,32*c6+31);c9++) {
                    A[c9][c8]=A[c9][c8]-A[c9][c7]*A[c7][c8] ;
                  }
                }
              }
            }
            
/*@ begin Loop(
transform Composite(
  regtile = (['c7', 'c8', 'c9'],[8, 8, 1]),
  permut = [(['c7'],['c8'],['c9'])],
  scalarreplace = (False, 'double'))
for (c7=max(32*c4,0);c7<=min(min(32*c6-1,32*c4+31),32*c5+30);c7++) {

              for (c8=max(c7+1,32*c5);c8<=min(32*c5+31,N-1);c8++) 
                for (c9=32*c6;c9<=min(N-1,32*c6+31);c9++) 
{
                  A[c9][c8]=A[c9][c8]-A[c9][c7]*A[c7][c8] ;
}

            }
) @*/{
  for (c7t=max(32*c4,0); c7t<=min(min(32*c6-1,32*c4+31),32*c5+30)-7; c7t=c7t+8) {
    for (c8t=max(c7t+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][c7t]*A[c7t][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][c7t]*A[c7t][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][c7t]*A[c7t][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][c7t]*A[c7t][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][c7t]*A[c7t][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][c7t]*A[c7t][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][c7t]*A[c7t][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][c7t]*A[c7t][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][c7t]*A[c7t][c8];
      }
    for (c8t=max((c7t+1)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+1)]*A[(c7t+1)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+1)]*A[(c7t+1)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+1)]*A[(c7t+1)][c8];
      }
    for (c8t=max((c7t+2)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+2)]*A[(c7t+2)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+2)]*A[(c7t+2)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+2)]*A[(c7t+2)][c8];
      }
    for (c8t=max((c7t+3)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+3)]*A[(c7t+3)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+3)]*A[(c7t+3)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+3)]*A[(c7t+3)][c8];
      }
    for (c8t=max((c7t+4)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+4)]*A[(c7t+4)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+4)]*A[(c7t+4)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+4)]*A[(c7t+4)][c8];
      }
    for (c8t=max((c7t+5)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+5)]*A[(c7t+5)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+5)]*A[(c7t+5)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+5)]*A[(c7t+5)][c8];
      }
    for (c8t=max((c7t+6)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+6)]*A[(c7t+6)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+6)]*A[(c7t+6)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+6)]*A[(c7t+6)][c8];
      }
    for (c8t=max((c7t+7)+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][(c7t+7)]*A[(c7t+7)][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][(c7t+7)]*A[(c7t+7)][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][(c7t+7)]*A[(c7t+7)][c8];
      }
  }
  for (c7=c7t; c7<=min(min(32*c6-1,32*c4+31),32*c5+30); c7=c7+1) {
    for (c8t=max(c7+1,32*c5); c8t<=min(32*c5+31,N-1)-7; c8t=c8t+8) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8t]=A[c9][c8t]-A[c9][c7]*A[c7][c8t];
        A[c9][(c8t+1)]=A[c9][(c8t+1)]-A[c9][c7]*A[c7][(c8t+1)];
        A[c9][(c8t+2)]=A[c9][(c8t+2)]-A[c9][c7]*A[c7][(c8t+2)];
        A[c9][(c8t+3)]=A[c9][(c8t+3)]-A[c9][c7]*A[c7][(c8t+3)];
        A[c9][(c8t+4)]=A[c9][(c8t+4)]-A[c9][c7]*A[c7][(c8t+4)];
        A[c9][(c8t+5)]=A[c9][(c8t+5)]-A[c9][c7]*A[c7][(c8t+5)];
        A[c9][(c8t+6)]=A[c9][(c8t+6)]-A[c9][c7]*A[c7][(c8t+6)];
        A[c9][(c8t+7)]=A[c9][(c8t+7)]-A[c9][c7]*A[c7][(c8t+7)];
      }
    for (c8=c8t; c8<=min(32*c5+31,N-1); c8=c8+1) 
      for (c9=32*c6; c9<=min(N-1,32*c6+31); c9++ ) {
        A[c9][c8]=A[c9][c8]-A[c9][c7]*A[c7][c8];
      }
  }
}
/*@ end @*/


            if ((c1 == c2+c3) && (-c4 == -c6) && (c4 <= min(floord(N-33,32),floord(32*c5-1,32)))) {
              for (c8=max(32*c5,32*c4+32);c8<=min(N-1,32*c5+31);c8++) {
                A[32*c4+31][c8]=A[32*c4+31][c8]/A[32*c4+31][32*c4+31] ;
              }
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */

/*@ end @*/


    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return ((int) A[0][0]); 

}
                                    
