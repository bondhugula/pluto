

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
/* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 2.33s. */
for (c1=-1;c1<=floord(3*N-5,128);c1++) {
	lb1=max(max(ceild(64*c1-N+2,64),ceild(32*c1-63,96)),0);
	ub1=min(floord(64*c1+63,64),floord(N-1,128));
#pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5,c6,c7,c8,c9)
	for (c2=lb1; c2<=ub1; c2++) {
    for (c3=max(ceild(32*c1-32*c2-1953,2016),ceild(32*c1-32*c2-31,32));c3<=floord(N-1,64);c3++) {
      for (c4=max(max(0,2*c1-2*c2-64*c3-62),2*c1-2*c2);c4<=min(min(min(min(2*c1-2*c2+1,floord(992*c3+961,16)),floord(N-2,32)),floord(64*c2+63,16)),floord(32*c3+31,16));c4++) {
        for (c5=max(max(ceild(16*c4-7,8),0),8*c2);c5<=min(8*c2+7,floord(N-1,16));c5++) {
          for (c6=max(max(max(max(ceild(16*c4-465,496),ceild(2*c1-2*c2-2*c3-c4-31,31)),ceild(-2*c1+2*c2+2*c3+c4-31,33)),2*c3),ceild(16*c4-15,16));c6<=min(2*c3+1,floord(N-1,32));c6++) {
            if ((c1 == c2+c3) && (c4 == c6)) {
              for (c7=max(0,32*c6);c7<=min(min(32*c6+30,N-2),16*c5+14);c7++) {
                for (c8=max(16*c5,c7+1);c8<=min(16*c5+15,N-1);c8++) {
                  A[c7][c8]=A[c7][c8]/A[c7][c7] ;
                  for (c9=c7+1;c9<=min(32*c6+31,N-1);c9++) {
                    A[c9][c8]=A[c9][c8]-A[c9][c7]*A[c7][c8] ;
                  }
                }
              }
            }
            
/*@ begin Loop(
transform Composite(
  regtile = (['c7', 'c8', 'c9'],[32, 8, 1]),
  permut = [(['c7'],['c8'],['c9'])],
  scalarreplace = (True, 'double'),
  vector = (True, ['ivdep','vector always']))
for (c7=max(32*c4,0);c7<=min(min(32*c6-1,16*c5+14),32*c4+31);c7++) {

              for (c8=max(c7+1,16*c5);c8<=min(16*c5+15,N-1);c8++) 
                for (c9=32*c6;c9<=min(N-1,32*c6+31);c9++) 
{
                  A[c9][c8]=A[c9][c8]-A[c9][c7]*A[c7][c8] ;
}

            }
) @*/{
  for (c7t=max(32*c4,0); c7t<=min(min(32*c6-1,16*c5+14),32*c4+31)-31; c7t=c7t+32) {
    newlb_c8=-2147483648;
    newub_c8=min(16*c5+15,N-1);
    register int cbv_1;
    cbv_1=c7t+31;
#pragma ivdep
#pragma vector always
    for (c7=c7t; c7<=cbv_1; c7=c7+1) {
      newlb_c8=max(newlb_c8,max(c7+1,16*c5));
    }
    for (c7=c7t; c7<=c7t+31; c7=c7+1) {
      for (c8=max(c7+1,16*c5); c8<=newlb_c8-1; c8=c8+1) {
        register int cbv_2, cbv_3;
        cbv_2=32*c6;
        cbv_3=min(N-1,32*c6+31);
#pragma ivdep
#pragma vector always
        for (c9=cbv_2; c9<=cbv_3; c9++ ) {
          double scv_1;
          scv_1=A[c9][c8];
          scv_1=scv_1-A[c9][c7]*A[c7][c8];
          A[c9][c8]=scv_1;
        }
      }
    }
    for (c8t=newlb_c8; c8t<=newub_c8-7; c8t=c8t+8) {
      register int cbv_4, cbv_5;
      cbv_4=32*c6;
      cbv_5=min(N-1,32*c6+31);
#pragma ivdep
#pragma vector always
      for (c9=cbv_4; c9<=cbv_5; c9++ ) {
        double scv_2, scv_3, scv_4, scv_5, scv_6, scv_7, scv_8, scv_9;
        double scv_10, scv_11, scv_12, scv_13, scv_14, scv_15, scv_16, scv_17;
        double scv_18, scv_19, scv_20, scv_21, scv_22, scv_23, scv_24, scv_25;
        double scv_26, scv_27, scv_28, scv_29, scv_30, scv_31, scv_32, scv_33;
        double scv_34, scv_35, scv_36, scv_37, scv_38, scv_39, scv_40, scv_41;
        scv_2=A[c9][(c7t+31)];
        scv_3=A[c9][c7t];
        scv_4=A[c9][(c7t+30)];
        scv_5=A[c9][(c8t+7)];
        scv_6=A[c9][(c7t+20)];
        scv_7=A[c9][(c8t+5)];
        scv_8=A[c9][(c7t+19)];
        scv_9=A[c9][(c7t+27)];
        scv_10=A[c9][(c7t+2)];
        scv_11=A[c9][(c8t+1)];
        scv_12=A[c9][(c7t+22)];
        scv_13=A[c9][(c7t+5)];
        scv_14=A[c9][(c7t+11)];
        scv_15=A[c9][(c7t+6)];
        scv_16=A[c9][(c7t+15)];
        scv_17=A[c9][(c7t+1)];
        scv_18=A[c9][(c7t+21)];
        scv_19=A[c9][(c7t+9)];
        scv_20=A[c9][(c7t+16)];
        scv_21=A[c9][(c7t+12)];
        scv_22=A[c9][(c7t+17)];
        scv_23=A[c9][(c8t+2)];
        scv_24=A[c9][c8t];
        scv_25=A[c9][(c7t+23)];
        scv_26=A[c9][(c7t+3)];
        scv_27=A[c9][(c8t+4)];
        scv_28=A[c9][(c7t+4)];
        scv_29=A[c9][(c7t+25)];
        scv_30=A[c9][(c8t+6)];
        scv_31=A[c9][(c7t+18)];
        scv_32=A[c9][(c7t+7)];
        scv_33=A[c9][(c7t+14)];
        scv_34=A[c9][(c7t+26)];
        scv_35=A[c9][(c7t+8)];
        scv_36=A[c9][(c7t+24)];
        scv_37=A[c9][(c7t+29)];
        scv_38=A[c9][(c7t+28)];
        scv_39=A[c9][(c7t+13)];
        scv_40=A[c9][(c8t+3)];
        scv_41=A[c9][(c7t+10)];
        scv_24=scv_24-scv_3*A[c7t][c8t];
        scv_11=scv_11-scv_3*A[c7t][(c8t+1)];
        scv_23=scv_23-scv_3*A[c7t][(c8t+2)];
        scv_40=scv_40-scv_3*A[c7t][(c8t+3)];
        scv_27=scv_27-scv_3*A[c7t][(c8t+4)];
        scv_7=scv_7-scv_3*A[c7t][(c8t+5)];
        scv_30=scv_30-scv_3*A[c7t][(c8t+6)];
        scv_5=scv_5-scv_3*A[c7t][(c8t+7)];
        scv_24=scv_24-scv_17*A[(c7t+1)][c8t];
        scv_11=scv_11-scv_17*A[(c7t+1)][(c8t+1)];
        scv_23=scv_23-scv_17*A[(c7t+1)][(c8t+2)];
        scv_40=scv_40-scv_17*A[(c7t+1)][(c8t+3)];
        scv_27=scv_27-scv_17*A[(c7t+1)][(c8t+4)];
        scv_7=scv_7-scv_17*A[(c7t+1)][(c8t+5)];
        scv_30=scv_30-scv_17*A[(c7t+1)][(c8t+6)];
        scv_5=scv_5-scv_17*A[(c7t+1)][(c8t+7)];
        scv_24=scv_24-scv_10*A[(c7t+2)][c8t];
        scv_11=scv_11-scv_10*A[(c7t+2)][(c8t+1)];
        scv_23=scv_23-scv_10*A[(c7t+2)][(c8t+2)];
        scv_40=scv_40-scv_10*A[(c7t+2)][(c8t+3)];
        scv_27=scv_27-scv_10*A[(c7t+2)][(c8t+4)];
        scv_7=scv_7-scv_10*A[(c7t+2)][(c8t+5)];
        scv_30=scv_30-scv_10*A[(c7t+2)][(c8t+6)];
        scv_5=scv_5-scv_10*A[(c7t+2)][(c8t+7)];
        scv_24=scv_24-scv_26*A[(c7t+3)][c8t];
        scv_11=scv_11-scv_26*A[(c7t+3)][(c8t+1)];
        scv_23=scv_23-scv_26*A[(c7t+3)][(c8t+2)];
        scv_40=scv_40-scv_26*A[(c7t+3)][(c8t+3)];
        scv_27=scv_27-scv_26*A[(c7t+3)][(c8t+4)];
        scv_7=scv_7-scv_26*A[(c7t+3)][(c8t+5)];
        scv_30=scv_30-scv_26*A[(c7t+3)][(c8t+6)];
        scv_5=scv_5-scv_26*A[(c7t+3)][(c8t+7)];
        scv_24=scv_24-scv_28*A[(c7t+4)][c8t];
        scv_11=scv_11-scv_28*A[(c7t+4)][(c8t+1)];
        scv_23=scv_23-scv_28*A[(c7t+4)][(c8t+2)];
        scv_40=scv_40-scv_28*A[(c7t+4)][(c8t+3)];
        scv_27=scv_27-scv_28*A[(c7t+4)][(c8t+4)];
        scv_7=scv_7-scv_28*A[(c7t+4)][(c8t+5)];
        scv_30=scv_30-scv_28*A[(c7t+4)][(c8t+6)];
        scv_5=scv_5-scv_28*A[(c7t+4)][(c8t+7)];
        scv_24=scv_24-scv_13*A[(c7t+5)][c8t];
        scv_11=scv_11-scv_13*A[(c7t+5)][(c8t+1)];
        scv_23=scv_23-scv_13*A[(c7t+5)][(c8t+2)];
        scv_40=scv_40-scv_13*A[(c7t+5)][(c8t+3)];
        scv_27=scv_27-scv_13*A[(c7t+5)][(c8t+4)];
        scv_7=scv_7-scv_13*A[(c7t+5)][(c8t+5)];
        scv_30=scv_30-scv_13*A[(c7t+5)][(c8t+6)];
        scv_5=scv_5-scv_13*A[(c7t+5)][(c8t+7)];
        scv_24=scv_24-scv_15*A[(c7t+6)][c8t];
        scv_11=scv_11-scv_15*A[(c7t+6)][(c8t+1)];
        scv_23=scv_23-scv_15*A[(c7t+6)][(c8t+2)];
        scv_40=scv_40-scv_15*A[(c7t+6)][(c8t+3)];
        scv_27=scv_27-scv_15*A[(c7t+6)][(c8t+4)];
        scv_7=scv_7-scv_15*A[(c7t+6)][(c8t+5)];
        scv_30=scv_30-scv_15*A[(c7t+6)][(c8t+6)];
        scv_5=scv_5-scv_15*A[(c7t+6)][(c8t+7)];
        scv_24=scv_24-scv_32*A[(c7t+7)][c8t];
        scv_11=scv_11-scv_32*A[(c7t+7)][(c8t+1)];
        scv_23=scv_23-scv_32*A[(c7t+7)][(c8t+2)];
        scv_40=scv_40-scv_32*A[(c7t+7)][(c8t+3)];
        scv_27=scv_27-scv_32*A[(c7t+7)][(c8t+4)];
        scv_7=scv_7-scv_32*A[(c7t+7)][(c8t+5)];
        scv_30=scv_30-scv_32*A[(c7t+7)][(c8t+6)];
        scv_5=scv_5-scv_32*A[(c7t+7)][(c8t+7)];
        scv_24=scv_24-scv_35*A[(c7t+8)][c8t];
        scv_11=scv_11-scv_35*A[(c7t+8)][(c8t+1)];
        scv_23=scv_23-scv_35*A[(c7t+8)][(c8t+2)];
        scv_40=scv_40-scv_35*A[(c7t+8)][(c8t+3)];
        scv_27=scv_27-scv_35*A[(c7t+8)][(c8t+4)];
        scv_7=scv_7-scv_35*A[(c7t+8)][(c8t+5)];
        scv_30=scv_30-scv_35*A[(c7t+8)][(c8t+6)];
        scv_5=scv_5-scv_35*A[(c7t+8)][(c8t+7)];
        scv_24=scv_24-scv_19*A[(c7t+9)][c8t];
        scv_11=scv_11-scv_19*A[(c7t+9)][(c8t+1)];
        scv_23=scv_23-scv_19*A[(c7t+9)][(c8t+2)];
        scv_40=scv_40-scv_19*A[(c7t+9)][(c8t+3)];
        scv_27=scv_27-scv_19*A[(c7t+9)][(c8t+4)];
        scv_7=scv_7-scv_19*A[(c7t+9)][(c8t+5)];
        scv_30=scv_30-scv_19*A[(c7t+9)][(c8t+6)];
        scv_5=scv_5-scv_19*A[(c7t+9)][(c8t+7)];
        scv_24=scv_24-scv_41*A[(c7t+10)][c8t];
        scv_11=scv_11-scv_41*A[(c7t+10)][(c8t+1)];
        scv_23=scv_23-scv_41*A[(c7t+10)][(c8t+2)];
        scv_40=scv_40-scv_41*A[(c7t+10)][(c8t+3)];
        scv_27=scv_27-scv_41*A[(c7t+10)][(c8t+4)];
        scv_7=scv_7-scv_41*A[(c7t+10)][(c8t+5)];
        scv_30=scv_30-scv_41*A[(c7t+10)][(c8t+6)];
        scv_5=scv_5-scv_41*A[(c7t+10)][(c8t+7)];
        scv_24=scv_24-scv_14*A[(c7t+11)][c8t];
        scv_11=scv_11-scv_14*A[(c7t+11)][(c8t+1)];
        scv_23=scv_23-scv_14*A[(c7t+11)][(c8t+2)];
        scv_40=scv_40-scv_14*A[(c7t+11)][(c8t+3)];
        scv_27=scv_27-scv_14*A[(c7t+11)][(c8t+4)];
        scv_7=scv_7-scv_14*A[(c7t+11)][(c8t+5)];
        scv_30=scv_30-scv_14*A[(c7t+11)][(c8t+6)];
        scv_5=scv_5-scv_14*A[(c7t+11)][(c8t+7)];
        scv_24=scv_24-scv_21*A[(c7t+12)][c8t];
        scv_11=scv_11-scv_21*A[(c7t+12)][(c8t+1)];
        scv_23=scv_23-scv_21*A[(c7t+12)][(c8t+2)];
        scv_40=scv_40-scv_21*A[(c7t+12)][(c8t+3)];
        scv_27=scv_27-scv_21*A[(c7t+12)][(c8t+4)];
        scv_7=scv_7-scv_21*A[(c7t+12)][(c8t+5)];
        scv_30=scv_30-scv_21*A[(c7t+12)][(c8t+6)];
        scv_5=scv_5-scv_21*A[(c7t+12)][(c8t+7)];
        scv_24=scv_24-scv_39*A[(c7t+13)][c8t];
        scv_11=scv_11-scv_39*A[(c7t+13)][(c8t+1)];
        scv_23=scv_23-scv_39*A[(c7t+13)][(c8t+2)];
        scv_40=scv_40-scv_39*A[(c7t+13)][(c8t+3)];
        scv_27=scv_27-scv_39*A[(c7t+13)][(c8t+4)];
        scv_7=scv_7-scv_39*A[(c7t+13)][(c8t+5)];
        scv_30=scv_30-scv_39*A[(c7t+13)][(c8t+6)];
        scv_5=scv_5-scv_39*A[(c7t+13)][(c8t+7)];
        scv_24=scv_24-scv_33*A[(c7t+14)][c8t];
        scv_11=scv_11-scv_33*A[(c7t+14)][(c8t+1)];
        scv_23=scv_23-scv_33*A[(c7t+14)][(c8t+2)];
        scv_40=scv_40-scv_33*A[(c7t+14)][(c8t+3)];
        scv_27=scv_27-scv_33*A[(c7t+14)][(c8t+4)];
        scv_7=scv_7-scv_33*A[(c7t+14)][(c8t+5)];
        scv_30=scv_30-scv_33*A[(c7t+14)][(c8t+6)];
        scv_5=scv_5-scv_33*A[(c7t+14)][(c8t+7)];
        scv_24=scv_24-scv_16*A[(c7t+15)][c8t];
        scv_11=scv_11-scv_16*A[(c7t+15)][(c8t+1)];
        scv_23=scv_23-scv_16*A[(c7t+15)][(c8t+2)];
        scv_40=scv_40-scv_16*A[(c7t+15)][(c8t+3)];
        scv_27=scv_27-scv_16*A[(c7t+15)][(c8t+4)];
        scv_7=scv_7-scv_16*A[(c7t+15)][(c8t+5)];
        scv_30=scv_30-scv_16*A[(c7t+15)][(c8t+6)];
        scv_5=scv_5-scv_16*A[(c7t+15)][(c8t+7)];
        scv_24=scv_24-scv_20*A[(c7t+16)][c8t];
        scv_11=scv_11-scv_20*A[(c7t+16)][(c8t+1)];
        scv_23=scv_23-scv_20*A[(c7t+16)][(c8t+2)];
        scv_40=scv_40-scv_20*A[(c7t+16)][(c8t+3)];
        scv_27=scv_27-scv_20*A[(c7t+16)][(c8t+4)];
        scv_7=scv_7-scv_20*A[(c7t+16)][(c8t+5)];
        scv_30=scv_30-scv_20*A[(c7t+16)][(c8t+6)];
        scv_5=scv_5-scv_20*A[(c7t+16)][(c8t+7)];
        scv_24=scv_24-scv_22*A[(c7t+17)][c8t];
        scv_11=scv_11-scv_22*A[(c7t+17)][(c8t+1)];
        scv_23=scv_23-scv_22*A[(c7t+17)][(c8t+2)];
        scv_40=scv_40-scv_22*A[(c7t+17)][(c8t+3)];
        scv_27=scv_27-scv_22*A[(c7t+17)][(c8t+4)];
        scv_7=scv_7-scv_22*A[(c7t+17)][(c8t+5)];
        scv_30=scv_30-scv_22*A[(c7t+17)][(c8t+6)];
        scv_5=scv_5-scv_22*A[(c7t+17)][(c8t+7)];
        scv_24=scv_24-scv_31*A[(c7t+18)][c8t];
        scv_11=scv_11-scv_31*A[(c7t+18)][(c8t+1)];
        scv_23=scv_23-scv_31*A[(c7t+18)][(c8t+2)];
        scv_40=scv_40-scv_31*A[(c7t+18)][(c8t+3)];
        scv_27=scv_27-scv_31*A[(c7t+18)][(c8t+4)];
        scv_7=scv_7-scv_31*A[(c7t+18)][(c8t+5)];
        scv_30=scv_30-scv_31*A[(c7t+18)][(c8t+6)];
        scv_5=scv_5-scv_31*A[(c7t+18)][(c8t+7)];
        scv_24=scv_24-scv_8*A[(c7t+19)][c8t];
        scv_11=scv_11-scv_8*A[(c7t+19)][(c8t+1)];
        scv_23=scv_23-scv_8*A[(c7t+19)][(c8t+2)];
        scv_40=scv_40-scv_8*A[(c7t+19)][(c8t+3)];
        scv_27=scv_27-scv_8*A[(c7t+19)][(c8t+4)];
        scv_7=scv_7-scv_8*A[(c7t+19)][(c8t+5)];
        scv_30=scv_30-scv_8*A[(c7t+19)][(c8t+6)];
        scv_5=scv_5-scv_8*A[(c7t+19)][(c8t+7)];
        scv_24=scv_24-scv_6*A[(c7t+20)][c8t];
        scv_11=scv_11-scv_6*A[(c7t+20)][(c8t+1)];
        scv_23=scv_23-scv_6*A[(c7t+20)][(c8t+2)];
        scv_40=scv_40-scv_6*A[(c7t+20)][(c8t+3)];
        scv_27=scv_27-scv_6*A[(c7t+20)][(c8t+4)];
        scv_7=scv_7-scv_6*A[(c7t+20)][(c8t+5)];
        scv_30=scv_30-scv_6*A[(c7t+20)][(c8t+6)];
        scv_5=scv_5-scv_6*A[(c7t+20)][(c8t+7)];
        scv_24=scv_24-scv_18*A[(c7t+21)][c8t];
        scv_11=scv_11-scv_18*A[(c7t+21)][(c8t+1)];
        scv_23=scv_23-scv_18*A[(c7t+21)][(c8t+2)];
        scv_40=scv_40-scv_18*A[(c7t+21)][(c8t+3)];
        scv_27=scv_27-scv_18*A[(c7t+21)][(c8t+4)];
        scv_7=scv_7-scv_18*A[(c7t+21)][(c8t+5)];
        scv_30=scv_30-scv_18*A[(c7t+21)][(c8t+6)];
        scv_5=scv_5-scv_18*A[(c7t+21)][(c8t+7)];
        scv_24=scv_24-scv_12*A[(c7t+22)][c8t];
        scv_11=scv_11-scv_12*A[(c7t+22)][(c8t+1)];
        scv_23=scv_23-scv_12*A[(c7t+22)][(c8t+2)];
        scv_40=scv_40-scv_12*A[(c7t+22)][(c8t+3)];
        scv_27=scv_27-scv_12*A[(c7t+22)][(c8t+4)];
        scv_7=scv_7-scv_12*A[(c7t+22)][(c8t+5)];
        scv_30=scv_30-scv_12*A[(c7t+22)][(c8t+6)];
        scv_5=scv_5-scv_12*A[(c7t+22)][(c8t+7)];
        scv_24=scv_24-scv_25*A[(c7t+23)][c8t];
        scv_11=scv_11-scv_25*A[(c7t+23)][(c8t+1)];
        scv_23=scv_23-scv_25*A[(c7t+23)][(c8t+2)];
        scv_40=scv_40-scv_25*A[(c7t+23)][(c8t+3)];
        scv_27=scv_27-scv_25*A[(c7t+23)][(c8t+4)];
        scv_7=scv_7-scv_25*A[(c7t+23)][(c8t+5)];
        scv_30=scv_30-scv_25*A[(c7t+23)][(c8t+6)];
        scv_5=scv_5-scv_25*A[(c7t+23)][(c8t+7)];
        scv_24=scv_24-scv_36*A[(c7t+24)][c8t];
        scv_11=scv_11-scv_36*A[(c7t+24)][(c8t+1)];
        scv_23=scv_23-scv_36*A[(c7t+24)][(c8t+2)];
        scv_40=scv_40-scv_36*A[(c7t+24)][(c8t+3)];
        scv_27=scv_27-scv_36*A[(c7t+24)][(c8t+4)];
        scv_7=scv_7-scv_36*A[(c7t+24)][(c8t+5)];
        scv_30=scv_30-scv_36*A[(c7t+24)][(c8t+6)];
        scv_5=scv_5-scv_36*A[(c7t+24)][(c8t+7)];
        scv_24=scv_24-scv_29*A[(c7t+25)][c8t];
        scv_11=scv_11-scv_29*A[(c7t+25)][(c8t+1)];
        scv_23=scv_23-scv_29*A[(c7t+25)][(c8t+2)];
        scv_40=scv_40-scv_29*A[(c7t+25)][(c8t+3)];
        scv_27=scv_27-scv_29*A[(c7t+25)][(c8t+4)];
        scv_7=scv_7-scv_29*A[(c7t+25)][(c8t+5)];
        scv_30=scv_30-scv_29*A[(c7t+25)][(c8t+6)];
        scv_5=scv_5-scv_29*A[(c7t+25)][(c8t+7)];
        scv_24=scv_24-scv_34*A[(c7t+26)][c8t];
        scv_11=scv_11-scv_34*A[(c7t+26)][(c8t+1)];
        scv_23=scv_23-scv_34*A[(c7t+26)][(c8t+2)];
        scv_40=scv_40-scv_34*A[(c7t+26)][(c8t+3)];
        scv_27=scv_27-scv_34*A[(c7t+26)][(c8t+4)];
        scv_7=scv_7-scv_34*A[(c7t+26)][(c8t+5)];
        scv_30=scv_30-scv_34*A[(c7t+26)][(c8t+6)];
        scv_5=scv_5-scv_34*A[(c7t+26)][(c8t+7)];
        scv_24=scv_24-scv_9*A[(c7t+27)][c8t];
        scv_11=scv_11-scv_9*A[(c7t+27)][(c8t+1)];
        scv_23=scv_23-scv_9*A[(c7t+27)][(c8t+2)];
        scv_40=scv_40-scv_9*A[(c7t+27)][(c8t+3)];
        scv_27=scv_27-scv_9*A[(c7t+27)][(c8t+4)];
        scv_7=scv_7-scv_9*A[(c7t+27)][(c8t+5)];
        scv_30=scv_30-scv_9*A[(c7t+27)][(c8t+6)];
        scv_5=scv_5-scv_9*A[(c7t+27)][(c8t+7)];
        scv_24=scv_24-scv_38*A[(c7t+28)][c8t];
        scv_11=scv_11-scv_38*A[(c7t+28)][(c8t+1)];
        scv_23=scv_23-scv_38*A[(c7t+28)][(c8t+2)];
        scv_40=scv_40-scv_38*A[(c7t+28)][(c8t+3)];
        scv_27=scv_27-scv_38*A[(c7t+28)][(c8t+4)];
        scv_7=scv_7-scv_38*A[(c7t+28)][(c8t+5)];
        scv_30=scv_30-scv_38*A[(c7t+28)][(c8t+6)];
        scv_5=scv_5-scv_38*A[(c7t+28)][(c8t+7)];
        scv_24=scv_24-scv_37*A[(c7t+29)][c8t];
        scv_11=scv_11-scv_37*A[(c7t+29)][(c8t+1)];
        scv_23=scv_23-scv_37*A[(c7t+29)][(c8t+2)];
        scv_40=scv_40-scv_37*A[(c7t+29)][(c8t+3)];
        scv_27=scv_27-scv_37*A[(c7t+29)][(c8t+4)];
        scv_7=scv_7-scv_37*A[(c7t+29)][(c8t+5)];
        scv_30=scv_30-scv_37*A[(c7t+29)][(c8t+6)];
        scv_5=scv_5-scv_37*A[(c7t+29)][(c8t+7)];
        scv_24=scv_24-scv_4*A[(c7t+30)][c8t];
        scv_11=scv_11-scv_4*A[(c7t+30)][(c8t+1)];
        scv_23=scv_23-scv_4*A[(c7t+30)][(c8t+2)];
        scv_40=scv_40-scv_4*A[(c7t+30)][(c8t+3)];
        scv_27=scv_27-scv_4*A[(c7t+30)][(c8t+4)];
        scv_7=scv_7-scv_4*A[(c7t+30)][(c8t+5)];
        scv_30=scv_30-scv_4*A[(c7t+30)][(c8t+6)];
        scv_5=scv_5-scv_4*A[(c7t+30)][(c8t+7)];
        scv_24=scv_24-scv_2*A[(c7t+31)][c8t];
        scv_11=scv_11-scv_2*A[(c7t+31)][(c8t+1)];
        scv_23=scv_23-scv_2*A[(c7t+31)][(c8t+2)];
        scv_40=scv_40-scv_2*A[(c7t+31)][(c8t+3)];
        scv_27=scv_27-scv_2*A[(c7t+31)][(c8t+4)];
        scv_7=scv_7-scv_2*A[(c7t+31)][(c8t+5)];
        scv_30=scv_30-scv_2*A[(c7t+31)][(c8t+6)];
        scv_5=scv_5-scv_2*A[(c7t+31)][(c8t+7)];
        A[c9][(c8t+7)]=scv_5;
        A[c9][(c8t+5)]=scv_7;
        A[c9][(c8t+1)]=scv_11;
        A[c9][(c8t+2)]=scv_23;
        A[c9][c8t]=scv_24;
        A[c9][(c8t+4)]=scv_27;
        A[c9][(c8t+6)]=scv_30;
        A[c9][(c8t+3)]=scv_40;
      }
    }
    for (c8=c8t; c8<=newub_c8; c8=c8+1) {
      register int cbv_6, cbv_7;
      cbv_6=32*c6;
      cbv_7=min(N-1,32*c6+31);
#pragma ivdep
#pragma vector always
      for (c9=cbv_6; c9<=cbv_7; c9++ ) {
        double scv_42;
        scv_42=A[c9][c8];
        scv_42=scv_42-A[c9][c7t]*A[c7t][c8];
        scv_42=scv_42-A[c9][(c7t+1)]*A[(c7t+1)][c8];
        scv_42=scv_42-A[c9][(c7t+2)]*A[(c7t+2)][c8];
        scv_42=scv_42-A[c9][(c7t+3)]*A[(c7t+3)][c8];
        scv_42=scv_42-A[c9][(c7t+4)]*A[(c7t+4)][c8];
        scv_42=scv_42-A[c9][(c7t+5)]*A[(c7t+5)][c8];
        scv_42=scv_42-A[c9][(c7t+6)]*A[(c7t+6)][c8];
        scv_42=scv_42-A[c9][(c7t+7)]*A[(c7t+7)][c8];
        scv_42=scv_42-A[c9][(c7t+8)]*A[(c7t+8)][c8];
        scv_42=scv_42-A[c9][(c7t+9)]*A[(c7t+9)][c8];
        scv_42=scv_42-A[c9][(c7t+10)]*A[(c7t+10)][c8];
        scv_42=scv_42-A[c9][(c7t+11)]*A[(c7t+11)][c8];
        scv_42=scv_42-A[c9][(c7t+12)]*A[(c7t+12)][c8];
        scv_42=scv_42-A[c9][(c7t+13)]*A[(c7t+13)][c8];
        scv_42=scv_42-A[c9][(c7t+14)]*A[(c7t+14)][c8];
        scv_42=scv_42-A[c9][(c7t+15)]*A[(c7t+15)][c8];
        scv_42=scv_42-A[c9][(c7t+16)]*A[(c7t+16)][c8];
        scv_42=scv_42-A[c9][(c7t+17)]*A[(c7t+17)][c8];
        scv_42=scv_42-A[c9][(c7t+18)]*A[(c7t+18)][c8];
        scv_42=scv_42-A[c9][(c7t+19)]*A[(c7t+19)][c8];
        scv_42=scv_42-A[c9][(c7t+20)]*A[(c7t+20)][c8];
        scv_42=scv_42-A[c9][(c7t+21)]*A[(c7t+21)][c8];
        scv_42=scv_42-A[c9][(c7t+22)]*A[(c7t+22)][c8];
        scv_42=scv_42-A[c9][(c7t+23)]*A[(c7t+23)][c8];
        scv_42=scv_42-A[c9][(c7t+24)]*A[(c7t+24)][c8];
        scv_42=scv_42-A[c9][(c7t+25)]*A[(c7t+25)][c8];
        scv_42=scv_42-A[c9][(c7t+26)]*A[(c7t+26)][c8];
        scv_42=scv_42-A[c9][(c7t+27)]*A[(c7t+27)][c8];
        scv_42=scv_42-A[c9][(c7t+28)]*A[(c7t+28)][c8];
        scv_42=scv_42-A[c9][(c7t+29)]*A[(c7t+29)][c8];
        scv_42=scv_42-A[c9][(c7t+30)]*A[(c7t+30)][c8];
        scv_42=scv_42-A[c9][(c7t+31)]*A[(c7t+31)][c8];
        A[c9][c8]=scv_42;
      }
    }
    for (c7=c7t; c7<=c7t+31; c7=c7+1) {
      for (c8=newub_c8+1; c8<=min(16*c5+15,N-1); c8=c8+1) {
        register int cbv_8, cbv_9;
        cbv_8=32*c6;
        cbv_9=min(N-1,32*c6+31);
#pragma ivdep
#pragma vector always
        for (c9=cbv_8; c9<=cbv_9; c9++ ) {
          double scv_43;
          scv_43=A[c9][c8];
          scv_43=scv_43-A[c9][c7]*A[c7][c8];
          A[c9][c8]=scv_43;
        }
      }
    }
  }
  for (c7=c7t; c7<=min(min(32*c6-1,16*c5+14),32*c4+31); c7=c7+1) {
    for (c8t=max(c7+1,16*c5); c8t<=min(16*c5+15,N-1)-7; c8t=c8t+8) {
      register int cbv_10, cbv_11;
      cbv_10=32*c6;
      cbv_11=min(N-1,32*c6+31);
#pragma ivdep
#pragma vector always
      for (c9=cbv_10; c9<=cbv_11; c9++ ) {
        double scv_44, scv_45, scv_46, scv_47, scv_48, scv_49, scv_50, scv_51;
        double scv_52;
        scv_44=A[c9][(c8t+6)];
        scv_45=A[c9][(c8t+4)];
        scv_46=A[c9][(c8t+5)];
        scv_47=A[c9][(c8t+2)];
        scv_48=A[c9][c8t];
        scv_49=A[c9][(c8t+3)];
        scv_50=A[c9][c7];
        scv_51=A[c9][(c8t+7)];
        scv_52=A[c9][(c8t+1)];
        scv_48=scv_48-scv_50*A[c7][c8t];
        scv_52=scv_52-scv_50*A[c7][(c8t+1)];
        scv_47=scv_47-scv_50*A[c7][(c8t+2)];
        scv_49=scv_49-scv_50*A[c7][(c8t+3)];
        scv_45=scv_45-scv_50*A[c7][(c8t+4)];
        scv_46=scv_46-scv_50*A[c7][(c8t+5)];
        scv_44=scv_44-scv_50*A[c7][(c8t+6)];
        scv_51=scv_51-scv_50*A[c7][(c8t+7)];
        A[c9][(c8t+6)]=scv_44;
        A[c9][(c8t+4)]=scv_45;
        A[c9][(c8t+5)]=scv_46;
        A[c9][(c8t+2)]=scv_47;
        A[c9][c8t]=scv_48;
        A[c9][(c8t+3)]=scv_49;
        A[c9][(c8t+7)]=scv_51;
        A[c9][(c8t+1)]=scv_52;
      }
    }
    for (c8=c8t; c8<=min(16*c5+15,N-1); c8=c8+1) {
      register int cbv_12, cbv_13;
      cbv_12=32*c6;
      cbv_13=min(N-1,32*c6+31);
#pragma ivdep
#pragma vector always
      for (c9=cbv_12; c9<=cbv_13; c9++ ) {
        double scv_53;
        scv_53=A[c9][c8];
        scv_53=scv_53-A[c9][c7]*A[c7][c8];
        A[c9][c8]=scv_53;
      }
    }
  }
}
/*@ end @*/


            if ((c1 == c2+c3) && (-c4 == -c6) && (c4 <= min(floord(N-33,32),floord(16*c5-17,32)))) {
              for (c8=max(16*c5,32*c4+32);c8<=min(N-1,16*c5+15);c8++) {
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
                                    
