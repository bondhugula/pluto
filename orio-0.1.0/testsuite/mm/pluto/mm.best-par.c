

#include <stdio.h>
#include <sys/time.h>
#include <math.h>


#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define M NCONT
#define N NCONT
#define K CONT
double A[M][K];
double B[K][N];
double C[M][N];

void init_arrays()
{
  int i1, i2;
  for (i1=0; i1<M; i1++)
   for (i2=0; i2<K; i2++)
    A[i1][i2] = (i1+i2) % 5 + 1;
  for (i1=0; i1<K; i1++)
   for (i2=0; i2<N; i2++)
    B[i1][i2] = (i1+i2) % 5 + 1;
  for (i1=0; i1<M; i1++)
   for (i2=0; i2<N; i2++)
    C[i1][i2] = 0;
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
    
  

int i, j, k;
int ii, jj, kk;
int iii, jjj, kkk;

{
  double C_copy[128][64];
  double B_copy;
  double A_copy[256][128];
  register int cbv_1;
  cbv_1=K-1;
#pragma omp parallel for private(iii,jjj,kkk,ii,jj,kk,i,j,k,A_copy,B_copy,C_copy)
  for (kkk=0; kkk<=cbv_1; kkk=kkk+512) {
    for (iii=0; iii<=M-1; iii=iii+256) {
      for (jjj=0; jjj<=N-1; jjj=jjj+1024) {
        for (kk=kkk; kk<=min(K-1,kkk+256); kk=kk+256) {
          for (ii=iii; ii<=min(M-1,iii+128); ii=ii+128) {
            for (k=kk; k<=min(K-1,kk+255); k=k+1) 
              for (i=ii; i<=min(M-1,ii+127); i=i+1) 
                A_copy[(k-kk)][(i-ii)]=A[i][k];
            for (jj=jjj; jj<=min(N-1,jjj+960); jj=jj+64) {
              for (i=ii; i<=min(M-1,ii+127); i=i+1) 
                for (j=jj; j<=min(N-1,jj+63); j=j+1) 
                  C_copy[(i-ii)][(j-jj)]=C[i][j];
              for (k=kk; k<=min(K-1,kk+255)-7; k=k+8) {
                for (i=ii; i<=min(M-1,ii+127)-7; i=i+8) {
                  register int cbv_2;
                  cbv_2=min(N-1,jj+63);
#pragma ivdep
#pragma vector always
                  for (j=jj; j<=cbv_2; j=j+1) {
                    double scv_1, scv_2, scv_3, scv_4, scv_5, scv_6, scv_7, scv_8;
                    double scv_9, scv_10, scv_11, scv_12, scv_13, scv_14, scv_15, scv_16;
                    scv_1=B[k][j];
                    scv_2=B[(k+6)][j];
                    scv_3=B[(k+5)][j];
                    scv_4=C_copy[(i-ii+4)][(j-jj)];
                    scv_5=C_copy[(i-ii+2)][(j-jj)];
                    scv_6=B[(k+4)][j];
                    scv_7=C_copy[(i-ii+3)][(j-jj)];
                    scv_8=C_copy[(i-ii+6)][(j-jj)];
                    scv_9=B[(k+3)][j];
                    scv_10=C_copy[(i-ii+5)][(j-jj)];
                    scv_11=C_copy[(i-ii+1)][(j-jj)];
                    scv_12=B[(k+1)][j];
                    scv_13=C_copy[(i-ii+7)][(j-jj)];
                    scv_14=B[(k+2)][j];
                    scv_15=B[(k+7)][j];
                    scv_16=C_copy[(i-ii)][(j-jj)];
                    scv_16=scv_16+A_copy[(k-kk)][(i-ii)]*scv_1;
                    scv_11=scv_11+A_copy[(k-kk)][(i-ii+1)]*scv_1;
                    scv_5=scv_5+A_copy[(k-kk)][(i-ii+2)]*scv_1;
                    scv_7=scv_7+A_copy[(k-kk)][(i-ii+3)]*scv_1;
                    scv_4=scv_4+A_copy[(k-kk)][(i-ii+4)]*scv_1;
                    scv_10=scv_10+A_copy[(k-kk)][(i-ii+5)]*scv_1;
                    scv_8=scv_8+A_copy[(k-kk)][(i-ii+6)]*scv_1;
                    scv_13=scv_13+A_copy[(k-kk)][(i-ii+7)]*scv_1;
                    scv_16=scv_16+A_copy[(k-kk+1)][(i-ii)]*scv_12;
                    scv_11=scv_11+A_copy[(k-kk+1)][(i-ii+1)]*scv_12;
                    scv_5=scv_5+A_copy[(k-kk+1)][(i-ii+2)]*scv_12;
                    scv_7=scv_7+A_copy[(k-kk+1)][(i-ii+3)]*scv_12;
                    scv_4=scv_4+A_copy[(k-kk+1)][(i-ii+4)]*scv_12;
                    scv_10=scv_10+A_copy[(k-kk+1)][(i-ii+5)]*scv_12;
                    scv_8=scv_8+A_copy[(k-kk+1)][(i-ii+6)]*scv_12;
                    scv_13=scv_13+A_copy[(k-kk+1)][(i-ii+7)]*scv_12;
                    scv_16=scv_16+A_copy[(k-kk+2)][(i-ii)]*scv_14;
                    scv_11=scv_11+A_copy[(k-kk+2)][(i-ii+1)]*scv_14;
                    scv_5=scv_5+A_copy[(k-kk+2)][(i-ii+2)]*scv_14;
                    scv_7=scv_7+A_copy[(k-kk+2)][(i-ii+3)]*scv_14;
                    scv_4=scv_4+A_copy[(k-kk+2)][(i-ii+4)]*scv_14;
                    scv_10=scv_10+A_copy[(k-kk+2)][(i-ii+5)]*scv_14;
                    scv_8=scv_8+A_copy[(k-kk+2)][(i-ii+6)]*scv_14;
                    scv_13=scv_13+A_copy[(k-kk+2)][(i-ii+7)]*scv_14;
                    scv_16=scv_16+A_copy[(k-kk+3)][(i-ii)]*scv_9;
                    scv_11=scv_11+A_copy[(k-kk+3)][(i-ii+1)]*scv_9;
                    scv_5=scv_5+A_copy[(k-kk+3)][(i-ii+2)]*scv_9;
                    scv_7=scv_7+A_copy[(k-kk+3)][(i-ii+3)]*scv_9;
                    scv_4=scv_4+A_copy[(k-kk+3)][(i-ii+4)]*scv_9;
                    scv_10=scv_10+A_copy[(k-kk+3)][(i-ii+5)]*scv_9;
                    scv_8=scv_8+A_copy[(k-kk+3)][(i-ii+6)]*scv_9;
                    scv_13=scv_13+A_copy[(k-kk+3)][(i-ii+7)]*scv_9;
                    scv_16=scv_16+A_copy[(k-kk+4)][(i-ii)]*scv_6;
                    scv_11=scv_11+A_copy[(k-kk+4)][(i-ii+1)]*scv_6;
                    scv_5=scv_5+A_copy[(k-kk+4)][(i-ii+2)]*scv_6;
                    scv_7=scv_7+A_copy[(k-kk+4)][(i-ii+3)]*scv_6;
                    scv_4=scv_4+A_copy[(k-kk+4)][(i-ii+4)]*scv_6;
                    scv_10=scv_10+A_copy[(k-kk+4)][(i-ii+5)]*scv_6;
                    scv_8=scv_8+A_copy[(k-kk+4)][(i-ii+6)]*scv_6;
                    scv_13=scv_13+A_copy[(k-kk+4)][(i-ii+7)]*scv_6;
                    scv_16=scv_16+A_copy[(k-kk+5)][(i-ii)]*scv_3;
                    scv_11=scv_11+A_copy[(k-kk+5)][(i-ii+1)]*scv_3;
                    scv_5=scv_5+A_copy[(k-kk+5)][(i-ii+2)]*scv_3;
                    scv_7=scv_7+A_copy[(k-kk+5)][(i-ii+3)]*scv_3;
                    scv_4=scv_4+A_copy[(k-kk+5)][(i-ii+4)]*scv_3;
                    scv_10=scv_10+A_copy[(k-kk+5)][(i-ii+5)]*scv_3;
                    scv_8=scv_8+A_copy[(k-kk+5)][(i-ii+6)]*scv_3;
                    scv_13=scv_13+A_copy[(k-kk+5)][(i-ii+7)]*scv_3;
                    scv_16=scv_16+A_copy[(k-kk+6)][(i-ii)]*scv_2;
                    scv_11=scv_11+A_copy[(k-kk+6)][(i-ii+1)]*scv_2;
                    scv_5=scv_5+A_copy[(k-kk+6)][(i-ii+2)]*scv_2;
                    scv_7=scv_7+A_copy[(k-kk+6)][(i-ii+3)]*scv_2;
                    scv_4=scv_4+A_copy[(k-kk+6)][(i-ii+4)]*scv_2;
                    scv_10=scv_10+A_copy[(k-kk+6)][(i-ii+5)]*scv_2;
                    scv_8=scv_8+A_copy[(k-kk+6)][(i-ii+6)]*scv_2;
                    scv_13=scv_13+A_copy[(k-kk+6)][(i-ii+7)]*scv_2;
                    scv_16=scv_16+A_copy[(k-kk+7)][(i-ii)]*scv_15;
                    scv_11=scv_11+A_copy[(k-kk+7)][(i-ii+1)]*scv_15;
                    scv_5=scv_5+A_copy[(k-kk+7)][(i-ii+2)]*scv_15;
                    scv_7=scv_7+A_copy[(k-kk+7)][(i-ii+3)]*scv_15;
                    scv_4=scv_4+A_copy[(k-kk+7)][(i-ii+4)]*scv_15;
                    scv_10=scv_10+A_copy[(k-kk+7)][(i-ii+5)]*scv_15;
                    scv_8=scv_8+A_copy[(k-kk+7)][(i-ii+6)]*scv_15;
                    scv_13=scv_13+A_copy[(k-kk+7)][(i-ii+7)]*scv_15;
                    C_copy[(i-ii+4)][(j-jj)]=scv_4;
                    C_copy[(i-ii+2)][(j-jj)]=scv_5;
                    C_copy[(i-ii+3)][(j-jj)]=scv_7;
                    C_copy[(i-ii+6)][(j-jj)]=scv_8;
                    C_copy[(i-ii+5)][(j-jj)]=scv_10;
                    C_copy[(i-ii+1)][(j-jj)]=scv_11;
                    C_copy[(i-ii+7)][(j-jj)]=scv_13;
                    C_copy[(i-ii)][(j-jj)]=scv_16;
                  }
                }
                for (; i<=min(M-1,ii+127); i=i+1) {
                  register int cbv_3;
                  cbv_3=min(N-1,jj+63);
#pragma ivdep
#pragma vector always
                  for (j=jj; j<=cbv_3; j=j+1) {
                    double scv_17;
                    scv_17=C_copy[(i-ii)][(j-jj)];
                    scv_17=scv_17+A_copy[(k-kk)][(i-ii)]*B[k][j];
                    scv_17=scv_17+A_copy[(k-kk+1)][(i-ii)]*B[(k+1)][j];
                    scv_17=scv_17+A_copy[(k-kk+2)][(i-ii)]*B[(k+2)][j];
                    scv_17=scv_17+A_copy[(k-kk+3)][(i-ii)]*B[(k+3)][j];
                    scv_17=scv_17+A_copy[(k-kk+4)][(i-ii)]*B[(k+4)][j];
                    scv_17=scv_17+A_copy[(k-kk+5)][(i-ii)]*B[(k+5)][j];
                    scv_17=scv_17+A_copy[(k-kk+6)][(i-ii)]*B[(k+6)][j];
                    scv_17=scv_17+A_copy[(k-kk+7)][(i-ii)]*B[(k+7)][j];
                    C_copy[(i-ii)][(j-jj)]=scv_17;
                  }
                }
              }
              for (; k<=min(K-1,kk+255); k=k+1) {
                for (i=ii; i<=min(M-1,ii+127)-7; i=i+8) {
                  register int cbv_4;
                  cbv_4=min(N-1,jj+63);
#pragma ivdep
#pragma vector always
                  for (j=jj; j<=cbv_4; j=j+1) {
                    double scv_18, scv_19, scv_20, scv_21, scv_22, scv_23, scv_24, scv_25;
                    double scv_26;
                    scv_18=C_copy[(i-ii+7)][(j-jj)];
                    scv_19=B[k][j];
                    scv_20=C_copy[(i-ii+5)][(j-jj)];
                    scv_21=C_copy[(i-ii+3)][(j-jj)];
                    scv_22=C_copy[(i-ii+4)][(j-jj)];
                    scv_23=C_copy[(i-ii+6)][(j-jj)];
                    scv_24=C_copy[(i-ii+1)][(j-jj)];
                    scv_25=C_copy[(i-ii)][(j-jj)];
                    scv_26=C_copy[(i-ii+2)][(j-jj)];
                    scv_25=scv_25+A_copy[(k-kk)][(i-ii)]*scv_19;
                    scv_24=scv_24+A_copy[(k-kk)][(i-ii+1)]*scv_19;
                    scv_26=scv_26+A_copy[(k-kk)][(i-ii+2)]*scv_19;
                    scv_21=scv_21+A_copy[(k-kk)][(i-ii+3)]*scv_19;
                    scv_22=scv_22+A_copy[(k-kk)][(i-ii+4)]*scv_19;
                    scv_20=scv_20+A_copy[(k-kk)][(i-ii+5)]*scv_19;
                    scv_23=scv_23+A_copy[(k-kk)][(i-ii+6)]*scv_19;
                    scv_18=scv_18+A_copy[(k-kk)][(i-ii+7)]*scv_19;
                    C_copy[(i-ii+7)][(j-jj)]=scv_18;
                    C_copy[(i-ii+5)][(j-jj)]=scv_20;
                    C_copy[(i-ii+3)][(j-jj)]=scv_21;
                    C_copy[(i-ii+4)][(j-jj)]=scv_22;
                    C_copy[(i-ii+6)][(j-jj)]=scv_23;
                    C_copy[(i-ii+1)][(j-jj)]=scv_24;
                    C_copy[(i-ii)][(j-jj)]=scv_25;
                    C_copy[(i-ii+2)][(j-jj)]=scv_26;
                  }
                }
                for (; i<=min(M-1,ii+127); i=i+1) {
                  register int cbv_5;
                  cbv_5=min(N-1,jj+63);
#pragma ivdep
#pragma vector always
                  for (j=jj; j<=cbv_5; j=j+1) {
                    double scv_27;
                    scv_27=C_copy[(i-ii)][(j-jj)];
                    scv_27=scv_27+A_copy[(k-kk)][(i-ii)]*B[k][j];
                    C_copy[(i-ii)][(j-jj)]=scv_27;
                  }
                }
              }
              for (i=ii; i<=min(M-1,ii+127); i=i+1) 
                for (j=jj; j<=min(N-1,jj+63); j=j+1) 
                  C[i][j]=C_copy[(i-ii)][(j-jj)];
            }
          }
        }
      }
    }
  }
}

    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return 1;
}
                                                                   
