

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

void init_arrays() {
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

double rtclock() {
    struct timezone tzp;
    struct timeval tp;
    int stat;
    gettimeofday (&tp, &tzp);
    return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main() {
    init_arrays();

    double annot_t_start=0, annot_t_end=0, annot_t_total=0;
    int annot_i;

    for (annot_i=0; annot_i<REPS; annot_i++) {
        annot_t_start = rtclock();



        int i, j, k;
        int ii, jj, kk;
        int iii, jjj, kkk;

        {
            double C_copy;
            double B_copy;
            double A_copy;
            register int cbv_1;
            cbv_1=M-1;
            #pragma omp parallel for private(iii,jjj,kkk,ii,jj,kk,i,j,k,A_copy,B_copy,C_copy)
            for (iii=0; iii<=cbv_1; iii=iii+256) {
                for (kkk=0; kkk<=K-1; kkk=kkk+256) {
                    for (ii=iii; ii<=min(M-1,iii+224); ii=ii+32) {
                        for (i=ii; i<=min(M-1,ii+31)-3; i=i+4) {
                            for (k=kkk; k<=min(K-1,kkk+255)-3; k=k+4) {
                                register int cbv_2;
                                cbv_2=N-1;
#pragma ivdep
#pragma vector always
                                for (j=0; j<=cbv_2; j=j+1) {
                                    double scv_1, scv_2, scv_3, scv_4, scv_5, scv_6, scv_7, scv_8;
                                    scv_1=B[k][j];
                                    scv_2=B[(k+1)][j];
                                    scv_3=C[(i+2)][j];
                                    scv_4=B[(k+3)][j];
                                    scv_5=C[(i+3)][j];
                                    scv_6=B[(k+2)][j];
                                    scv_7=C[(i+1)][j];
                                    scv_8=C[i][j];
                                    scv_8=scv_8+A[i][k]*scv_1;
                                    scv_7=scv_7+A[(i+1)][k]*scv_1;
                                    scv_3=scv_3+A[(i+2)][k]*scv_1;
                                    scv_5=scv_5+A[(i+3)][k]*scv_1;
                                    scv_8=scv_8+A[i][(k+1)]*scv_2;
                                    scv_7=scv_7+A[(i+1)][(k+1)]*scv_2;
                                    scv_3=scv_3+A[(i+2)][(k+1)]*scv_2;
                                    scv_5=scv_5+A[(i+3)][(k+1)]*scv_2;
                                    scv_8=scv_8+A[i][(k+2)]*scv_6;
                                    scv_7=scv_7+A[(i+1)][(k+2)]*scv_6;
                                    scv_3=scv_3+A[(i+2)][(k+2)]*scv_6;
                                    scv_5=scv_5+A[(i+3)][(k+2)]*scv_6;
                                    scv_8=scv_8+A[i][(k+3)]*scv_4;
                                    scv_7=scv_7+A[(i+1)][(k+3)]*scv_4;
                                    scv_3=scv_3+A[(i+2)][(k+3)]*scv_4;
                                    scv_5=scv_5+A[(i+3)][(k+3)]*scv_4;
                                    C[(i+2)][j]=scv_3;
                                    C[(i+3)][j]=scv_5;
                                    C[(i+1)][j]=scv_7;
                                    C[i][j]=scv_8;
                                }
                            }
                            for (; k<=min(K-1,kkk+255); k=k+1) {
                                register int cbv_3;
                                cbv_3=N-1;
#pragma ivdep
#pragma vector always
                                for (j=0; j<=cbv_3; j=j+1) {
                                    double scv_9, scv_10, scv_11, scv_12, scv_13;
                                    scv_9=B[k][j];
                                    scv_10=C[(i+1)][j];
                                    scv_11=C[i][j];
                                    scv_12=C[(i+3)][j];
                                    scv_13=C[(i+2)][j];
                                    scv_11=scv_11+A[i][k]*scv_9;
                                    scv_10=scv_10+A[(i+1)][k]*scv_9;
                                    scv_13=scv_13+A[(i+2)][k]*scv_9;
                                    scv_12=scv_12+A[(i+3)][k]*scv_9;
                                    C[(i+1)][j]=scv_10;
                                    C[i][j]=scv_11;
                                    C[(i+3)][j]=scv_12;
                                    C[(i+2)][j]=scv_13;
                                }
                            }
                        }
                        for (; i<=min(M-1,ii+31); i=i+1) {
                            for (k=kkk; k<=min(K-1,kkk+255)-3; k=k+4) {
                                register int cbv_4;
                                cbv_4=N-1;
#pragma ivdep
#pragma vector always
                                for (j=0; j<=cbv_4; j=j+1) {
                                    double scv_14;
                                    scv_14=C[i][j];
                                    scv_14=scv_14+A[i][k]*B[k][j];
                                    scv_14=scv_14+A[i][(k+1)]*B[(k+1)][j];
                                    scv_14=scv_14+A[i][(k+2)]*B[(k+2)][j];
                                    scv_14=scv_14+A[i][(k+3)]*B[(k+3)][j];
                                    C[i][j]=scv_14;
                                }
                            }
                            for (; k<=min(K-1,kkk+255); k=k+1) {
                                register int cbv_5;
                                cbv_5=N-1;
#pragma ivdep
#pragma vector always
                                for (j=0; j<=cbv_5; j=j+1) {
                                    double scv_15;
                                    scv_15=C[i][j];
                                    scv_15=scv_15+A[i][k]*B[k][j];
                                    C[i][j]=scv_15;
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

