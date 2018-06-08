#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

double A[N][N+20];
double B[N][N+20];

void init_arrays() {
    int i,j;
    for (i=0; i<N; i++)
        for (j=0; j<N; j++) {
            B[i][j] = (i+j) % 5 + 1;
            if (i < j)
                A[i][j] = (i+j) % 5 + 1;
            else if (i == j)
                A[i][j] = 1;
            else
                A[i][j] = -1;
        }
}

void print_array() {
    int i, j;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fprintf(stderr, "%lf ", round(B[i][j]));
            if (j%80 == 79) fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
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




        register int i,j,k,t;
        register int c1t, c2t, c3t, c4t, c5t, c6t, c7t, c8t, c9t, c10t, c11t, c12t;
        register int newlb_c1, newlb_c2, newlb_c3, newlb_c4, newlb_c5, newlb_c6,
                 newlb_c7, newlb_c8, newlb_c9, newlb_c10, newlb_c11, newlb_c12;
        register int newub_c1, newub_c2, newub_c3, newub_c4, newub_c5, newub_c6,
                 newub_c7, newub_c8, newub_c9, newub_c10, newub_c11, newub_c12;



#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))


        int c1, c2, c3, c4, c5, c6, c7, c8, c9;

        register int lb, ub, lb1, ub1, lb2, ub2;
        /* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 0.26s. */
        lb1=0;
        ub1=floord(N-1,1024);
        #pragma omp parallel for shared(lb1,ub1) private(c1,c2,c3,c4,c5,c6,c7,c8,c9)
        for (c1=lb1; c1<=ub1; c1++) {
            {
                for (c2=0; c2<=floord(N-2,512); c2++ ) {
                    for (c3=max(0,ceild(256*c2-127,128)); c3<=floord(N-1,256); c3++ ) {
                        for (c4=max(8*c1,0); c4<=min(floord(N-1,128),8*c1+7); c4++ ) {
                            for (c5=max(0,4*c2); c5<=min(min(4*c2+3,floord(N-2,128)),floord(128*c3+127,64)); c5++ ) {
                                for (c8t=max(0,128*c5); c8t<=min(min(128*c5+127,N-2),256*c3+254)-10; c8t=c8t+11) {
                                    newlb_c9=-2147483648;
                                    newub_c9=min(N-1,256*c3+255);
                                    register int cbv_1;
                                    cbv_1=c8t+10;
#pragma ivdep
#pragma vector always
                                    for (c8=c8t; c8<=cbv_1; c8=c8+1) {
                                        newlb_c9=max(newlb_c9,max(c8+1,256*c3));
                                    }
                                    for (c8=c8t; c8<=c8t+10; c8=c8+1) {
                                        for (c9=max(c8+1,256*c3); c9<=newlb_c9-1; c9=c9+1) {
                                            register int cbv_2, cbv_3;
                                            cbv_2=max(0,128*c4);
                                            cbv_3=min(N-1,128*c4+127);
#pragma ivdep
#pragma vector always
                                            for (c7=cbv_2; c7<=cbv_3; c7++ ) {
                                                double scv_1;
                                                scv_1=B[c8][c7];
                                                scv_1=alpha*A[c8][c9]*B[c9][c7]+scv_1;
                                                B[c8][c7]=scv_1;
                                            }
                                        }
                                    }
                                    for (c9t=newlb_c9; c9t<=newub_c9-7; c9t=c9t+8) {
                                        register int cbv_4, cbv_5;
                                        cbv_4=max(0,128*c4);
                                        cbv_5=min(N-1,128*c4+127);
#pragma ivdep
#pragma vector always
                                        for (c7=cbv_4; c7<=cbv_5; c7++ ) {
                                            double scv_2, scv_3, scv_4, scv_5, scv_6, scv_7, scv_8, scv_9;
                                            double scv_10, scv_11, scv_12, scv_13, scv_14, scv_15, scv_16, scv_17;
                                            double scv_18, scv_19, scv_20;
                                            scv_2=B[(c9t+1)][c7];
                                            scv_3=B[(c9t+2)][c7];
                                            scv_4=B[(c9t+3)][c7];
                                            scv_5=B[(c9t+6)][c7];
                                            scv_6=B[(c8t+8)][c7];
                                            scv_7=B[(c8t+2)][c7];
                                            scv_8=B[(c8t+4)][c7];
                                            scv_9=B[(c8t+1)][c7];
                                            scv_10=B[(c8t+10)][c7];
                                            scv_11=B[(c9t+5)][c7];
                                            scv_12=B[(c8t+5)][c7];
                                            scv_13=B[c9t][c7];
                                            scv_14=B[(c9t+7)][c7];
                                            scv_15=B[(c8t+9)][c7];
                                            scv_16=B[(c8t+7)][c7];
                                            scv_17=B[(c9t+4)][c7];
                                            scv_18=B[c8t][c7];
                                            scv_19=B[(c8t+3)][c7];
                                            scv_20=B[(c8t+6)][c7];
                                            scv_18=alpha*A[c8t][c9t]*scv_13+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+1)]*scv_2+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+2)]*scv_3+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+3)]*scv_4+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+4)]*scv_17+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+5)]*scv_11+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+6)]*scv_5+scv_18;
                                            scv_18=alpha*A[c8t][(c9t+7)]*scv_14+scv_18;
                                            scv_9=alpha*A[(c8t+1)][c9t]*scv_13+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+1)]*scv_2+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+2)]*scv_3+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+3)]*scv_4+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+4)]*scv_17+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+5)]*scv_11+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+6)]*scv_5+scv_9;
                                            scv_9=alpha*A[(c8t+1)][(c9t+7)]*scv_14+scv_9;
                                            scv_7=alpha*A[(c8t+2)][c9t]*scv_13+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+1)]*scv_2+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+2)]*scv_3+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+3)]*scv_4+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+4)]*scv_17+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+5)]*scv_11+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+6)]*scv_5+scv_7;
                                            scv_7=alpha*A[(c8t+2)][(c9t+7)]*scv_14+scv_7;
                                            scv_19=alpha*A[(c8t+3)][c9t]*scv_13+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+1)]*scv_2+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+2)]*scv_3+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+3)]*scv_4+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+4)]*scv_17+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+5)]*scv_11+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+6)]*scv_5+scv_19;
                                            scv_19=alpha*A[(c8t+3)][(c9t+7)]*scv_14+scv_19;
                                            scv_8=alpha*A[(c8t+4)][c9t]*scv_13+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+1)]*scv_2+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+2)]*scv_3+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+3)]*scv_4+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+4)]*scv_17+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+5)]*scv_11+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+6)]*scv_5+scv_8;
                                            scv_8=alpha*A[(c8t+4)][(c9t+7)]*scv_14+scv_8;
                                            scv_12=alpha*A[(c8t+5)][c9t]*scv_13+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+1)]*scv_2+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+2)]*scv_3+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+3)]*scv_4+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+4)]*scv_17+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+5)]*scv_11+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+6)]*scv_5+scv_12;
                                            scv_12=alpha*A[(c8t+5)][(c9t+7)]*scv_14+scv_12;
                                            scv_20=alpha*A[(c8t+6)][c9t]*scv_13+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+1)]*scv_2+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+2)]*scv_3+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+3)]*scv_4+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+4)]*scv_17+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+5)]*scv_11+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+6)]*scv_5+scv_20;
                                            scv_20=alpha*A[(c8t+6)][(c9t+7)]*scv_14+scv_20;
                                            scv_16=alpha*A[(c8t+7)][c9t]*scv_13+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+1)]*scv_2+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+2)]*scv_3+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+3)]*scv_4+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+4)]*scv_17+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+5)]*scv_11+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+6)]*scv_5+scv_16;
                                            scv_16=alpha*A[(c8t+7)][(c9t+7)]*scv_14+scv_16;
                                            scv_6=alpha*A[(c8t+8)][c9t]*scv_13+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+1)]*scv_2+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+2)]*scv_3+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+3)]*scv_4+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+4)]*scv_17+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+5)]*scv_11+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+6)]*scv_5+scv_6;
                                            scv_6=alpha*A[(c8t+8)][(c9t+7)]*scv_14+scv_6;
                                            scv_15=alpha*A[(c8t+9)][c9t]*scv_13+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+1)]*scv_2+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+2)]*scv_3+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+3)]*scv_4+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+4)]*scv_17+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+5)]*scv_11+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+6)]*scv_5+scv_15;
                                            scv_15=alpha*A[(c8t+9)][(c9t+7)]*scv_14+scv_15;
                                            scv_10=alpha*A[(c8t+10)][c9t]*scv_13+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+1)]*scv_2+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+2)]*scv_3+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+3)]*scv_4+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+4)]*scv_17+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+5)]*scv_11+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+6)]*scv_5+scv_10;
                                            scv_10=alpha*A[(c8t+10)][(c9t+7)]*scv_14+scv_10;
                                            B[(c8t+8)][c7]=scv_6;
                                            B[(c8t+2)][c7]=scv_7;
                                            B[(c8t+4)][c7]=scv_8;
                                            B[(c8t+1)][c7]=scv_9;
                                            B[(c8t+10)][c7]=scv_10;
                                            B[(c8t+5)][c7]=scv_12;
                                            B[(c8t+9)][c7]=scv_15;
                                            B[(c8t+7)][c7]=scv_16;
                                            B[c8t][c7]=scv_18;
                                            B[(c8t+3)][c7]=scv_19;
                                            B[(c8t+6)][c7]=scv_20;
                                        }
                                    }
                                    for (c9=c9t; c9<=newub_c9; c9=c9+1) {
                                        register int cbv_6, cbv_7;
                                        cbv_6=max(0,128*c4);
                                        cbv_7=min(N-1,128*c4+127);
#pragma ivdep
#pragma vector always
                                        for (c7=cbv_6; c7<=cbv_7; c7++ ) {
                                            double scv_21, scv_22, scv_23, scv_24, scv_25, scv_26, scv_27, scv_28;
                                            double scv_29, scv_30, scv_31, scv_32;
                                            scv_21=B[(c8t+7)][c7];
                                            scv_22=B[c8t][c7];
                                            scv_23=B[(c8t+10)][c7];
                                            scv_24=B[(c8t+5)][c7];
                                            scv_25=B[(c8t+8)][c7];
                                            scv_26=B[(c8t+3)][c7];
                                            scv_27=B[(c8t+6)][c7];
                                            scv_28=B[(c8t+2)][c7];
                                            scv_29=B[(c8t+1)][c7];
                                            scv_30=B[c9][c7];
                                            scv_31=B[(c8t+9)][c7];
                                            scv_32=B[(c8t+4)][c7];
                                            scv_22=alpha*A[c8t][c9]*scv_30+scv_22;
                                            scv_29=alpha*A[(c8t+1)][c9]*scv_30+scv_29;
                                            scv_28=alpha*A[(c8t+2)][c9]*scv_30+scv_28;
                                            scv_26=alpha*A[(c8t+3)][c9]*scv_30+scv_26;
                                            scv_32=alpha*A[(c8t+4)][c9]*scv_30+scv_32;
                                            scv_24=alpha*A[(c8t+5)][c9]*scv_30+scv_24;
                                            scv_27=alpha*A[(c8t+6)][c9]*scv_30+scv_27;
                                            scv_21=alpha*A[(c8t+7)][c9]*scv_30+scv_21;
                                            scv_25=alpha*A[(c8t+8)][c9]*scv_30+scv_25;
                                            scv_31=alpha*A[(c8t+9)][c9]*scv_30+scv_31;
                                            scv_23=alpha*A[(c8t+10)][c9]*scv_30+scv_23;
                                            B[(c8t+7)][c7]=scv_21;
                                            B[c8t][c7]=scv_22;
                                            B[(c8t+10)][c7]=scv_23;
                                            B[(c8t+5)][c7]=scv_24;
                                            B[(c8t+8)][c7]=scv_25;
                                            B[(c8t+3)][c7]=scv_26;
                                            B[(c8t+6)][c7]=scv_27;
                                            B[(c8t+2)][c7]=scv_28;
                                            B[(c8t+1)][c7]=scv_29;
                                            B[(c8t+9)][c7]=scv_31;
                                            B[(c8t+4)][c7]=scv_32;
                                        }
                                    }
                                    for (c8=c8t; c8<=c8t+10; c8=c8+1) {
                                        for (c9=newub_c9+1; c9<=min(N-1,256*c3+255); c9=c9+1) {
                                            register int cbv_8, cbv_9;
                                            cbv_8=max(0,128*c4);
                                            cbv_9=min(N-1,128*c4+127);
#pragma ivdep
#pragma vector always
                                            for (c7=cbv_8; c7<=cbv_9; c7++ ) {
                                                double scv_33;
                                                scv_33=B[c8][c7];
                                                scv_33=alpha*A[c8][c9]*B[c9][c7]+scv_33;
                                                B[c8][c7]=scv_33;
                                            }
                                        }
                                    }
                                }
                                for (c8=c8t; c8<=min(min(128*c5+127,N-2),256*c3+254); c8=c8+1) {
                                    for (c9t=max(c8+1,256*c3); c9t<=min(N-1,256*c3+255)-7; c9t=c9t+8) {
                                        register int cbv_10, cbv_11;
                                        cbv_10=max(0,128*c4);
                                        cbv_11=min(N-1,128*c4+127);
#pragma ivdep
#pragma vector always
                                        for (c7=cbv_10; c7<=cbv_11; c7++ ) {
                                            double scv_34;
                                            scv_34=B[c8][c7];
                                            scv_34=alpha*A[c8][c9t]*B[c9t][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+1)]*B[(c9t+1)][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+2)]*B[(c9t+2)][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+3)]*B[(c9t+3)][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+4)]*B[(c9t+4)][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+5)]*B[(c9t+5)][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+6)]*B[(c9t+6)][c7]+scv_34;
                                            scv_34=alpha*A[c8][(c9t+7)]*B[(c9t+7)][c7]+scv_34;
                                            B[c8][c7]=scv_34;
                                        }
                                    }
                                    for (c9=c9t; c9<=min(N-1,256*c3+255); c9=c9+1) {
                                        register int cbv_12, cbv_13;
                                        cbv_12=max(0,128*c4);
                                        cbv_13=min(N-1,128*c4+127);
#pragma ivdep
#pragma vector always
                                        for (c7=cbv_12; c7<=cbv_13; c7++ ) {
                                            double scv_35;
                                            scv_35=B[c8][c7];
                                            scv_35=alpha*A[c8][c9]*B[c9][c7]+scv_35;
                                            B[c8][c7]=scv_35;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
        /* End of CLooG code */


        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    //print_array();

    return ((int) B[0][0]);
}
