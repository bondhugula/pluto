
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

double A[N][N+17];

void init_arrays() {
    int i, j;
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            A[i][j] = i*i+j*j;
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
        /* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 0.18s. */
        for (c1=-1; c1<=floord(2*T+N-4,32); c1++) {
            lb1=max(max(0,ceild(32*c1-T+1,32)),ceild(16*c1-15,32));
            ub1=min(min(floord(32*c1+31,32),floord(32*c1+N+29,64)),floord(T+N-3,32));
            #pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5,c6,c7,c8,c9)
            for (c2=lb1; c2<=ub1; c2++) {
                {
                    for (c3=max(max(max(max(ceild(64*c2-N-508,512),ceild(16*c1-255,256)),ceild(16*c2-255,256)),0),ceild(64*c1-64*c2-509,512)); c3<=min(min(min(min(floord(T+N-3,256),floord(32*c2+T+N+28,512)),floord(32*c1+N+60,512)),floord(32*c1-32*c2+N+29,256)),floord(64*c2+N+59,512)); c3++ ) {
                        for (c6=max(max(max(max(max(ceild(16*c1-31,32),0),ceild(16*c2-31,32)),8*c3),ceild(64*c2-N-60,64)),ceild(64*c1-64*c2-61,64)); c6<=min(min(min(min(min(8*c3+7,floord(32*c1-32*c2+N+29,32)),floord(32*c2+T+N+28,64)),floord(64*c2+N+59,64)),floord(32*c1+N+60,64)),floord(T+N-3,32)); c6++ ) {
                            for (c7t=max(max(max(max(32*c6-N+2,-32*c2+64*c6-N-29),32*c2-N+2),0),32*c1-32*c2); c7t<=min(min(min(min(32*c1-32*c2+31,T-1),-32*c2+64*c6+62),32*c2+30),floord(64*c6+61,2))-7; c7t=c7t+8) {
                                for (c7=c7t; c7<=c7t+7; c7=c7+1) {
                                    for (c8=max(max(64*c6-c7-N+2,c7+1),32*c2); c8<=min(min(32*c2+31,c7+N-2),64*c6-c7+62); c8++ ) {
                                        register int cbv_1, cbv_2;
                                        cbv_1=max(64*c6,c7+c8+1);
                                        cbv_2=min(c7+c8+N-2,64*c6+63)-7;
#pragma ivdep
#pragma vector always
                                        for (c9t=cbv_1; c9t<=cbv_2; c9t=c9t+8) {
                                            double scv_1, scv_2, scv_3, scv_4, scv_5, scv_6, scv_7, scv_8;
                                            scv_1=A[-c7+c8][-c7-c8+(c9t+7)];
                                            scv_2=A[-c7+c8][-c7-c8+(c9t+1)];
                                            scv_3=A[-c7+c8][-c7-c8+(c9t+3)];
                                            scv_4=A[-c7+c8][-c7-c8+(c9t+2)];
                                            scv_5=A[-c7+c8][-c7-c8+(c9t+4)];
                                            scv_6=A[-c7+c8][-c7-c8+c9t];
                                            scv_7=A[-c7+c8][-c7-c8+(c9t+5)];
                                            scv_8=A[-c7+c8][-c7-c8+(c9t+6)];
                                            scv_6=(A[1+-c7+c8][1+-c7-c8+c9t]+A[1+-c7+c8][-c7-c8+c9t]+A[1+-c7+c8][-c7-c8+c9t-1]+A[-c7+c8][1+-c7-c8+c9t]+scv_6+A[-c7+c8][-c7-c8+c9t-1]+A[-c7+c8-1][1+-c7-c8+c9t]+A[-c7+c8-1][-c7-c8+c9t]+A[-c7+c8-1][-c7-c8+c9t-1])/9;
                                            scv_2=(A[1+-c7+c8][1+-c7-c8+(c9t+1)]+A[1+-c7+c8][-c7-c8+(c9t+1)]+A[1+-c7+c8][-c7-c8+(c9t+1)-1]+A[-c7+c8][1+-c7-c8+(c9t+1)]+scv_2+A[-c7+c8][-c7-c8+(c9t+1)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+1)]+A[-c7+c8-1][-c7-c8+(c9t+1)]+A[-c7+c8-1][-c7-c8+(c9t+1)-1])/9;
                                            scv_4=(A[1+-c7+c8][1+-c7-c8+(c9t+2)]+A[1+-c7+c8][-c7-c8+(c9t+2)]+A[1+-c7+c8][-c7-c8+(c9t+2)-1]+A[-c7+c8][1+-c7-c8+(c9t+2)]+scv_4+A[-c7+c8][-c7-c8+(c9t+2)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+2)]+A[-c7+c8-1][-c7-c8+(c9t+2)]+A[-c7+c8-1][-c7-c8+(c9t+2)-1])/9;
                                            scv_3=(A[1+-c7+c8][1+-c7-c8+(c9t+3)]+A[1+-c7+c8][-c7-c8+(c9t+3)]+A[1+-c7+c8][-c7-c8+(c9t+3)-1]+A[-c7+c8][1+-c7-c8+(c9t+3)]+scv_3+A[-c7+c8][-c7-c8+(c9t+3)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+3)]+A[-c7+c8-1][-c7-c8+(c9t+3)]+A[-c7+c8-1][-c7-c8+(c9t+3)-1])/9;
                                            scv_5=(A[1+-c7+c8][1+-c7-c8+(c9t+4)]+A[1+-c7+c8][-c7-c8+(c9t+4)]+A[1+-c7+c8][-c7-c8+(c9t+4)-1]+A[-c7+c8][1+-c7-c8+(c9t+4)]+scv_5+A[-c7+c8][-c7-c8+(c9t+4)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+4)]+A[-c7+c8-1][-c7-c8+(c9t+4)]+A[-c7+c8-1][-c7-c8+(c9t+4)-1])/9;
                                            scv_7=(A[1+-c7+c8][1+-c7-c8+(c9t+5)]+A[1+-c7+c8][-c7-c8+(c9t+5)]+A[1+-c7+c8][-c7-c8+(c9t+5)-1]+A[-c7+c8][1+-c7-c8+(c9t+5)]+scv_7+A[-c7+c8][-c7-c8+(c9t+5)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+5)]+A[-c7+c8-1][-c7-c8+(c9t+5)]+A[-c7+c8-1][-c7-c8+(c9t+5)-1])/9;
                                            scv_8=(A[1+-c7+c8][1+-c7-c8+(c9t+6)]+A[1+-c7+c8][-c7-c8+(c9t+6)]+A[1+-c7+c8][-c7-c8+(c9t+6)-1]+A[-c7+c8][1+-c7-c8+(c9t+6)]+scv_8+A[-c7+c8][-c7-c8+(c9t+6)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+6)]+A[-c7+c8-1][-c7-c8+(c9t+6)]+A[-c7+c8-1][-c7-c8+(c9t+6)-1])/9;
                                            scv_1=(A[1+-c7+c8][1+-c7-c8+(c9t+7)]+A[1+-c7+c8][-c7-c8+(c9t+7)]+A[1+-c7+c8][-c7-c8+(c9t+7)-1]+A[-c7+c8][1+-c7-c8+(c9t+7)]+scv_1+A[-c7+c8][-c7-c8+(c9t+7)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+7)]+A[-c7+c8-1][-c7-c8+(c9t+7)]+A[-c7+c8-1][-c7-c8+(c9t+7)-1])/9;
                                            A[-c7+c8][-c7-c8+(c9t+7)]=scv_1;
                                            A[-c7+c8][-c7-c8+(c9t+1)]=scv_2;
                                            A[-c7+c8][-c7-c8+(c9t+3)]=scv_3;
                                            A[-c7+c8][-c7-c8+(c9t+2)]=scv_4;
                                            A[-c7+c8][-c7-c8+(c9t+4)]=scv_5;
                                            A[-c7+c8][-c7-c8+c9t]=scv_6;
                                            A[-c7+c8][-c7-c8+(c9t+5)]=scv_7;
                                            A[-c7+c8][-c7-c8+(c9t+6)]=scv_8;
                                        }
                                        register int cbv_3;
                                        cbv_3=min(c7+c8+N-2,64*c6+63);
#pragma ivdep
#pragma vector always
                                        for (c9=c9t; c9<=cbv_3; c9=c9+1) {
                                            double scv_9;
                                            scv_9=A[-c7+c8][-c7-c8+c9];
                                            scv_9=(A[1+-c7+c8][1+-c7-c8+c9]+A[1+-c7+c8][-c7-c8+c9]+A[1+-c7+c8][-c7-c8+c9-1]+A[-c7+c8][1+-c7-c8+c9]+scv_9+A[-c7+c8][-c7-c8+c9-1]+A[-c7+c8-1][1+-c7-c8+c9]+A[-c7+c8-1][-c7-c8+c9]+A[-c7+c8-1][-c7-c8+c9-1])/9;
                                            A[-c7+c8][-c7-c8+c9]=scv_9;
                                        }
                                    }
                                }
                            }
                            for (c7=c7t; c7<=min(min(min(min(32*c1-32*c2+31,T-1),-32*c2+64*c6+62),32*c2+30),floord(64*c6+61,2)); c7=c7+1) {
                                for (c8=max(max(64*c6-c7-N+2,c7+1),32*c2); c8<=min(min(32*c2+31,c7+N-2),64*c6-c7+62); c8++ ) {
                                    register int cbv_4, cbv_5;
                                    cbv_4=max(64*c6,c7+c8+1);
                                    cbv_5=min(c7+c8+N-2,64*c6+63)-7;
#pragma ivdep
#pragma vector always
                                    for (c9t=cbv_4; c9t<=cbv_5; c9t=c9t+8) {
                                        double scv_10, scv_11, scv_12, scv_13, scv_14, scv_15, scv_16, scv_17;
                                        scv_10=A[-c7+c8][-c7-c8+(c9t+7)];
                                        scv_11=A[-c7+c8][-c7-c8+(c9t+1)];
                                        scv_12=A[-c7+c8][-c7-c8+(c9t+3)];
                                        scv_13=A[-c7+c8][-c7-c8+(c9t+2)];
                                        scv_14=A[-c7+c8][-c7-c8+(c9t+4)];
                                        scv_15=A[-c7+c8][-c7-c8+c9t];
                                        scv_16=A[-c7+c8][-c7-c8+(c9t+5)];
                                        scv_17=A[-c7+c8][-c7-c8+(c9t+6)];
                                        scv_15=(A[1+-c7+c8][1+-c7-c8+c9t]+A[1+-c7+c8][-c7-c8+c9t]+A[1+-c7+c8][-c7-c8+c9t-1]+A[-c7+c8][1+-c7-c8+c9t]+scv_15+A[-c7+c8][-c7-c8+c9t-1]+A[-c7+c8-1][1+-c7-c8+c9t]+A[-c7+c8-1][-c7-c8+c9t]+A[-c7+c8-1][-c7-c8+c9t-1])/9;
                                        scv_11=(A[1+-c7+c8][1+-c7-c8+(c9t+1)]+A[1+-c7+c8][-c7-c8+(c9t+1)]+A[1+-c7+c8][-c7-c8+(c9t+1)-1]+A[-c7+c8][1+-c7-c8+(c9t+1)]+scv_11+A[-c7+c8][-c7-c8+(c9t+1)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+1)]+A[-c7+c8-1][-c7-c8+(c9t+1)]+A[-c7+c8-1][-c7-c8+(c9t+1)-1])/9;
                                        scv_13=(A[1+-c7+c8][1+-c7-c8+(c9t+2)]+A[1+-c7+c8][-c7-c8+(c9t+2)]+A[1+-c7+c8][-c7-c8+(c9t+2)-1]+A[-c7+c8][1+-c7-c8+(c9t+2)]+scv_13+A[-c7+c8][-c7-c8+(c9t+2)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+2)]+A[-c7+c8-1][-c7-c8+(c9t+2)]+A[-c7+c8-1][-c7-c8+(c9t+2)-1])/9;
                                        scv_12=(A[1+-c7+c8][1+-c7-c8+(c9t+3)]+A[1+-c7+c8][-c7-c8+(c9t+3)]+A[1+-c7+c8][-c7-c8+(c9t+3)-1]+A[-c7+c8][1+-c7-c8+(c9t+3)]+scv_12+A[-c7+c8][-c7-c8+(c9t+3)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+3)]+A[-c7+c8-1][-c7-c8+(c9t+3)]+A[-c7+c8-1][-c7-c8+(c9t+3)-1])/9;
                                        scv_14=(A[1+-c7+c8][1+-c7-c8+(c9t+4)]+A[1+-c7+c8][-c7-c8+(c9t+4)]+A[1+-c7+c8][-c7-c8+(c9t+4)-1]+A[-c7+c8][1+-c7-c8+(c9t+4)]+scv_14+A[-c7+c8][-c7-c8+(c9t+4)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+4)]+A[-c7+c8-1][-c7-c8+(c9t+4)]+A[-c7+c8-1][-c7-c8+(c9t+4)-1])/9;
                                        scv_16=(A[1+-c7+c8][1+-c7-c8+(c9t+5)]+A[1+-c7+c8][-c7-c8+(c9t+5)]+A[1+-c7+c8][-c7-c8+(c9t+5)-1]+A[-c7+c8][1+-c7-c8+(c9t+5)]+scv_16+A[-c7+c8][-c7-c8+(c9t+5)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+5)]+A[-c7+c8-1][-c7-c8+(c9t+5)]+A[-c7+c8-1][-c7-c8+(c9t+5)-1])/9;
                                        scv_17=(A[1+-c7+c8][1+-c7-c8+(c9t+6)]+A[1+-c7+c8][-c7-c8+(c9t+6)]+A[1+-c7+c8][-c7-c8+(c9t+6)-1]+A[-c7+c8][1+-c7-c8+(c9t+6)]+scv_17+A[-c7+c8][-c7-c8+(c9t+6)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+6)]+A[-c7+c8-1][-c7-c8+(c9t+6)]+A[-c7+c8-1][-c7-c8+(c9t+6)-1])/9;
                                        scv_10=(A[1+-c7+c8][1+-c7-c8+(c9t+7)]+A[1+-c7+c8][-c7-c8+(c9t+7)]+A[1+-c7+c8][-c7-c8+(c9t+7)-1]+A[-c7+c8][1+-c7-c8+(c9t+7)]+scv_10+A[-c7+c8][-c7-c8+(c9t+7)-1]+A[-c7+c8-1][1+-c7-c8+(c9t+7)]+A[-c7+c8-1][-c7-c8+(c9t+7)]+A[-c7+c8-1][-c7-c8+(c9t+7)-1])/9;
                                        A[-c7+c8][-c7-c8+(c9t+7)]=scv_10;
                                        A[-c7+c8][-c7-c8+(c9t+1)]=scv_11;
                                        A[-c7+c8][-c7-c8+(c9t+3)]=scv_12;
                                        A[-c7+c8][-c7-c8+(c9t+2)]=scv_13;
                                        A[-c7+c8][-c7-c8+(c9t+4)]=scv_14;
                                        A[-c7+c8][-c7-c8+c9t]=scv_15;
                                        A[-c7+c8][-c7-c8+(c9t+5)]=scv_16;
                                        A[-c7+c8][-c7-c8+(c9t+6)]=scv_17;
                                    }
                                    register int cbv_6;
                                    cbv_6=min(c7+c8+N-2,64*c6+63);
#pragma ivdep
#pragma vector always
                                    for (c9=c9t; c9<=cbv_6; c9=c9+1) {
                                        double scv_18;
                                        scv_18=A[-c7+c8][-c7-c8+c9];
                                        scv_18=(A[1+-c7+c8][1+-c7-c8+c9]+A[1+-c7+c8][-c7-c8+c9]+A[1+-c7+c8][-c7-c8+c9-1]+A[-c7+c8][1+-c7-c8+c9]+scv_18+A[-c7+c8][-c7-c8+c9-1]+A[-c7+c8-1][1+-c7-c8+c9]+A[-c7+c8-1][-c7-c8+c9]+A[-c7+c8-1][-c7-c8+c9-1])/9;
                                        A[-c7+c8][-c7-c8+c9]=scv_18;
                                    }
                                }
                            }
                        }
                    }
                }
                /*@ end @*/

            }
        }
        /* End of CLooG code */


        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    return ((int) A[0][0]);

}
