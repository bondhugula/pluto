

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

double a[T][N];

void init_input_vars() {
    int i, j;
    for (i=0; i<T; i++)
        for (j=0; j<N; j++)
            a[i][j] = i+((double)j)/N;
}

double rtclock() {
    struct timezone tzp;
    struct timeval tp;
    int stat;
    gettimeofday (&tp, &tzp);
    return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main() {
    init_input_vars();


    double orio_t_start=0, orio_t_end=0, orio_t_total=0;
    int orio_i;
    int t,i,j;

    for (orio_i=0; orio_i<REPS; orio_i++) {
        orio_t_start = rtclock();





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


        int c1, c2, c3, c4, c5, c6;

        register int lb, ub, lb1, ub1, lb2, ub2;
        /* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 0.01s. */
        for (c1=-1; c1<=floord(33*T+N-35,1024); c1++) {
            lb1=max(max(0,ceild(16*c1-511,528)),ceild(32*c1-T+1,32));
            ub1=min(min(floord(16*c1+15,16),floord(32*c1+N+29,1056)),floord(T+N-3,1024));
            #pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5,c6)
            for (c2=lb1; c2<=ub1; c2++) {
                /*@ begin Loop(
                transform Composite(
                permut = [['c5', 'c6']],
                  regtile = (['c5', 'c6'],[2, 8]),
                  scalarreplace = (False, 'double'),
                  vector = (True, ['ivdep','vector always']))

                    for (c4=max(max(0,ceild(16*c1-16*c2-63,64)),8*c2);c4<=min(min(8*c2+7,floord(32*c1-32*c2+N+29,128)),floord(T+N-3,128));c4++) {
                      for (c5=max(max(32*c1-32*c2,1),128*c4-N+2);c5<=min(min(128*c4+126,T-1),32*c1-32*c2+31);c5++) {
                        for (c6=max(c5+1,128*c4);c6<=min(128*c4+127,c5+N-2);c6++) {
                          a[c5][-c5+c6]=((double)(333))/1000*(a[c5-1][1+-c5+c6]+a[c5-1][-c5+c6]+a[c5-1][-c5+c6-1]) ;
                        }
                      }
                    }

                ) @*/for (c4=max(max(0,ceild(16*c1-16*c2-63,64)),8*c2); c4<=min(min(8*c2+7,floord(32*c1-32*c2+N+29,128)),floord(T+N-3,128)); c4++ ) {
                    for (c5t=max(max(32*c1-32*c2,1),128*c4-N+2); c5t<=min(min(128*c4+126,T-1),32*c1-32*c2+31)-1; c5t=c5t+2) {
                        newlb_c6=-2147483648;
                        newub_c6=2147483647;
                        register int cbv_1;
                        cbv_1=c5t+1;
#pragma ivdep
#pragma vector always
                        for (c5=c5t; c5<=cbv_1; c5=c5+1) {
                            newlb_c6=max(newlb_c6,max(c5+1,128*c4));
                            newub_c6=min(newub_c6,min(128*c4+127,c5+N-2));
                        }
                        for (c5=c5t; c5<=c5t+1; c5=c5+1) {
                            register int cbv_2, cbv_3;
                            cbv_2=max(c5+1,128*c4);
                            cbv_3=newlb_c6-1;
#pragma ivdep
#pragma vector always
                            for (c6=cbv_2; c6<=cbv_3; c6=c6+1) {
                                a[c5][-c5+c6]=((double)(333))/1000*(a[c5-1][1+-c5+c6]+a[c5-1][-c5+c6]+a[c5-1][-c5+c6-1]);
                            }
                        }
                        register int cbv_4;
                        cbv_4=newub_c6-7;
#pragma ivdep
#pragma vector always
                        for (c6t=newlb_c6; c6t<=cbv_4; c6t=c6t+8) {
                            a[c5t][-c5t+c6t]=((double)(333))/1000*(a[c5t-1][1+-c5t+c6t]+a[c5t-1][-c5t+c6t]+a[c5t-1][-c5t+c6t-1]);
                            a[c5t][-c5t+(c6t+1)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+1)]+a[c5t-1][-c5t+(c6t+1)]+a[c5t-1][-c5t+(c6t+1)-1]);
                            a[c5t][-c5t+(c6t+2)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+2)]+a[c5t-1][-c5t+(c6t+2)]+a[c5t-1][-c5t+(c6t+2)-1]);
                            a[c5t][-c5t+(c6t+3)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+3)]+a[c5t-1][-c5t+(c6t+3)]+a[c5t-1][-c5t+(c6t+3)-1]);
                            a[c5t][-c5t+(c6t+4)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+4)]+a[c5t-1][-c5t+(c6t+4)]+a[c5t-1][-c5t+(c6t+4)-1]);
                            a[c5t][-c5t+(c6t+5)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+5)]+a[c5t-1][-c5t+(c6t+5)]+a[c5t-1][-c5t+(c6t+5)-1]);
                            a[c5t][-c5t+(c6t+6)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+6)]+a[c5t-1][-c5t+(c6t+6)]+a[c5t-1][-c5t+(c6t+6)-1]);
                            a[c5t][-c5t+(c6t+7)]=((double)(333))/1000*(a[c5t-1][1+-c5t+(c6t+7)]+a[c5t-1][-c5t+(c6t+7)]+a[c5t-1][-c5t+(c6t+7)-1]);
                            a[(c5t+1)][-(c5t+1)+c6t]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+c6t]+a[(c5t+1)-1][-(c5t+1)+c6t]+a[(c5t+1)-1][-(c5t+1)+c6t-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+1)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+1)]+a[(c5t+1)-1][-(c5t+1)+(c6t+1)]+a[(c5t+1)-1][-(c5t+1)+(c6t+1)-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+2)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+2)]+a[(c5t+1)-1][-(c5t+1)+(c6t+2)]+a[(c5t+1)-1][-(c5t+1)+(c6t+2)-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+3)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+3)]+a[(c5t+1)-1][-(c5t+1)+(c6t+3)]+a[(c5t+1)-1][-(c5t+1)+(c6t+3)-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+4)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+4)]+a[(c5t+1)-1][-(c5t+1)+(c6t+4)]+a[(c5t+1)-1][-(c5t+1)+(c6t+4)-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+5)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+5)]+a[(c5t+1)-1][-(c5t+1)+(c6t+5)]+a[(c5t+1)-1][-(c5t+1)+(c6t+5)-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+6)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+6)]+a[(c5t+1)-1][-(c5t+1)+(c6t+6)]+a[(c5t+1)-1][-(c5t+1)+(c6t+6)-1]);
                            a[(c5t+1)][-(c5t+1)+(c6t+7)]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+(c6t+7)]+a[(c5t+1)-1][-(c5t+1)+(c6t+7)]+a[(c5t+1)-1][-(c5t+1)+(c6t+7)-1]);
                        }
#pragma ivdep
#pragma vector always
                        for (c6=c6t; c6<=newub_c6; c6=c6+1) {
                            a[c5t][-c5t+c6]=((double)(333))/1000*(a[c5t-1][1+-c5t+c6]+a[c5t-1][-c5t+c6]+a[c5t-1][-c5t+c6-1]);
                            a[(c5t+1)][-(c5t+1)+c6]=((double)(333))/1000*(a[(c5t+1)-1][1+-(c5t+1)+c6]+a[(c5t+1)-1][-(c5t+1)+c6]+a[(c5t+1)-1][-(c5t+1)+c6-1]);
                        }
                        for (c5=c5t; c5<=c5t+1; c5=c5+1) {
                            register int cbv_5, cbv_6;
                            cbv_5=newub_c6+1;
                            cbv_6=min(128*c4+127,c5+N-2);
#pragma ivdep
#pragma vector always
                            for (c6=cbv_5; c6<=cbv_6; c6=c6+1) {
                                a[c5][-c5+c6]=((double)(333))/1000*(a[c5-1][1+-c5+c6]+a[c5-1][-c5+c6]+a[c5-1][-c5+c6-1]);
                            }
                        }
                    }
                    for (c5=c5t; c5<=min(min(128*c4+126,T-1),32*c1-32*c2+31); c5=c5+1) {
                        register int cbv_7, cbv_8;
                        cbv_7=max(c5+1,128*c4);
                        cbv_8=min(128*c4+127,c5+N-2)-7;
#pragma ivdep
#pragma vector always
                        for (c6t=cbv_7; c6t<=cbv_8; c6t=c6t+8) {
                            a[c5][-c5+c6t]=((double)(333))/1000*(a[c5-1][1+-c5+c6t]+a[c5-1][-c5+c6t]+a[c5-1][-c5+c6t-1]);
                            a[c5][-c5+(c6t+1)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+1)]+a[c5-1][-c5+(c6t+1)]+a[c5-1][-c5+(c6t+1)-1]);
                            a[c5][-c5+(c6t+2)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+2)]+a[c5-1][-c5+(c6t+2)]+a[c5-1][-c5+(c6t+2)-1]);
                            a[c5][-c5+(c6t+3)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+3)]+a[c5-1][-c5+(c6t+3)]+a[c5-1][-c5+(c6t+3)-1]);
                            a[c5][-c5+(c6t+4)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+4)]+a[c5-1][-c5+(c6t+4)]+a[c5-1][-c5+(c6t+4)-1]);
                            a[c5][-c5+(c6t+5)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+5)]+a[c5-1][-c5+(c6t+5)]+a[c5-1][-c5+(c6t+5)-1]);
                            a[c5][-c5+(c6t+6)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+6)]+a[c5-1][-c5+(c6t+6)]+a[c5-1][-c5+(c6t+6)-1]);
                            a[c5][-c5+(c6t+7)]=((double)(333))/1000*(a[c5-1][1+-c5+(c6t+7)]+a[c5-1][-c5+(c6t+7)]+a[c5-1][-c5+(c6t+7)-1]);
                        }
                        register int cbv_9;
                        cbv_9=min(128*c4+127,c5+N-2);
#pragma ivdep
#pragma vector always
                        for (c6=c6t; c6<=cbv_9; c6=c6+1) {
                            a[c5][-c5+c6]=((double)(333))/1000*(a[c5-1][1+-c5+c6]+a[c5-1][-c5+c6]+a[c5-1][-c5+c6-1]);
                        }
                    }
                }
                /*@ end @*/

            }
        }
        /* End of CLooG code */


        orio_t_end = rtclock();
        orio_t_total += orio_t_end - orio_t_start;
    }

    orio_t_total = orio_t_total / REPS;
    printf("%f\n", orio_t_total);

    return a[0][0];
}
