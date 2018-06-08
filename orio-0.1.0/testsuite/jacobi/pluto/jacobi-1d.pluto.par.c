

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



#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S1(zT0,zT1,zT2,zT3,t,i)	{a[t][i]=((double)(333))/1000*(a[t-1][1+i]+a[t-1][i]+a[t-1][i-1]);}

        int c1, c2, c3, c4, c5, c6;

        register int lb, ub, lb1, ub1, lb2, ub2;
        register int lbv, ubv;

        /* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 0.02s. */
        for (c1=-1; c1<=floord(2*T+N-4,256); c1++) {
            lb1=max(max(0,ceild(128*c1-127,256)),ceild(256*c1-T+1,256));
            ub1=min(min(floord(128*c1+127,128),floord(256*c1+N+253,512)),floord(T+N-3,256));
            #pragma omp parallel for shared(c1,lb1,ub1) private(c2,c3,c4,c5,c6)
            for (c2=lb1; c2<=ub1; c2++) {
                for (c3=max(max(ceild(256*c2-N-29,32),8*c1-8*c2),0); c3<=min(min(8*c1-8*c2+7,floord(128*c2+127,16)),floord(T-1,32)); c3++) {
                    for (c4=max(max(0,8*c2),ceild(16*c3-15,16)); c4<=min(min(8*c2+7,floord(T+N-3,32)),floord(32*c3+N+29,32)); c4++) {
                        /*@ begin Loop(
                        transform UnrollJam(ufactor=32)
                                for (c5=max(max(32*c3,1),32*c4-N+2);c5<=min(min(T-1,32*c3+31),32*c4+30);c5++)
                        {
                        {
                        	lbv=max(32*c4,c5+1);
                        	ubv=min(32*c4+31,c5+N-2);
                        #pragma ivdep
                        #pragma vector always
                        	for (c6=lbv; c6<=ubv; c6++) {
                                    S1(c1-c2,-c1+2*c2,c3,-c3+c4,c5,-c5+c6) ;
                                  }
                        }
                        }
                        ) @*/{

                            for (c5 = max(max(32 * c3, 1), 32 * c4 - N + 2); c5 <= min(min(T - 1, 32 * c3 + 31), 32 * c4 + 30) - 31; c5 = c5 + 32)     {

                                {
                                    lbv=max(32*c4,c5+1);
                                    ubv=min(32*c4+31,c5+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, c5, -c5 + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+1)+1);
                                    ubv=min(32*c4+31,(c5+1)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 1), -(c5 + 1) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+2)+1);
                                    ubv=min(32*c4+31,(c5+2)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 2), -(c5 + 2) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+3)+1);
                                    ubv=min(32*c4+31,(c5+3)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 3), -(c5 + 3) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+4)+1);
                                    ubv=min(32*c4+31,(c5+4)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 4), -(c5 + 4) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+5)+1);
                                    ubv=min(32*c4+31,(c5+5)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 5), -(c5 + 5) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+6)+1);
                                    ubv=min(32*c4+31,(c5+6)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 6), -(c5 + 6) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+7)+1);
                                    ubv=min(32*c4+31,(c5+7)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 7), -(c5 + 7) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+8)+1);
                                    ubv=min(32*c4+31,(c5+8)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 8), -(c5 + 8) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+9)+1);
                                    ubv=min(32*c4+31,(c5+9)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 9), -(c5 + 9) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+10)+1);
                                    ubv=min(32*c4+31,(c5+10)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 10), -(c5 + 10) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+11)+1);
                                    ubv=min(32*c4+31,(c5+11)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 11), -(c5 + 11) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+12)+1);
                                    ubv=min(32*c4+31,(c5+12)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 12), -(c5 + 12) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+13)+1);
                                    ubv=min(32*c4+31,(c5+13)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 13), -(c5 + 13) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+14)+1);
                                    ubv=min(32*c4+31,(c5+14)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 14), -(c5 + 14) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+15)+1);
                                    ubv=min(32*c4+31,(c5+15)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 15), -(c5 + 15) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+16)+1);
                                    ubv=min(32*c4+31,(c5+16)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 16), -(c5 + 16) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+17)+1);
                                    ubv=min(32*c4+31,(c5+17)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 17), -(c5 + 17) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+18)+1);
                                    ubv=min(32*c4+31,(c5+18)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 18), -(c5 + 18) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+19)+1);
                                    ubv=min(32*c4+31,(c5+19)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 19), -(c5 + 19) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+20)+1);
                                    ubv=min(32*c4+31,(c5+20)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 20), -(c5 + 20) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+21)+1);
                                    ubv=min(32*c4+31,(c5+21)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 21), -(c5 + 21) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+22)+1);
                                    ubv=min(32*c4+31,(c5+22)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 22), -(c5 + 22) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+23)+1);
                                    ubv=min(32*c4+31,(c5+23)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 23), -(c5 + 23) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+24)+1);
                                    ubv=min(32*c4+31,(c5+24)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 24), -(c5 + 24) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+25)+1);
                                    ubv=min(32*c4+31,(c5+25)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 25), -(c5 + 25) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+26)+1);
                                    ubv=min(32*c4+31,(c5+26)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 26), -(c5 + 26) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+27)+1);
                                    ubv=min(32*c4+31,(c5+27)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 27), -(c5 + 27) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+28)+1);
                                    ubv=min(32*c4+31,(c5+28)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 28), -(c5 + 28) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+29)+1);
                                    ubv=min(32*c4+31,(c5+29)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 29), -(c5 + 29) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+30)+1);
                                    ubv=min(32*c4+31,(c5+30)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 30), -(c5 + 30) + c6);
                                    }
                                }

                                {
                                    lbv=max(32*c4,(c5+31)+1);
                                    ubv=min(32*c4+31,(c5+31)+N-2);
#pragma ivdep
#pragma vector always
                                    for (c6=lbv; c6<=ubv; c6++) {
                                        S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, (c5 + 31), -(c5 + 31) + c6);
                                    }
                                }
                            }

                            for (; c5 <= min(min(T - 1, 32 * c3 + 31), 32 * c4 + 30); c5 = c5 + 1) {
                                lbv=max(32*c4,c5+1);
                                ubv=min(32*c4+31,c5+N-2);
#pragma ivdep
#pragma vector always
                                for (c6=lbv; c6<=ubv; c6++) {
                                    S1(c1 - c2, -c1 + 2 * c2, c3, -c3 + c4, c5, -c5 + c6);
                                }
                            }
                        }
                        /*@ end @*/
                    }
                }
            }
        }
        /* End of CLooG code */

        orio_t_end = rtclock();
        orio_t_total += orio_t_end - orio_t_start;
    }

    orio_t_total = orio_t_total / REPS;
    printf("%f\n", orio_t_total);

    return 0;
}
