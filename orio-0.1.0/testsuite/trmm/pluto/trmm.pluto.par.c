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


#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define S1(zT0,zT1,zT2,zT3,zT4,zT5,i,j,k)	{B[i][j]=alpha*A[i][k]*B[k][j]+B[i][j];}

        int c1, c2, c3, c4, c5, c6, c7, c8, c9;

        register int lb, ub, lb1, ub1, lb2, ub2;
        register int lbv, ubv;

        /* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 0.91s. */
        lb1=0;
        ub1=floord(N-1,256);
        #pragma omp parallel for shared(lb1,ub1) private(c1,c2,c3,c4,c5,c6,c7,c8,c9)
        for (c1=lb1; c1<=ub1; c1++) {
            for (c2=0; c2<=floord(N-2,256); c2++) {
                for (c3=max(ceild(128*c2-127,128),0); c3<=floord(N-1,256); c3++) {
                    for (c4=max(0,8*c1); c4<=min(8*c1+7,floord(N-1,32)); c4++) {
                        for (c5=max(8*c2,0); c5<=min(min(8*c2+7,floord(N-2,32)),floord(128*c3+127,16)); c5++) {
                            for (c6=max(max(8*c3,ceild(16*c5-15,16)),0); c6<=min(8*c3+7,floord(N-1,32)); c6++) {
                                /*@ begin Loop(
                                transform UnrollJam(ufactor=8)
                                            for (c7=max(max(32*c6,1),32*c5+1);c7<=min(N-1,32*c6+31);c7++)
                                transform UnrollJam(ufactor=8)
                                              for (c8=max(32*c5,0);c8<=min(c7-1,32*c5+31);c8++)
                                {
                                {
                                	lbv=max(32*c4,0);
                                	ubv=min(N-1,32*c4+31);
                                #pragma ivdep
                                #pragma vector always
                                	for (c9=lbv; c9<=ubv; c9++) {
                                                  S1(c2,c1,c3,c5,c4,c6,c8,c9,c7) ;
                                                }
                                }
                                }
                                ) @*/{

                                    for (c7 = max(max(32 * c6, 1), 32 * c5 + 1); c7 <= min(N - 1, 32 * c6 + 31) - 7; c7 = c7 + 8)     {

                                        for (c8 = max(32 * c5, 0); c8 <= min(c7 - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, c7);
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 1) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 1));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 1));
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 2) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 2));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 2));
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 3) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 3));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 3));
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 4) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 4));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 4));
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 5) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 5));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 5));
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 6) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 6));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 6));
                                            }
                                        }

                                        for (c8 = max(32 * c5, 0); c8 <= min((c7 + 7) - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, (c7 + 7));
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, (c7 + 7));
                                            }
                                        }

                                        for (; c8 <= min(c7 - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, c7);
                                            }
                                        }

                                        for (; c8 <= min((c7 + 1) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 1));
                                            }
                                        }

                                        for (; c8 <= min((c7 + 2) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 2));
                                            }
                                        }

                                        for (; c8 <= min((c7 + 3) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 3));
                                            }
                                        }

                                        for (; c8 <= min((c7 + 4) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 4));
                                            }
                                        }

                                        for (; c8 <= min((c7 + 5) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 5));
                                            }
                                        }

                                        for (; c8 <= min((c7 + 6) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 6));
                                            }
                                        }

                                        for (; c8 <= min((c7 + 7) - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, (c7 + 7));
                                            }
                                        }
                                    }

                                    for (; c7 <= min(N - 1, 32 * c6 + 31); c7 = c7 + 1)     {

                                        for (c8 = max(32 * c5, 0); c8 <= min(c7 - 1, 32 * c5 + 31) - 7; c8 = c8 + 8) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 1), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 2), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 3), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 4), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 5), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 6), c9, c7);
                                                S1(c2, c1, c3, c5, c4, c6, (c8 + 7), c9, c7);
                                            }
                                        }

                                        for (; c8 <= min(c7 - 1, 32 * c5 + 31); c8 = c8 + 1) {
                                            lbv=max(32*c4,0);
                                            ubv=min(N-1,32*c4+31);
#pragma ivdep
#pragma vector always
                                            for (c9=lbv; c9<=ubv; c9++) {
                                                S1(c2, c1, c3, c5, c4, c6, c8, c9, c7);
                                            }
                                        }
                                    }
                                }
                                /*@ end @*/
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

    return ((int) B[0][0]);
}



