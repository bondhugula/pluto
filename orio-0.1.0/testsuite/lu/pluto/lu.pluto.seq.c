
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
double A[N][N+13];

void init_arrays() {
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

        register int i,j,k;




#define S1(zT0,zT1,zT2,zT3,k,j)	{A[k][j]=A[k][j]/A[k][k];}
#define S2(zT0,zT1,zT2,zT3,zT4,zT5,k,i,j)	{A[i][j]=A[i][j]-A[i][k]*A[k][j];}

        int c1, c2, c3, c4, c5, c6, c7, c8, c9;

        register int lbv, ubv;

        /* Generated from PLuTo-produced CLooG file by CLooG v0.14.1 64 bits in 1.93s. */
        for (c1=0; c1<=floord(N-2,256); c1++) {
            for (c2=max(0,ceild(128*c1-127,128)); c2<=floord(N-1,256); c2++) {
                for (c3=max(ceild(128*c1-32385,32640),ceild(128*c1-127,128)); c3<=floord(N-1,256); c3++) {
                    for (c4=max(max(8*c1,8*c1-1792*c3-1778),0); c4<=min(min(min(min(floord(128*c2+127,16),floord(3968*c3+3937,16)),floord(128*c3+127,16)),8*c1+7),floord(N-2,32)); c4++) {
                        for (c5=max(max(ceild(16*c4-15,16),0),8*c2); c5<=min(8*c2+7,floord(N-1,32)); c5++) {
                            for (c6=max(max(max(max(ceild(16*c4-465,496),ceild(8*c1-8*c3-c4-217,223)),ceild(-8*c1+8*c3+c4-217,225)),8*c3),ceild(16*c4-15,16)); c6<=min(8*c3+7,floord(N-1,32)); c6++) {
                                if ((c1 == c3) && (c4 == c6)) {
                                    for (c7=max(32*c6,0); c7<=min(min(32*c6+30,N-2),32*c5+30); c7++) {
                                        for (c8=max(c7+1,32*c5); c8<=min(32*c5+31,N-1); c8++) {
                                            S1(c1,c2,c4,c5,c7,c8) ;
                                            for (c9=c7+1; c9<=min(32*c6+31,N-1); c9++) {
                                                S2(c1,c1,c2,c4,c4,c5,c7,c9,c8) ;
                                            }
                                        }
                                    }
                                }
                                for (c7=max(0,32*c4); c7<=min(min(32*c6-1,32*c4+31),32*c5+30); c7++) {
                                    /*@ begin Loop(
                                    transform UnrollJam(ufactor=8)
                                                  for (c8=max(32*c5,c7+1);c8<=min(N-1,32*c5+31);c8++)
                                    transform Unroll(ufactor=8)
                                                    for (c9=32*c6;c9<=min(N-1,32*c6+31);c9++)
                                    {
                                                      S2(c1,c3,c2,c4,c6,c5,c7,c9,c8) ;
                                    }
                                    ) @*/{

                                        for (c8 = max(32 * c5, c7 + 1); c8 <= min(N - 1, 32 * c5 + 31) - 7; c8 = c8 + 8)     {

                                            for (c9 = 32 * c6; c9 <= min(N - 1, 32 * c6 + 31) - 7; c9 = c9 + 8)         {
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), (c8 + 7));
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), (c8 + 7));
                                            }

                                            for (; c9 <= min(N - 1, 32 * c6 + 31); c9 = c9 + 1)         {
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 1));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 2));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 3));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 4));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 5));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 6));
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, (c8 + 7));
                                            }
                                        }

                                        for (; c8 <= min(N - 1, 32 * c5 + 31); c8 = c8 + 1)     {

                                            for (c9 = 32 * c6; c9 <= min(N - 1, 32 * c6 + 31) - 7; c9 = c9 + 8)         {
                                                S2(c1, c3, c2, c4, c6, c5, c7, c9, c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 1), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 2), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 3), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 4), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 5), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 6), c8);
                                                S2(c1, c3, c2, c4, c6, c5, c7, (c9 + 7), c8);
                                            }

                                            for (; c9 <= min(N - 1, 32 * c6 + 31); c9 = c9 + 1)         S2(c1, c3, c2, c4, c6, c5, c7, c9, c8);
                                        }
                                    }
                                    /*@ end @*/
                                }
                                if ((c1 == c3) && (-c4 == -c6) && (c4 <= min(floord(N-33,32),floord(32*c5-1,32)))) {
                                    for (c8=max(32*c5,32*c4+32); c8<=min(N-1,32*c5+31); c8++) {
                                        S1(c1,c2,c4,c5,32*c4+31,c8) ;
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

    return ((int) A[0][0]);

}









