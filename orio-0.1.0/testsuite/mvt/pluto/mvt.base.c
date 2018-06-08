

#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define M MSIZE
#define N NSIZE
double a[M][N];
double y_1[N];
double y_2[M];
double x1[M];
double x2[N];

void init_arrays() {
    int i1, i2;
    for (i1=0; i1<M; i1++)
        for (i2=0; i2<N; i2++)
            a[i1][i2] = (i1+i2) % 5 + 1;
    for (i1=0; i1<N; i1++)
        y_1[i1] = (i1) % 5 + 1;
    for (i1=0; i1<M; i1++)
        y_2[i1] = (i1) % 5 + 1;
    for (i1=0; i1<M; i1++)
        x1[i1] = 0;
    for (i1=0; i1<N; i1++)
        x2[i1] = 0;
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

        int i, j;
        int ii, jj;
        int iii, jjj;


        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                x1[i] = x1[i] + a[i][j] * y_1[j];
            }
        }

        for (i=0; i<N; i++) {
            for (j=0; j<N; j++) {
                x2[i] = x2[i] + a[j][i] * y_2[j];
            }
        }

        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    return 1;
}

