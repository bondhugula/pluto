
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

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

void print_array() {
    int i, j;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fprintf(stderr, "%lf ", C[i][j]);
            if (j%80 == 79) fprintf(stderr, "\n");
        }
        fprintf(stderr, "\n");
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

        int i, j, k;

        for (i=0; i<=M-1; i++ )
            for (j=0; j<=N-1; j++ )
                for (k=0; k<=K-1; k++ )
                    C[i][j]=C[i][j]+A[i][k]*B[k][j];

        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    //print_array();

    return 1;
}

