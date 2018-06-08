#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define V VSIZE
#define O OSIZE
double A2[V][O];
double T[V][O][O][O];
double R[V][V][O][O];

void init_arrays() {
    int i1, i2, i3, i4;
    for (i1=0; i1<V; i1++)
        for (i2=0; i2<O; i2++)
            A2[i1][i2] = (i1+i2) % 5 + 1;
    for (i1=0; i1<V; i1++)
        for (i2=0; i2<O; i2++)
            for (i3=0; i3<O; i3++)
                for (i4=0; i4<O; i4++)
                    T[i1][i2][i3][i4] = (i1+i2+i3+i4) % 5 + 1;
    for (i1=0; i1<V; i1++)
        for (i2=0; i2<V; i2++)
            for (i3=0; i3<O; i3++)
                for (i4=0; i4<O; i4++)
                    R[i1][i2][i3][i4] = 0;
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

    int v1,v2,o1,o2,ox;
    int tv1,tv2,to1,to2,tox;

    for (annot_i=0; annot_i<REPS; annot_i++) {
        annot_t_start = rtclock();

        for (v1=0; v1<=V-1; v1=v1+1)
            for (v2=0; v2<=V-1; v2=v2+1)
                for (o1=0; o1<=O-1; o1=o1+1)
                    for (o2=0; o2<=O-1; o2=o2+1)
                        for (ox=0; ox<=O-1; ox=ox+1)
                            R[v1][v2][o1][o2]=R[v1][v2][o1][o2]+T[v1][ox][o1][o2]*A2[v2][ox];

        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    return 1;
}

