

#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define T TVAL
#define N NVAL
#define coeff1 0.5
#define coeff2 0.7
double h[N];
double e[N+1];

void init_arrays() {
    int i1;
    for (i1=0; i1<N; i1++)
        h[i1] = (i1) % 5 + 1;
    for (i1=0; i1<N+1; i1++)
        e[i1] = (i1) % 5 + 1;
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



        int t, i, j, k, l;

        for (t=1; t<=T; t++) {
            for (i=1; i<=N-1; i++)
                e[i] = e[i] - coeff1*(h[i]-h[i-1]);
            for (i=0; i<=N-1; i++)
                h[i] = h[i] - coeff2*(e[i+1]-e[i]);
        }




        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    return (e[0] + h[0]);
}

