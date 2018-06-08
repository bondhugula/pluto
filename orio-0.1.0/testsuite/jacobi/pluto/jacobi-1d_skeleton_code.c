#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/*@ global @*/

double rtclock() {
    struct timezone tzp;
    struct timeval tp;
    int stat;
    gettimeofday (&tp, &tzp);
    return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main() {
    /*@ prologue @*/

    double orio_t_start=0, orio_t_end=0, orio_t_total=0;
    int orio_i;

    for (orio_i=0; orio_i<REPS; orio_i++) {
        orio_t_start = rtclock();

        /*@ tested code @*/

        orio_t_end = rtclock();
        orio_t_total += orio_t_end - orio_t_start;
    }

    orio_t_total = orio_t_total / REPS;
    printf("%f\n", orio_t_total);

    /*@ epilogue @*/

    return a[0][0];
}
