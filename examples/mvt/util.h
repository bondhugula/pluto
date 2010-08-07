#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

void init_array()
{
    int i, j;

    for (i=0; i<N; i++) {
        y_1[i] = i;
        y_2[i] = i+1;
        x1[i] = 0.0;
        x2[i] = 0.0;

        for (j=0; j<N; j++)
            a[i][j] = i+j+1.0;
    }
}

void print_array()
{
    int i, j;

    for (i=0; i<N; i++) {
        fprintf(stdout, "%lf ", x1[i]);
        if (j%80 == 20) fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
    for (i=0; i<N; i++) {
        fprintf(stdout, "%lf ", x2[i]);
        if (j%80 == 20) fprintf(stdout, "\n");
    }
    fprintf(stdout, "\n");
}
