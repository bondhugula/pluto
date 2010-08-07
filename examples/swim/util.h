#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void print()
{
    int i, j;

    for (i=0; i<N1-1; i++) {
        for (j=0; j<N2-1; j++) {
            fprintf(stdout, "%lf ", (UNEW[i][j]));
            fprintf(stdout, "%lf ", (VNEW[i][j]));
            fprintf(stdout, "%lf ", (PNEW[i][j]));
            if (j%80 == 20) fprintf(stdout, "\n");
        }
    }
    fprintf(stdout, "\n");
}

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}
