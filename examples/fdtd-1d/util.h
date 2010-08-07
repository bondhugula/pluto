#include <unistd.h>
#include <sys/time.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{
    int i, j;

        for (j=0; j<N; j++) {
            h[j] = ((double)j)/N;
            e[j] = ((double)j)/N;
        }
}


void print_array()
{
    int i, j;

    for (j=0; j<N; j++) {
	    fprintf(stdout, "%lf ", h[j]);
	    if (j%80 == 79) fprintf(stdout, "\n");
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
