#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double L[N][N];
double U[N][N];

void init_array()
{
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
                a[i][j] += L[i][k]*U[k][j];//i==j?1:0;
            }
        }
    }
}


void print_array()
{
    int i, j;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            fprintf(stdout, "%lf ", round(a[i][j]));
            if (j%80 == 79) fprintf(stdout, "\n");
        }
        fprintf(stdout, "\n");
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
