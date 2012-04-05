#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <sys/time.h>
#include <assert.h>

#define nx 300
#define ny 300
#define nz 300
// #define T 20
// #define 9 9

#define f60 0.2 
#define f61 0.5
#define f62 0.3

#define halfdtbydx 0.5
#define thirddtbydz 0.3
#define thirddtbydx 0.3
#define thirddtbydy 0.3

double a[nx+10][ny+10][nz+10];
double af[nx+10][ny+10][nz+10];
double ab[nx+10][ny+10][nz+10];
double al[nx+10][ny+10][nz+10];
double athird[nx+10][ny+10][nz+10];
double uxl[nx+10][ny+10][nz+10];
double uzf[nx+10][ny+10][nz+10];
double uyb[nx+10][ny+10][nz+10];


#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{

}


void print_array()
{
    int i, j, k;

    for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
    for (k=0; k<nz; k++) {
        fprintf(stdout, "%lf ", athird[i][j][k]);
        if (j%80 == 20) fprintf(stdout, "\n");
    }
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

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

double t_start, t_end;

int main()
{
	assert(nx >= 100);
	assert(ny >= 100);
	assert(nz >= 100);
	int i, j, k, l, m, n, t;

	init_array() ;

#ifdef PERFCTR
	PERF_INIT; 
#endif

	IF_TIME(t_start = rtclock());
#define reps 10

#pragma scop
    for (j = 4; j <= ny+9-2; j++)
        for (i = 4; i <= nx+9-3; i++)
            for (k = 4; k <= nz+9-3; k++)
                ab[j][i][k] = (f60 * (a[j-1][i][k] + a[j][i][k]) + f61 
                        * (a[j-2][i][k] + a[j+1][i][k]) + f62 * (a[j-3][i][k] 
                            + a[j+2][i][k])) * thirddtbydy * uyb[j][i][k];

    for (j = 4; j <= ny+9-3; j++)
        for (i = 4; i <= nx+9-2; i++)
            for (k = 4; k <= nz+9-3; k++)
                al[j][i][k] = (f60 * (a[j][i-1][k] + a[j][i][k]) + f61 
                        * (a[j][i-2][k] + a[j][i+1][k]) + f62 * (a[j][i-3][k] + a[j][i+2][k]))
                    * thirddtbydx * uxl[j][i][k];

    for (j = 4; j <= ny+9-3; j++)
        for (i = 4; i <= nx+9-3; i++)
            for (k = 4; k <= nz+9-2; k++)
                af[j][i][k] = (f60 * (a[j][i][k-1] + a[j][i][k]) + f61 
                        * (a[j][i][k-2] + a[j][i][k+1]) + f62 * (a[j][i][k-3] + a[j][i][k+2]))
                    * thirddtbydz * uzf[j][i][k];


    for (j = 4; j <= ny+9-3; j++)
        for (i = 4; i <= nx+9-3; i++)
            for (k = 4; k <= nz+9-3; k++)
                athird[j][i][k] = a[j][i][k] + (al[j][i+1][k] - al[j][i][k])
                    + (ab[j+1][i][k] - ab[j][i][k]) + (af[j][i][k+1] - af[j][i][k]);
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef PERFCTR
    PERF_EXIT; 
#endif

#ifdef TEST
    print_array();
#endif
    return 0;
}
