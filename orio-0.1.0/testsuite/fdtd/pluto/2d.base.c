

#include <stdio.h>
#include <sys/time.h>
#include <math.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define tmax TVAL
#define nx NXVAL
#define ny NYVAL
double ex[nx][ny+1];
double ey[nx+1][ny];
double hz[nx][ny];

void init_arrays() {
    int i1, i2;
    for (i1=0; i1<nx; i1++)
        for (i2=0; i2<ny+1; i2++)
            ex[i1][i2] = (i1+i2) % 5 + 1;
    for (i1=0; i1<nx+1; i1++)
        for (i2=0; i2<ny; i2++)
            ey[i1][i2] = (i1+i2) % 5 + 1;
    for (i1=0; i1<nx; i1++)
        for (i2=0; i2<ny; i2++)
            hz[i1][i2] = (i1+i2) % 5 + 1;
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

        int t, i, j, k, l, m, n;

        for(t=0; t<tmax; t++)  {
            for (j=0; j<ny; j++)
                ey[0][j] = t;
            for (i=1; i<nx; i++)
                for (j=0; j<ny; j++)
                    ey[i][j] = ey[i][j] - 0.5*(hz[i][j]-hz[i-1][j]);
            for (i=0; i<nx; i++)
                for (j=1; j<ny; j++)
                    ex[i][j] = ex[i][j] - 0.5*(hz[i][j]-hz[i][j-1]);
            for (i=0; i<nx; i++)
                for (j=0; j<ny; j++)
                    hz[i][j]=hz[i][j]-0.7*(ex[i][j+1]-ex[i][j]+ey[i+1][j]-ey[i][j]);
        }

        annot_t_end = rtclock();
        annot_t_total += annot_t_end - annot_t_start;
    }

    annot_t_total = annot_t_total / REPS;
    printf("%f\n", annot_t_total);

    return 1;
}

