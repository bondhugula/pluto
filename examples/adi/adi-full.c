#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TMAX 1000
#define NMAX 1000

static double X[NMAX][NMAX], A[NMAX][NMAX], B[NMAX][NMAX];

void adi(long T, long N) {

    int t, i1, i2;

    for (t = 0; t < T; t++) {

        // Column Sweep
        for (i1=0; i1<N; i1++) {
            for (i2 = 1; i2 < N; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1][i2-1] * A[i1][i2] / B[i1][i2-1];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1][i2-1];
            }
        }

        // Row Sweep
        for (i1=1; i1<N; i1++) {
            for (i2 = 0; i2 < N; i2++) {
                X[i1][i2] = X[i1][i2] - X[i1-1][i2] * A[i1][i2] / B[i1-1][i2];
                B[i1][i2] = B[i1][i2] - A[i1][i2] * A[i1][i2] / B[i1-1][i2];
            }
        }

    } // end of t loop
}

int main()
{
    long T, N;
    T = TMAX; N=NMAX;
    adi(T,N);

    return 0;
}
