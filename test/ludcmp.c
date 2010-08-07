#include <stdio.h>


double t_start, t_end;

#include "util.h"

#define n 2000

float a[n+1][n+1];
float x[n+1];
float y[n+1];
float b[n+1];
float w;

int main()
{
    int i, j, k;


    b[0] = 1.0;

    IF_TIME(t_start = rtclock());

    /* pluto start (n) */

    for (i = 0; i < n; i++)
    {
        for (j = i+1; j <= n; j++)
        {
            w = a[j][i];
            for (k = 0; k < i; k++)
                w = w- a[j][k] * a[k][i];
            a[j][i] = w / a[i][i];
        }
        for (j = i+1; j <= n; j++)
        {
            w = a[i+1][j];
            for (k = 0; k <= i; k++)
                w = w  - a[i+1][k] * a[k][j];
            a[i+1][j] = w;
        }
    }
    y[0] = b[0];
    for (i = 1; i <= n; i++)
    {
        w = b[i];
        for (j = 0; j < i; j++)
            w = w - a[i][j] * y[j];
        y[i] = w;
    }
    x[n] = y[n] / a[n][n];
    for (i = 0; i <= n - 1; i++)
    {
        w = y[n - 1 - (i)];
        for (j = n - i; j <= n; j++)
            w = w- a[n - 1 - i][j] * x[j];
        x[n - 1 - i] = w / a[n - 1 - (i)][n - 1-(i)];
    }

    /* pluto end */

    IF_TIME(t_end = rtclock());
    IF_TIME(printf("%0.6lfs\n", t_end - t_start));

}
