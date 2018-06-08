
/* pluto start (T,N) */

for (t=1; t<=T-1; t++)
    for (i=1; i<=N-2; i++)
        a[t][i] = 0.333 * (a[t-1][i-1] + a[t-1][i] + a[t-1][i+1]);

/* pluto end */
