#pragma scop
for (i=0; i<n; i++)
    for (j=2; j<i; j++)
        a[i] = a[i] + 1;

#pragma endscop

/* pluto end */
