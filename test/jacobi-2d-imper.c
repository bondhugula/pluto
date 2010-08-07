

#pragma scop
    for (t=0; t<T; t++) {
        for (i=2; i<N-1; i++) {
            for (j=2; j<N-1; j++) {
                b[i][j]= 0.2*(a[i][j]+a[i][j-1]+a[i][1+j]+a[1+i][j]+a[i-1][j]);
            }
        }
        for (i=2; i<N-1; i++) {
            for (j=2; j<N-1; j++)   {
                a[i][j]=b[i][j];
            }
        }
    }

#pragma endscop

