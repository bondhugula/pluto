
for (i=0; i<=N-1; i++) {
    for (j=0; j<=N-1; j++) {
        x1[i]=x1[i]+a[i][j]*y_1[j];
        x2[j]=x2[j]+a[i][j]*y_2[i];
    }
}

