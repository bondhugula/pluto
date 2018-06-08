
for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
        x1[i] = x1[i] + a[i][j] * y_1[j];
    }
}

for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
        x2[i] = x2[i] + a[j][i] * y_2[j];
    }
}


