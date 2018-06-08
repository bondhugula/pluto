void init_input_vars() {
    int i, j;
    for (i=0; i<T; i++)
        for (j=0; j<N; j++)
            a[i][j] = i+((double)j)/N;
}
