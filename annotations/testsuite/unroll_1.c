void ten_reciprocal_roots(double* x, double* f) {
    int i;
    /*@ begin Loop(
      transform Unroll(ufactor=8)
      for (i = 0; i <= N-1; i++)
        f[i] = 1.0 / sqrt(x[i]);
    ) @*/
    for (i = 0; i <= N-1; i++)
        f[i] = 1.0 / sqrt(x[i]);
    /*@ end @*/
}

void matrix_add(double** x, double** a, double** b, int m, int n) {
    int i, j;
    /*@ begin Loop(
      transform Unroll(ufactor=3)
      for (i = 0; i <= m-1; i++)
        transform Unroll(ufactor=3)
        for (j = 0; j <= n-1; j++)
          S(i,j);
    ) @*/
    for (i = 0; i <= m-1; i++)
        for (j = 0; j <= n-1; j++)
            S(i,j);
    /*@ end @*/
}


