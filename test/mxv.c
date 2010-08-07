CONSTANT n;

DO i = 1, n
    y[i] = 0;
    DO j = 1, n
        y[i] = y[i] + a[i,j]*x[j];
    END DO
END DO
