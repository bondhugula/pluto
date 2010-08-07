/* pluto start (n) */

DO i = 1, n
    DO j = 1, n
        y[i] = y[i] + a[i,j]*x[j];
    END DO
END DO

DO i = 1, n
    DO j = 1, n
        z[i] = z[i] + b[i,j]*y[j];
    END DO
END DO

/* pluto end */
