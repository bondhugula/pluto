
/* pluto start (N) */
DO r = 1, N
    DO q = 1, N
        DO p = 1, N
            sum[r][q][p] = 0;
            DO s = 1, N
                sum[r][q][p] = sum[r][q][p] + A[r,q,s]*C4[p,s]
            END DO
        END DO
        DO p = 1, N
            a[r,q,p] = sum[r][q][p];
        END DO
    END DO
END DO

/* pluto end */
