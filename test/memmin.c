constant N;

DO i = 1, N a[i] = c[i];
e[i] = d[i];
end do

    DO i = 1,
       N - 2 a[i] = 0.33 * (c[i - 1] + c[i] + c[i + 1]);
END DO

    DO i = 2,
       N - 3 d[i] = 0.33 * (a[i - 1] + a[i] + a[i + 1]);
END DO
