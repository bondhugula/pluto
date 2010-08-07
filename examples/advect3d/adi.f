constant N;

DO t = 1, N

DO j = 2, N
DO i = 1, N
        c(i, j) = c(i, j) - c(i, j - 1) * a(i, j) / b(i, j - 1)
        b(i, j) = b(i, j) - a(i, j) * a(i, j) / b(i, j - 1)
END DO
END DO

DO i = 1, N 
c(i, N) = c(i, N) / b(i, N)
END DO

DO j = 1, N-1
DO i = 2, N
c(i, N-j) = ( c(i, N-j) - a(i, N-j + 1) * c(i, N-j + 1) ) / b(i, N-j)
END DO
END DO


DO j = 1, N
DO i = 2, N
c(i, j) = c(i, j) - c(i - 1, j) * a(i, j) / b(i - 1, j)
b(i, j) = b(i, j) - a(i, j) * a(i, j) / b(i - 1, j)
END DO
END DO

DO j = 1, N
c(N, j) = c(N, j) / b(N, j)
END DO

DO j = 1, N 
DO i = 1,N-1
c(N-i, j) = ( c(N-i, j) - a(N-i + 1, j) * c(N-i + 1, j) ) / b(N-i, j)
END DO
END DO

END DO
