constant n;
/* pluto start (n) */

DO i=0, n
  A[2*i] = 42
END DO

DO j=0,n
  B[j] = A[3*j]
END DO
/* pluto end */
