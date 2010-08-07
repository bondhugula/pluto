CONSTANT n, T;

DO t=1, T
  DO i=1,n
    DO j=1,n
      u[t,i,j] = u[t-1,i-1,j] + u[t-1,i,j+1] + u[t-1,i+1,j] + u[t-1,i,j-1]
    END DO
  END DO
END DO
