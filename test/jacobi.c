
/* pluto start (T,N) */

DO t=1,T-1
  DO i=1,N-2
      u[t][i] = 0.333*(u[t-1][i-1] + u[t-1][i] + u[t-1][i+1]);
  END DO
END DO

/* pluto end */
