! Generated from ../../../git/cloog/test/levenshtein-1-2-3f.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.11s.
S1(0,0)
S2(1,0)
S3(1,1)
DO i=2, N
  S2(i,0)
  DO j=1, i-1
    S6(i,j)
  END DO
  S3(i,i)
END DO
S7(N+1,0)
DO j=1, N
  S6(N+1,j)
  S8(N+1,j)
END DO
DO i=N+2, 2*M-N-2
  j = FLOOR(REAL(i-N-1)/REAL(2))
  S7(i,j)
  IF (MOD(i+N, 2) == 0) THEN
    S5(i,(i-N)/2)
    S8(i,(i-N)/2)
  END IF
  DO j=CEILING(REAL(i-N+1)/REAL(2)), FLOOR(REAL(i+N-1)/REAL(2))
    S6(i,j)
    S8(i,j)
  END DO
  IF (MOD(i+N, 2) == 0) THEN
    S4(i,(i+N)/2)
    S8(i,(i+N)/2)
  END IF
END DO
DO i=2*M-N-1, 2*M-2
  DO j=i-M+1, M-1
    S6(i,j)
  END DO
END DO
