! Generated from ./test/dartef.cloog by CLooG v1.0.0 64 bits in 0.32s.
! CLooG asked for 164 KBytes.
IF (n >= 1) THEN
  t1 = -n+1
  t2 = n+1
  DO t3=n+3, 3*n+1
    IF (MOD(t3-n-1,2) == 0) THEN
      k = (t3-n-1)/2
      S1(i = 1,j = n)
    END IF
  END DO
END IF
IF ((n >= 2) .AND. (n <= 2)) THEN
  t1 = -n+2
  DO t2=-n+4, 3*n-2
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  t2 = n+3
  DO t3=1, n
    IF (MOD(t1+t2-3,2) == 0) THEN
      i = (t1+t2-3)/2
      j = -t1+i+1
      S2(k = t3)
    END IF
  END DO
END IF
IF (n >= 3) THEN
  t1 = -n+2
  DO t2=n, n+2
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  t2 = n+3
  DO t3=1, n
    i = t1+n-1
    S2(j = n,k = t3)
  END DO
END IF
DO t1=CEILING(REAL(-2*n+5)/REAL(2)), MIN(-n+6,-1)
  DO t2=-t1+2, -t1+4
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=-t1+5, t1+2*n
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  t2 = t1+2*n+1
  DO t3=1, n
    i = t1+n-1
    S2(j = n,k = t3)
  END DO
END DO
IF (n == 2) THEN
  DO t3=3+2, 3+2*n
    IF (MOD(1+3,2) == 0) THEN
      i = (1+3)/2
      j = -(1)+i
      IF (MOD(t3-i-j,2) == 0) THEN
        k = (t3-i-j)/2
        S1
      END IF
    END IF
  END DO
  DO t2=4, 6
    DO t3=1, n
      IF (MOD(1+t2-3,2) == 0) THEN
        i = (1+t2-3)/2
        j = -(1)+i+1
        S2(k = t3)
      END IF
    END DO
  END DO
END IF
DO t1=-n+7, -1
  DO t2=-t1+2, -t1+4
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=-t1+5, n-2
    DO t3=1, t2+1
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=n+1, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
    DO t3=t2+2, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=n-1, t1+2*n
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  t2 = t1+2*n+1
  DO t3=1, n
    i = t1+n-1
    S2(j = n,k = t3)
  END DO
END DO
IF (n >= 3) THEN
  DO t1=0, MIN(1,-n+6)
    DO t2=t1+2, -t1+4
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2,2) == 0) THEN
          i = (t1+t2)/2
          j = -t1+i
          IF (MOD(t3-i-j,2) == 0) THEN
            k = (t3-i-j)/2
            S1
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+5, -t1+2*n
      DO t3=1, n
        IF (MOD(t1+t2-3,2) == 0) THEN
          i = (t1+t2-3)/2
          j = -t1+i+1
          S2(k = t3)
        END IF
      END DO
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2,2) == 0) THEN
          i = (t1+t2)/2
          j = -t1+i
          IF (MOD(t3-i-j,2) == 0) THEN
            k = (t3-i-j)/2
            S1
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+2*n+1, t1+2*n+1
      DO t3=1, n
        IF (MOD(t1+t2-3,2) == 0) THEN
          i = (t1+t2-3)/2
          j = -t1+i+1
          S2(k = t3)
        END IF
      END DO
    END DO
  END DO
END IF
DO t1=MAX(-n+7,0), 1
  DO t2=t1+2, -t1+4
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=-t1+5, n-2
    DO t3=1, t2+1
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=n+1, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
    DO t3=t2+2, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=n-1, -t1+2*n
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=-t1+2*n+1, t1+2*n+1
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
  END DO
END DO
DO t1=2, n-5
  t2 = t1+2
  DO t3=t2+2, t2+2*n
    i = t1+1
    IF (MOD(t3-i-1,2) == 0) THEN
      k = (t3-i-1)/2
      S1(j = 1)
    END IF
  END DO
  DO t2=t1+3, n-2
    DO t3=1, t2+1
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=n+1, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
    DO t3=t2+2, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=n-1, -t1+2*n
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=-t1+2*n+1, -t1+2*n+3
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
  END DO
END DO
DO t1=MAX(2,n-4), FLOOR(REAL(2*n-3)/REAL(2))
  t2 = t1+2
  DO t3=t2+2, t2+2*n
    i = t1+1
    IF (MOD(t3-i-1,2) == 0) THEN
      k = (t3-i-1)/2
      S1(j = 1)
    END IF
  END DO
  DO t2=t1+3, -t1+2*n
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
    DO t3=t2+2, t2+2*n
      IF (MOD(t1+t2,2) == 0) THEN
        i = (t1+t2)/2
        j = -t1+i
        IF (MOD(t3-i-j,2) == 0) THEN
          k = (t3-i-j)/2
          S1
        END IF
      END IF
    END DO
  END DO
  DO t2=-t1+2*n+1, -t1+2*n+3
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
  END DO
END DO
IF (n >= 3) THEN
  t1 = n-1
  t2 = n+1
  DO t3=t2+2, t2+2*n
    i = t1+1
    IF (MOD(t3-i-1,2) == 0) THEN
      k = (t3-i-1)/2
      S1(j = 1)
    END IF
  END DO
  DO t2=n+2, n+4
    DO t3=1, n
      IF (MOD(t1+t2-3,2) == 0) THEN
        i = (t1+t2-3)/2
        j = -t1+i+1
        S2(k = t3)
      END IF
    END DO
  END DO
END IF
IF (n >= 1) THEN
  t2 = n+3
  DO t3=1, n
    S2(i = n,j = 1,k = t3)
  END DO
END IF
