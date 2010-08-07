! Generated from ../../../git/cloog/test/dartef.cloog by CLooG 0.14.0-136-gb91ef26 gmp bits in 0.29s.
IF (n >= 1) THEN
  DO t3=n+3, 3*n+1
    IF (MOD(t3+n+1, 2) == 0) THEN
      S1(1,n,(t3-n-1)/2)
    END IF
  END DO
  IF ((n <= 2) .AND. (n >= 2)) THEN
    DO t2=-n+4, 3*n-2
      DO t3=t2+2, t2+2*n
        IF (MOD(t2+n, 2) == 0) THEN
          IF (MOD(t3+n, 2) == 0) THEN
            S1((t2-n+2)/2,(t2+n-2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t3=1, n
      S2(1,n,t3)
    END DO
  END IF
  IF (n >= 3) THEN
    DO t2=n, n+2
      DO t3=t2+2, t2+2*n
        IF (MOD(t2+n, 2) == 0) THEN
          IF (MOD(t3+n, 2) == 0) THEN
            S1((t2-n+2)/2,(t2+n-2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t3=1, n
      S2(1,n,t3)
    END DO
  END IF
  DO t1=CEILING(REAL(-2*n+5)/REAL(2)), MIN(-1,-n+6)
    DO t2=-t1+2, -t1+4
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+5, t1+2*n
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t3=1, n
      S2(t1+n-1,n,t3)
    END DO
  END DO
  IF (n == 2) THEN
    DO t3=5, 7
      IF (MOD(t3+1, 2) == 0) THEN
        S1(2,1,(t3-3)/2)
      END IF
    END DO
    DO t2=4, 6
      DO t3=1, 2
        IF (MOD(t2, 2) == 0) THEN
          S2((t2-2)/2,(t2-2)/2,t3)
        END IF
      END DO
    END DO
  END IF
  DO t1=-n+7, -1
    DO t2=-t1+2, -t1+4
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+5, n-2
      DO t3=1, t2+1
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
      DO t3=n+1, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=n-1, t1+2*n
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t3=1, n
      S2(t1+n-1,n,t3)
    END DO
  END DO
  IF (n >= 3) THEN
    DO t1=0, MIN(1,-n+6)
      DO t2=t1+2, -t1+4
        DO t3=t2+2, t2+2*n
          IF (MOD(t1+t2, 2) == 0) THEN
            IF (MOD(t1+t3, 2) == 0) THEN
              S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
            END IF
          END IF
        END DO
      END DO
      DO t2=-t1+5, -t1+2*n
        DO t3=1, n
          IF (MOD(t1+t2+1, 2) == 0) THEN
            S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
          END IF
        END DO
        DO t3=t2+2, t2+2*n
          IF (MOD(t1+t2, 2) == 0) THEN
            IF (MOD(t1+t3, 2) == 0) THEN
              S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
            END IF
          END IF
        END DO
      END DO
      DO t2=-t1+2*n+1, t1+2*n+1
        DO t3=1, n
          IF (MOD(t1+t2+1, 2) == 0) THEN
            S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
          END IF
        END DO
      END DO
    END DO
  END IF
  DO t1=MAX(0,-n+7), 1
    DO t2=t1+2, -t1+4
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+5, n-2
      DO t3=1, t2+1
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
      DO t3=n+1, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=n-1, -t1+2*n
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+2*n+1, t1+2*n+1
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
    END DO
  END DO
  DO t1=2, n-5
    DO t3=t1+4, t1+2*n+2
      IF (MOD(t1+t3, 2) == 0) THEN
        S1(t1+1,1,(-t1+t3-2)/2)
      END IF
    END DO
    DO t2=t1+3, n-2
      DO t3=1, t2+1
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
      DO t3=n+1, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=n-1, -t1+2*n
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+2*n+1, -t1+2*n+3
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
    END DO
  END DO
  DO t1=MAX(2,n-4), FLOOR(REAL(2*n-3)/REAL(2))
    DO t3=t1+4, t1+2*n+2
      IF (MOD(t1+t3, 2) == 0) THEN
        S1(t1+1,1,(-t1+t3-2)/2)
      END IF
    END DO
    DO t2=t1+3, -t1+2*n
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
      DO t3=t2+2, t2+2*n
        IF (MOD(t1+t2, 2) == 0) THEN
          IF (MOD(t1+t3, 2) == 0) THEN
            S1((t1+t2)/2,(-t1+t2)/2,(-t2+t3)/2)
          END IF
        END IF
      END DO
    END DO
    DO t2=-t1+2*n+1, -t1+2*n+3
      DO t3=1, n
        IF (MOD(t1+t2+1, 2) == 0) THEN
          S2((t1+t2-3)/2,(-t1+t2-1)/2,t3)
        END IF
      END DO
    END DO
  END DO
  IF (n >= 3) THEN
    DO t3=n+3, 3*n+1
      IF (MOD(t3+n+1, 2) == 0) THEN
        S1(n,1,(t3-n-1)/2)
      END IF
    END DO
    DO t2=n+2, n+4
      DO t3=1, n
        IF (MOD(t2+n, 2) == 0) THEN
          S2((t2+n-4)/2,(t2-n)/2,t3)
        END IF
      END DO
    END DO
  END IF
  DO t3=1, n
    S2(n,1,t3)
  END DO
END IF
