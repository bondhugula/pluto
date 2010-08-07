      IMPLICIT INTEGER	(I-N)
      IMPLICIT REAL*8	(A-H, O-Z)

      PARAMETER (N1=1335, N2=1335)

      COMMON  U(N1,N2), V(N1,N2), P(N1,N2),
     *        UNEW(N1,N2), VNEW(N1,N2),
     *        PNEW(N1,N2), UOLD(N1,N2),
     *        VOLD(N1,N2), POLD(N1,N2),
     *        CU(N1,N2), CV(N1,N2),
     *        Z(N1,N2), H(N1,N2), PSI(N1,N2)

      TIME = 0
      NCYCLE = 0
   90 NCYCLE = NCYCLE + 1
C     COMPUTE CAPITAL  U, CAPITAL V, Z AND H
c     CALL CALC1 (inlined)
      FSDX = 4.D0/DX
      FSDY = 4.D0/DY
      DO 100 J=1,N
      DO 100 I=1,M
      CU(I+1,J) = .5D0*(P(I+1,J)+P(I,J))*U(I+1,J)
      CV(I,J+1) = .5D0*(P(I,J+1)+P(I,J))*V(I,J+1)
      Z(I+1,J+1) = (FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1)
     1          -U(I+1,J)))/(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
      H(I,J) = P(I,J)+.25D0*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)
     1               +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
  100 CONTINUE
C     PERIODIC CONTINUATION
      DO 110 J=1,N
      CU(1,J) = CU(M+1,J)
      CV(M+1,J+1) = CV(1,J+1)
      Z(1,J+1) = Z(M+1,J+1)
      H(M+1,J) = H(1,J)
  110 CONTINUE
      DO 115 I=1,M
      CU(I+1,N+1) = CU(I+1,1)
      CV(I,1) = CV(I,N+1)
      Z(I+1,1) = Z(I+1,N+1)
      H(I,N+1) = H(I,1)
  115 CONTINUE
      CU(1,N+1) = CU(M+1,1)
      CV(M+1,1) = CV(1,N+1)
      Z(1,1) = Z(M+1,N+1)
      H(M+1,N+1) = H(1,1)
C     COMPUTE NEW VALUES U,V AND P
c     CALL CALC2 (inlined)
      TDTS8 = TDT/8.D0
      TDTSDX = TDT/DX
      TDTSDY = TDT/DY
      DO 200 J=1,N
      DO 200 I=1,M
      UNEW(I+1,J) = UOLD(I+1,J)+
     1    TDTS8*(Z(I+1,J+1)+Z(I+1,J))*(CV(I+1,J+1)+CV(I,J+1)+CV(I,J)
     2       +CV(I+1,J))-TDTSDX*(H(I+1,J)-H(I,J))
      VNEW(I,J+1) = VOLD(I,J+1)-TDTS8*(Z(I+1,J+1)+Z(I,J+1))
     1       *(CU(I+1,J+1)+CU(I,J+1)+CU(I,J)+CU(I+1,J))
     2       -TDTSDY*(H(I,J+1)-H(I,J))
      PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))
     1       -TDTSDY*(CV(I,J+1)-CV(I,J))
  200 CONTINUE
C     PERIODIC CONTINUATION
      DO 210 J=1,N
      UNEW(1,J) = UNEW(M+1,J)
      VNEW(M+1,J+1) = VNEW(1,J+1)
      PNEW(M+1,J) = PNEW(1,J)
  210 CONTINUE
      DO 215 I=1,M
      UNEW(I+1,N+1) = UNEW(I+1,1)
      VNEW(I,1) = VNEW(I,N+1)
      PNEW(I,N+1) = PNEW(I,1)
  215 CONTINUE
      UNEW(1,N+1) = UNEW(M+1,1)
      VNEW(M+1,1) = VNEW(1,N+1)
      PNEW(M+1,N+1) = PNEW(1,1)
      TIME = TIME + DT
C *** Later addition - we could remove these reductions if desired; does not
correspond to the original SWIM code, but was added to SPEC vesion
        PCHECK = 0.0D0
        UCHECK = 0.0D0
        VCHECK = 0.0D0
        
        DO 3500 ICHECK = 1, MNMIN
         DO 4500 JCHECK = 1, MNMIN
         PCHECK = PCHECK + ABS(PNEW(ICHECK,JCHECK))
         UCHECK = UCHECK + ABS(UNEW(ICHECK,JCHECK))
         VCHECK = VCHECK + ABS(VNEW(ICHECK,JCHECK))
 4500   CONTINUE
	UNEW(ICHECK,ICHECK) = UNEW(ICHECK,ICHECK) 
     1  * ( MOD (ICHECK, 100) /100.)
 3500   CONTINUE

370   CONTINUE
C        TEST FOR END OF RUN
      IF(NCYCLE .GE. ITMAX) THEN
      STOP
      ENDIF
C
C     TIME SMOOTHING AND UPDATE FOR NEXT CYCLE
C
c     IF(NCYCLE .LE. 1) THEN
c        CALL CALC3Z
c     ELSE
c        CALL CALC3
c     ENDIF
C Got rid of this control conditional (for now); just look at
the iterations after the first one, which is handled differently

c     SUBROUTINE CALC3 Inlined
C         TIME SMOOTHER
      DO 300 J=1,N
      DO 300 I=1,M
      UOLD(I,J) = U(I,J)+ALPHA*(UNEW(I,J)-2.*U(I,J)+UOLD(I,J))
      VOLD(I,J) = V(I,J)+ALPHA*(VNEW(I,J)-2.*V(I,J)+VOLD(I,J))
      POLD(I,J) = P(I,J)+ALPHA*(PNEW(I,J)-2.*P(I,J)+POLD(I,J))
      U(I,J) = UNEW(I,J)
      V(I,J) = VNEW(I,J)
      P(I,J) = PNEW(I,J)
  300 CONTINUE

C     PERIODIC CONTINUATION
      DO 320 J=1,N
      UOLD(M+1,J) = UOLD(1,J)
      VOLD(M+1,J) = VOLD(1,J)
      POLD(M+1,J) = POLD(1,J)
      U(M+1,J) = U(1,J)
      V(M+1,J) = V(1,J)
      P(M+1,J) = P(1,J)
  320 CONTINUE
      DO 325 I=1,M
      UOLD(I,N+1) = UOLD(I,1)
      VOLD(I,N+1) = VOLD(I,1)
      POLD(I,N+1) = POLD(I,1)
      U(I,N+1) = U(I,1)
      V(I,N+1) = V(I,1)
      P(I,N+1) = P(I,1)
  325 CONTINUE
      UOLD(M+1,N+1) = UOLD(1,1)
      VOLD(M+1,N+1) = VOLD(1,1)
      POLD(M+1,N+1) = POLD(1,1)
      U(M+1,N+1) = U(1,1)
      V(M+1,N+1) = V(1,1)
      P(M+1,N+1) = P(1,1)
C

      GO TO 90
      END
