/* pluto start (N1,N2,N3,M,N) */

DO T = 1, N3

    FSDX = 4/DX
    FSDY = 4/DY
    DO J = 1, N
        DO I=1,M
            CU(I+1,J) = 0.5*(P(I+1,J)+P(I,J))*U(I+1,J)
            CV(I,J+1) = 0.5*(P(I,J+1)+P(I,J))*V(I,J+1)
            Z(I+1,J+1) = (FSDX*(V(I+1,J+1)-V(I,J+1))-FSDY*(U(I+1,J+1)
                      -U(I+1,J)))/(P(I,J)+P(I+1,J)+P(I+1,J+1)+P(I,J+1))
            H(I,J) = P(I,J)+0.25*(U(I+1,J)*U(I+1,J)+U(I,J)*U(I,J)
               +V(I,J+1)*V(I,J+1)+V(I,J)*V(I,J))
        END DO
    END DO


    DO J=1,N
        CU(1,J) = CU(M+1,J)
        CV(M+1,J+1) = CV(1,J+1)
        Z(1,J+1) = Z(M+1,J+1)
        H(M+1,J) = H(1,J)
    END DO

    DO I=1,M
        CU(I+1,N+1) = CU(I+1,1)
        CV(I,1) = CV(I,N+1)
        Z(I+1,1) = Z(I+1,N+1)
        H(I,N+1) = H(I,1)
    END DO


    CU(1,N+1) = CU(M+1,1)
    CV(M+1,1) = CV(1,N+1)
    Z(1,1) = Z(M+1,N+1)
    H(M+1,N+1) = H(1,1)

    TDTS8 = TDT/8
    TDTSDX = TDT/DX
    TDTSDY = TDT/DY

    DO J=1,N
        DO I=1,M
            UNEW(I+1,J) = UOLD(I+1,J)+
            TDTS8*(Z(I+1,J+1)+Z(I+1,J))*(CV(I+1,J+1)+CV(I,J+1)+CV(I,J)
               +CV(I+1,J))-TDTSDX*(H(I+1,J)-H(I,J))
            VNEW(I,J+1) = VOLD(I,J+1)-TDTS8*(Z(I+1,J+1)+Z(I,J+1))
               *(CU(I+1,J+1)+CU(I,J+1)+CU(I,J)+CU(I+1,J))
                  -TDTSDY*(H(I,J+1)-H(I,J))
            PNEW(I,J) = POLD(I,J)-TDTSDX*(CU(I+1,J)-CU(I,J))
                   -TDTSDY*(CV(I,J+1)-CV(I,J))
        END DO
    END DO

    DO J=1,N
        UNEW(1,J) = UNEW(M+1,J)
        VNEW(M+1,J+1) = VNEW(1,J+1)
        PNEW(M+1,J) = PNEW(1,J)
    END DO 

    DO I=1,M
        UNEW(I+1,N+1) = UNEW(I+1,1)
        VNEW(I,1) = VNEW(I,N+1)
        PNEW(I,N+1) = PNEW(I,1)
    END DO 

    UNEW(1,N+1) = UNEW(M+1,1)
    VNEW(M+1,1) = VNEW(1,N+1)
    PNEW(M+1,N+1) = PNEW(1,1)
    TIME = TIME + DT


    DO J=1,N
        DO I=1,M
            UOLD(I,J) = U(I,J)+ALPHA*(UNEW(I,J)-2*U(I,J)+UOLD(I,J))
            VOLD(I,J) = V(I,J)+ALPHA*(VNEW(I,J)-2*V(I,J)+VOLD(I,J))
            POLD(I,J) = P(I,J)+ALPHA*(PNEW(I,J)-2*P(I,J)+POLD(I,J))
            U(I,J) = UNEW(I,J)
            V(I,J) = VNEW(I,J)
            P(I,J) = PNEW(I,J)
        END DO 
    END DO

    DO J=1,N
        UOLD(M+1,J) = UOLD(1,J)
        VOLD(M+1,J) = VOLD(1,J)
        POLD(M+1,J) = POLD(1,J)
        U(M+1,J) = U(1,J)
        V(M+1,J) = V(1,J)
        P(M+1,J) = P(1,J)
    END DO

    DO I=1,M
        UOLD(I,N+1) = UOLD(I,1)
        VOLD(I,N+1) = VOLD(I,1)
        POLD(I,N+1) = POLD(I,1)
        U(I,N+1) = U(I,1)
        V(I,N+1) = V(I,1)
        P(I,N+1) = P(I,1)
    END DO

    UOLD(M+1,N+1) = UOLD(1,1)
    VOLD(M+1,N+1) = VOLD(1,1)
    POLD(M+1,N+1) = POLD(1,1)
    U(M+1,N+1) = U(1,1)
    V(M+1,N+1) = V(1,1)
    P(M+1,N+1) = P(1,1)

END DO
/* pluto end */
