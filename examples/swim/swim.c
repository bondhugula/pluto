#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#include "decls.h"
#include "util.h"

#ifdef PERFCTR
#include <papi.h>
#include "papi_defs.h"
#endif

double t_start, t_end;

main ()
{
    int t, i, j;
    int TDT;

    int M, N, MP1, NP1;
    double A, EL, DI, DJ, TPI, PCF, PI;
    double DX, DY, DT, ALPHA;

    DT = 20;


    M = N1 -1;
    N = N2 -1;
    DX = 0.25E5;
    DY = 0.25E5;
    A = 1.0E6;
    ALPHA = 0.001;

    TDT = DT;
    MP1 = M+1;
    NP1 = N+1;
    EL = N*DX;
    PI = 4.0*atan(1.0);
    TPI = PI+PI;

    DI = TPI/M;
    DJ = TPI/N;
    PCF = PI*PI*A*A/(EL*EL);

    /* Initial VALUES OF THE STREAM FUNCTION AND P */
    for (i=0; i<MP1; i++)   {
        for (j=0; j<NP1; j++)   {
            PSI[i][j] = A*sin((i-0.5)*DI)*sin((j-0.5)*DJ);
            P[i][j] = PCF*(cos(2.0*(i-1)*DI)
                    +cos(2.0*(j-1)*DJ))+50000.0;
        }
    }

    /* Initialize Velocities */
    for (i=0; i<M; i++) {
        for (j=0; j<N; j++)   {
            U[i+1][j] = -(PSI[i+1][j+1]-PSI[i+1][j])/DY;
            V[i][j+1] = (PSI[i+1][j+1]-PSI[i][j+1])/DX;
        }
    }

    /* Periodic Continuation */
    for (j=0; j<N; j++){ 
        U[0][j] = U[M][j];
        V[M][j+1] = V[0][j+1];
    }
    for (i=0; i<M; i++) {
        U[i+1][N] = U[i+1][0];
        V[i][0] = V[i][N];
    }

    U[0][N] = U[M][0];
    V[M][0] = V[0][N];
    for (i=0; i<MP1; i++)   {
        for (j=0; j<NP1; j++)   {
            UOLD[i][j] = U[i][j];
            VOLD[i][j] = V[i][j];
            POLD[i][j] = P[i][j];
        }
    }
    /* END OF initialization */

#ifdef PERFCTR
    PERF_INIT; 
#endif

    IF_TIME(t_start = rtclock());

    TDT=DT;

    /* pluto start (M,N,N3) */

    for (t=0; t<N3; t++) {

        FSDX = 4.0/DX;
        FSDY = 4.0/DY;

        for (i=0; i<M; i++) {
            for (j=0; j<N; j++) {
                CU[i+1][j] = 0.5*(P[i+1][j]+P[i][j])*U[i+1][j];
                CV[i][j+1] = 0.5*(P[i][j+1]+P[i][j])*V[i][j+1];
                Z[i+1][j+1] = (FSDX*(V[i+1][j+1]-V[i][j+1])-FSDY*(U[i+1][j+1]
                            -U[i+1][j]))/(P[i][j]+P[i+1][j]+P[i+1][j+1]+P[i][j+1] + 0.5);
                H[i][j] = P[i][j]+0.25*(U[i+1][j]*U[i+1][j]+U[i][j]*U[i][j]
                        +V[i][j+1]*V[i][j+1]+V[i][j]*V[i][j]);
            }
        }


        for (j=0; j<N; j++)  {
            CU[0][j] = CU[M+1][j];
            CV[M][j+1] = CV[0][j+1];
            Z[0][j+1] = Z[M][j+1];
            H[M][j] = H[0][j];
        }

        for (i=0; i<M; i++) {
            CU[i+1][N] = CU[i+1][0];
            CV[i][0] = CV[i][N];
            Z[i+1][0] = Z[i+1][N];
            H[i][N] = H[i][0];
        }


        CU[0][N] = CU[M][0];
        CV[M][0] = CV[0][N];
        Z[0][0] = Z[M][N];
        H[M][N] = H[0][0];

        TDTS8 = TDT/8.0;
        TDTSDX = TDT/DX;
        TDTSDY = TDT/DY;


        for (i=0; i<M; i++) {
            for (j=0; j<N; j++) {
                UNEW[i+1][j] = UOLD[i+1][j]+
                    TDTS8*(Z[i+1][j+1]+Z[i+1][j])*(CV[i+1][j+1]+CV[i][j+1]+CV[i][j]
                            +CV[i+1][j])-TDTSDX*(H[i+1][j]-H[i][j]);
                VNEW[i][j+1] = VOLD[i][j+1]-TDTS8*(Z[i+1][j+1]+Z[i][j+1])
                    *(CU[i+1][j+1]+CU[i][j+1]+CU[i][j]+CU[i+1][j])
                    -TDTSDY*(H[i][j+1]-H[i][j]);
                PNEW[i][j] = POLD[i][j]-TDTSDX*(CU[i+1][j]-CU[i][j])
                    -TDTSDY*(CV[i][j+1]-CV[i][j]);
            }
        }

        for (j=0; j<N; j++) {
            UNEW[0][j] = UNEW[M][j];
            VNEW[M][j+1] = VNEW[0][j+1];
            PNEW[M][j] = PNEW[0][j];
        } 

        for (i=0; i<M; i++) {
            UNEW[i+1][N] = UNEW[i+1][0];
            VNEW[i][0] = VNEW[i][N];
            PNEW[i][N] = PNEW[i][0];
        } 

        UNEW[0][N] = UNEW[M][0];
        VNEW[M][0] = VNEW[0][N];
        PNEW[M][N] = PNEW[0][0];

        for (i=0; i<M; i++) {
            for (j=0; j<N; j++) {
                UOLD[i][j] = U[i][j]+ALPHA*(UNEW[i][j]-2*U[i][j]+UOLD[i][j]);
                VOLD[i][j] = V[i][j]+ALPHA*(VNEW[i][j]-2*V[i][j]+VOLD[i][j]);
                POLD[i][j] = P[i][j]+ALPHA*(PNEW[i][j]-2*P[i][j]+POLD[i][j]);
                U[i][j] = UNEW[i][j];
                V[i][j] = VNEW[i][j];
                P[i][j] = PNEW[i][j];
            } 
        }

        for (j=0; j<N; j++) {
            UOLD[M][j] = UOLD[0][j];
            VOLD[M][j] = VOLD[0][j];
            POLD[M][j] = POLD[0][j];
            U[M][j] = U[0][j];
            V[M][j] = V[0][j];
            P[M][j] = P[0][j];
        }

        for (i=0; i<M; i++) {
            UOLD[i][N] = UOLD[i][0];
            VOLD[i][N] = VOLD[i][0];
            POLD[i][N] = POLD[i][0];
            U[i][N] = U[i][0];
            V[i][N] = V[i][0];
            P[i][N] = P[i][0];
        }

        UOLD[M][N] = UOLD[0][0];
        VOLD[M][N] = VOLD[0][0];
        POLD[M][N] = POLD[0][0];
        U[M][N] = U[0][0];
        V[M][N] = V[0][0];
        P[M][N] = P[0][0];

    }
    /* pluto end */

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));

#ifdef TEST
    print();
#endif
}

