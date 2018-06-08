

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

typedef struct {
    int testid;
    char coord[1024];
    double tm;
} TimingInfo;



#ifndef REPS
#define REPS 1000
#endif



#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define SIZE 5000000
#define n 5000000
#define a1 4.232
#define a2 134531.2145
#define a3 43.24141
#define a4 241.24314
double x1[n];
double x2[n];
double x3[n];
double x4[n];
double y[n];

void init_arrays() {
    int i1;
    for (i1=0; i1<n; i1++)
        x1[i1] = (i1) % 5 + 1;
    for (i1=0; i1<n; i1++)
        x2[i1] = (i1) % 5 + 1;
    for (i1=0; i1<n; i1++)
        x3[i1] = (i1) % 5 + 1;
    for (i1=0; i1<n; i1++)
        x4[i1] = (i1) % 5 + 1;
    for (i1=0; i1<n; i1++)
        y[i1] = 0;
}

int main(int argc, char *argv[]) {
    int numprocs, myid, _i;
    TimingInfo mytimeinfo;
    TimingInfo *timevec;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

    /* Construct the MPI type for the timing info (what a pain!) */
    MPI_Datatype TimingInfoMPIType;
    MPI_Datatype type[3] = {MPI_INT, MPI_CHAR, MPI_DOUBLE};
    int blocklen[3] = {1,1024,1};
    MPI_Aint disp[3];
    int base;
    MPI_Address( &mytimeinfo.testid, disp);
    MPI_Address( &mytimeinfo.coord, disp+1);
    MPI_Address( &mytimeinfo.tm, disp+2);
    base = disp[0];
    for (_i=0; _i <3; _i++) disp[_i] -= base;
    MPI_Type_struct( 3, blocklen, disp, type, &TimingInfoMPIType);
    MPI_Type_commit( &TimingInfoMPIType);
    /* end of MPI type construction */

    if (myid == 0) timevec = (TimingInfo*) malloc(numprocs * sizeof(TimingInfo));
    init_arrays();

    double annot_t_start=0, annot_t_end=0, annot_t_total=0;
    int annot_i;

    if (myid == 0) {
        strcpy(mytimeinfo.coord,"[13]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-14; i=i+14) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-14; i=i+14) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 1) {
        strcpy(mytimeinfo.coord,"[15]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-16; i=i+16) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-16; i=i+16) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 2) {
        strcpy(mytimeinfo.coord,"[8]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-9; i=i+9) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-9; i=i+9) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 3) {
        strcpy(mytimeinfo.coord,"[9]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-10; i=i+10) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-10; i=i+10) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 4) {
        strcpy(mytimeinfo.coord,"[18]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-19; i=i+19) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                        y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
                        y[(i+17)]=y[(i+17)]+a1*x1[(i+17)]+a2*x2[(i+17)]+a3*x3[(i+17)]+a4*x4[(i+17)];
                        y[(i+18)]=y[(i+18)]+a1*x1[(i+18)]+a2*x2[(i+18)]+a3*x3[(i+18)]+a4*x4[(i+18)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-19; i=i+19) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                        y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
                        y[(i+17)]=y[(i+17)]+a1*x1[(i+17)]+a2*x2[(i+17)]+a3*x3[(i+17)]+a4*x4[(i+17)];
                        y[(i+18)]=y[(i+18)]+a1*x1[(i+18)]+a2*x2[(i+18)]+a3*x3[(i+18)]+a4*x4[(i+18)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 5) {
        strcpy(mytimeinfo.coord,"[10]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-11; i=i+11) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-11; i=i+11) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 6) {
        strcpy(mytimeinfo.coord,"[11]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-12; i=i+12) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-12; i=i+12) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 7) {
        strcpy(mytimeinfo.coord,"[6]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-7; i=i+7) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-7; i=i+7) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 8) {
        strcpy(mytimeinfo.coord,"[7]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-8; i=i+8) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-8; i=i+8) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 9) {
        strcpy(mytimeinfo.coord,"[12]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-13; i=i+13) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-13; i=i+13) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 10) {
        strcpy(mytimeinfo.coord,"[4]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-5; i=i+5) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-5; i=i+5) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 11) {
        strcpy(mytimeinfo.coord,"[5]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-6; i=i+6) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-6; i=i+6) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 12) {
        strcpy(mytimeinfo.coord,"[2]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-3; i=i+3) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-3; i=i+3) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 13) {
        strcpy(mytimeinfo.coord,"[3]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-4; i=i+4) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-4; i=i+4) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 14) {
        strcpy(mytimeinfo.coord,"[16]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-17; i=i+17) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                        y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-17; i=i+17) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                        y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 15) {
        strcpy(mytimeinfo.coord,"[17]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-18; i=i+18) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                        y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
                        y[(i+17)]=y[(i+17)]+a1*x1[(i+17)]+a2*x2[(i+17)]+a3*x3[(i+17)]+a4*x4[(i+17)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-18; i=i+18) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                        y[(i+15)]=y[(i+15)]+a1*x1[(i+15)]+a2*x2[(i+15)]+a3*x3[(i+15)]+a4*x4[(i+15)];
                        y[(i+16)]=y[(i+16)]+a1*x1[(i+16)]+a2*x2[(i+16)]+a3*x3[(i+16)]+a4*x4[(i+16)];
                        y[(i+17)]=y[(i+17)]+a1*x1[(i+17)]+a2*x2[(i+17)]+a3*x3[(i+17)]+a4*x4[(i+17)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 16) {
        strcpy(mytimeinfo.coord,"[1]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-2; i=i+2) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-2; i=i+2) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    } else if (myid == 17) {
        strcpy(mytimeinfo.coord,"[14]");
        mytimeinfo.testid = myid;
        for (annot_i=0; annot_i<REPS; annot_i++) {
            annot_t_start = MPI_Wtime();



            int i;

            /*@ begin Align (x1[],x2[],x3[],x4[],y[]) @*/
#pragma disjoint (*x1,*x2,*x3,*x4,*y)
            if ((((int)(x1)|(int)(x2)|(int)(x3)|(int)(x4)|(int)(y)) & 0xF) == 0) {
                __alignx(16,x1);
                __alignx(16,x2);
                __alignx(16,x3);
                __alignx(16,x4);
                __alignx(16,y);

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-15; i=i+15) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            } else {

                /*@ begin Loop (
                    transform Unroll(ufactor=UF)
                      for (i = 0; i <= n-1; i++)
                        y[i] = y[i] + a1*x1[i] + a2*x2[i] + a3*x3[i] + a4*x4[i];
                ) @*/  {
                    for (i=0; i<=n-15; i=i+15) {
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                        y[(i+1)]=y[(i+1)]+a1*x1[(i+1)]+a2*x2[(i+1)]+a3*x3[(i+1)]+a4*x4[(i+1)];
                        y[(i+2)]=y[(i+2)]+a1*x1[(i+2)]+a2*x2[(i+2)]+a3*x3[(i+2)]+a4*x4[(i+2)];
                        y[(i+3)]=y[(i+3)]+a1*x1[(i+3)]+a2*x2[(i+3)]+a3*x3[(i+3)]+a4*x4[(i+3)];
                        y[(i+4)]=y[(i+4)]+a1*x1[(i+4)]+a2*x2[(i+4)]+a3*x3[(i+4)]+a4*x4[(i+4)];
                        y[(i+5)]=y[(i+5)]+a1*x1[(i+5)]+a2*x2[(i+5)]+a3*x3[(i+5)]+a4*x4[(i+5)];
                        y[(i+6)]=y[(i+6)]+a1*x1[(i+6)]+a2*x2[(i+6)]+a3*x3[(i+6)]+a4*x4[(i+6)];
                        y[(i+7)]=y[(i+7)]+a1*x1[(i+7)]+a2*x2[(i+7)]+a3*x3[(i+7)]+a4*x4[(i+7)];
                        y[(i+8)]=y[(i+8)]+a1*x1[(i+8)]+a2*x2[(i+8)]+a3*x3[(i+8)]+a4*x4[(i+8)];
                        y[(i+9)]=y[(i+9)]+a1*x1[(i+9)]+a2*x2[(i+9)]+a3*x3[(i+9)]+a4*x4[(i+9)];
                        y[(i+10)]=y[(i+10)]+a1*x1[(i+10)]+a2*x2[(i+10)]+a3*x3[(i+10)]+a4*x4[(i+10)];
                        y[(i+11)]=y[(i+11)]+a1*x1[(i+11)]+a2*x2[(i+11)]+a3*x3[(i+11)]+a4*x4[(i+11)];
                        y[(i+12)]=y[(i+12)]+a1*x1[(i+12)]+a2*x2[(i+12)]+a3*x3[(i+12)]+a4*x4[(i+12)];
                        y[(i+13)]=y[(i+13)]+a1*x1[(i+13)]+a2*x2[(i+13)]+a3*x3[(i+13)]+a4*x4[(i+13)];
                        y[(i+14)]=y[(i+14)]+a1*x1[(i+14)]+a2*x2[(i+14)]+a3*x3[(i+14)]+a4*x4[(i+14)];
                    }
                    for (; i<=n-1; i=i+1)
                        y[i]=y[i]+a1*x1[i]+a2*x2[i]+a3*x3[i]+a4*x4[i];
                }
                /*@ end @*/

            }
            /*@ end @*/

        }
        annot_t_end = MPI_Wtime();
        annot_t_total += annot_t_end - annot_t_start;
        mytimeinfo.tm = annot_t_total / REPS;


    }


    MPI_Gather(&mytimeinfo, 1, TimingInfoMPIType, timevec, 1, TimingInfoMPIType, 0, MPI_COMM_WORLD);

    if (myid == 0) {
        printf("{'%s' : %g", mytimeinfo.coord, mytimeinfo.tm);
        for (_i = 1; _i < numprocs; _i++) {
            printf(", '%s' : %g", timevec[_i].coord, timevec[_i].tm);
        }
        printf("}\n");
    }

    MPI_Finalize();
    return 0;
}

