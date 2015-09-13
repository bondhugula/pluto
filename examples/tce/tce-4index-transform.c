#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <assert.h>

#define N 55
double A[N][N][N][N];
double T1[N][N][N][N];
double T2[N][N][N][N];
double T3[N][N][N][N];
double B[N][N][N][N];
double C1[N][N];
double C2[N][N];
double C3[N][N];
double C4[N][N];
void init_array()
{
    int i,j,k,l;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            for (k=0; k<N; k++) {
                for (l=0; l<N; l++) {
                    A[i][j][k][l] = 1 + ((double)k)/N;
                }
            }
        }
    }
}

void init()
{
    init_array();
}

void print_array()
{
    int i, j, k, l;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            for (k=0; k<N; k++) {
                for (l=0; l<N; l++) {
                    fprintf(stdout, "%lf ", B[i][j][k][l]);
                }
            }
        }
    }

    fprintf(stdout, "\n");
}


void print_results()
{
    print_array();

}


int main()
{
    int a, q, r, s, p, b, c, d;

    init_array();

#pragma scop
    for(a=0; a<N; a++)
        for(q=0; q<N; q++)
            for(r=0; r<N; r++)
                for(s=0; s<N; s++)
                    for(p=0; p<N; p++)
                        T1[a][q][r][s] = T1[a][q][r][s] + A[p][q][r][s]*C4[p][a];

    for(a=0; a<N; a++)
        for(b=0; b<N; b++)
            for(r=0; r<N; r++)
                for(s=0; s<N; s++)
                    for(q=0; q<N; q++)
                        T2[a][b][r][s] = T2[a][b][r][s] + T1[a][q][r][s]*C3[q][b];

    for(a=0; a<N; a++)
        for(b=0; b<N; b++)
            for(c=0; c<N; c++)
                for(s=0; s<N; s++)
                    for(r=0; r<N; r++)
                        T3[a][b][c][s] = T3[a][b][c][s] + T2[a][b][r][s]*C2[r][c];

    for(a=0; a<N; a++)
        for(b=0; b<N; b++)
            for(c=0; c<N; c++)
                for(d=0; d<N; d++)
                    for(s=0; s<N; s++)
                        B[a][b][c][d] = B[a][b][c][d] + T3[a][b][c][s]*C1[s][d];
#pragma endscop

#ifdef TEST
    print_array();
#endif
    return 0;
}
