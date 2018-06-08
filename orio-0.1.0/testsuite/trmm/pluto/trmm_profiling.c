#include <stdio.h>
#include <stdlib.h>

#define N 700
#define alpha 1
double A[N][N+20];
double B[N][N+20];

void init_array() {
    int i,j;
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            B[i][j] = (i+j) % 5 + 1;
    for (i=0; i<N; i++)
        for (j=0; j<N; j++)
            if (i < j)
                A[i][j] = (i+j) % 5 + 1;
            else if (i == j)
                A[i][j] = 1;
            else
                A[i][j] = -1;
}

int main() {
    init_array();

    /*@ profiled code @*/

    return B[0][0];   //needed to avoid dead code elimination
}
