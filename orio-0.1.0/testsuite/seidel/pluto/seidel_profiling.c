#include <stdio.h>
#include <stdlib.h>

#define N 600
#define T 300
double A[N][N];

void init_array() {
    int i, j;

    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = i*i+j*j;
        }
    }
}

int main() {
    init_array();

    /*@ profiled code @*/

    return A[0][0];   //needed to avoid dead code elimination
}
