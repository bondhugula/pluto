#include <stdio.h>
#include <stdlib.h>

#define T 100
#define N 100000
double a[T][N];

void init_array() {
    int i, j;
    for (i=0; i<T; i++)
        for (j=0; j<N; j++)
            a[i][j] = i+((double)j)/N;
}

int main() {
    init_array();

    /*@ profiled code @*/

    return a[0][0];   //needed to avoid dead code elimination
}
