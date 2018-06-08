#include <stdio.h>
#include <stdlib.h>

#define tmax 400
#define nx 800
#define ny 800
double ex[nx][ny+1];
double ey[nx+1][ny];
double hz[nx][ny];

void init_array() {
    int i, j;

    for (i=0; i<nx+1; i++)  {
        for (j=0; j<ny; j++)  {
            ey[i][j] = 0;
        }
    }

    for (i=0; i<nx; i++)  {
        for (j=0; j<ny+1; j++)  {
            ex[i][j] = 0;
        }
    }

    for (j=0; j<ny; j++)  {
        ey[0][j] = ((double)j)/ny;
    }

    for (i=0; i<nx; i++)    {
        for (j=0; j<ny; j++)  {
            hz[i][j] = 0;
        }
    }
}

int main() {
    init_array();

    /*@ profiled code @*/

    return hz[0][0];   //needed to avoid dead code elimination
}
