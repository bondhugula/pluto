#define N 10

//int N;

int main(void) {
    int i,j;
    int **a,**b;


#pragma scop
    for(i=0; i<3*N-1; i++)
        for(j=0; j<N; j++) {
            if((i+j>=N-1) && (i+j<=3*N-2))
                a[i][j] = 0;
            if((i+j>=2*N-1) && (i+j<=4*N-2))
                b[i][j] = a[i-N][j];
        }
#pragma endscop
    /*
    for(i=0; i<2; i++)
        {
    if((i>=0) && (i<=1))
      a[i][0] = 0;
    if((i>=1))
      b[i][0] = a[i-1][0];
        }
    */
    return 0;
}
