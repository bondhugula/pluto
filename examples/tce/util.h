
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

