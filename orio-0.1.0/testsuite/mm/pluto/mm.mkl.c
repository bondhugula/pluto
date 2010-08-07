
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mkl.h"

#define alpha 1
#define beta 1

#define M NCONT
#define N NCONT
#define K CONT
double A[M][K];
double B[K][N];
double C[M][N];

void init_arrays()
{
  int i, j;
  for (i=0; i<M; i++)
   for (j=0; j<K; j++)
    A[i][j] = (i+j) % 5 + 1;
  for (i=0; i<K; i++)
   for (j=0; j<N; j++)
    B[i][j] = (i+j) % 5 + 1;
  for (i=0; i<M; i++)
   for (j=0; j<N; j++)
    C[i][j] = 0;
}

void print_array()
{
  int i, j;
  for (i=0; i<M; i++) {
    for (j=0; j<N; j++) {
      fprintf(stderr, "%lf ", C[i][j]);
      if (j%80 == 79) fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
}

double rtclock()
{
  struct timezone tzp;
  struct timeval tp;
  int stat;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main()
{
  int i, j, k;
  int LDA, LDB, LDC;
  double *a, *b, *c;
  LDA=K;
  LDB=N;
  LDC=N;

  init_arrays();

  a = (double*) malloc(sizeof(double)* M*K);
  b = (double*) malloc(sizeof(double)* K*N);
  c = (double*) malloc(sizeof(double)* M*N);

  for(i=0; i<M; i++)
    for(k=0; k<K; k++)
      a[i*K+k]= A[i][k];
  for(k=0; k<K; k++)
    for(j=0; j<N; j++)
      b[k*N+j]= B[k][j];
  for(i=0; i<M; i++)
    for(j=0; j<N; j++)
      c[i*N+j]= C[i][j];

  double annot_t_start=0, annot_t_end=0, annot_t_total=0;
  int annot_i;

  for (annot_i=0; annot_i<REPS; annot_i++)
  {
    annot_t_start = rtclock();

    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,M,N,K,alpha,a,LDA,b,LDB,beta,c,LDC);
    
    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }

  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
    
  //for( i=0; i<M; i++)
  //for( j=0; j<N; j++)
  //C[i][j] = c[i*N+j];
  //print_array();

  return 1;
}
                                                                   
