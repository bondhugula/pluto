/**
 * gesummv.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 4000. */
#include "gesummv.h"

  DATA_TYPE alpha;
  DATA_TYPE beta;

/* Array initialization. */
static
void init_array(int n,
		DATA_TYPE *alpha,
		DATA_TYPE *beta,
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		DATA_TYPE POLYBENCH_1D(x,N,n),
        DATA_TYPE POLYBENCH_1D(tmp,N,n),
        DATA_TYPE POLYBENCH_1D(y,N,n)
        )
{
  int i, j;

  *alpha = 43532;
  *beta = 12313;
  for (i = 0; i < n; i++)
    {
      x[i] = ((DATA_TYPE) i) / n;
      tmp[i] = 0;
      y[i] = 0;
      for (j = 0; j < n; j++) {
	A[i][j] = ((DATA_TYPE) i*j) / n;
	B[i][j] = ((DATA_TYPE) i*j) / n;
      }
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(y,N,n))

{
  int i;

  for (i = 0; i < n; i++) {
    fprintf (stderr, DATA_PRINTF_MODIFIER, y[i]);
    if (i % 20 == 0) fprintf (stderr, "\n");
  }
}

#define __DECLARATION_OF_alpha DATA_TYPE alpha
#define __DECLARATION_OF_beta DATA_TYPE beta
#define __DECLARATION_OF_A DATA_TYPE POLYBENCH_2D(A,N,N,n,n)
#define __DECLARATION_OF_B DATA_TYPE POLYBENCH_2D(B,N,N,n,n)
#define __DECLARATION_OF_tmp DATA_TYPE POLYBENCH_1D(tmp,N,n)
#define __DECLARATION_OF_x DATA_TYPE POLYBENCH_1D(x,N,n)
#define __DECLARATION_OF_y DATA_TYPE POLYBENCH_1D(y,N,n)

#define __PLACE_TO_INSERT_FORWARD_DECLARATIONS

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_gesummv(int n,
		    DATA_TYPE alpha,
		    DATA_TYPE beta
#ifndef USE_LOCAL_ARRAYS
		    ,DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		    DATA_TYPE POLYBENCH_2D(B,N,N,n,n),
		    DATA_TYPE POLYBENCH_1D(tmp,N,n),
		    DATA_TYPE POLYBENCH_1D(x,N,n),
		    DATA_TYPE POLYBENCH_1D(y,N,n)
#endif
    )
{
  int i, j;

#pragma scop
  for (i = 0; i < _PB_N; i++)
    {
      for (j = 0; j < _PB_N; j++)
	{
	  tmp[i] = A[i][j] * x[j] + tmp[i];
	  y[i] = B[i][j] * x[j] + y[i];
	}
      y[i] = alpha * tmp[i] + beta * y[i];
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

#ifndef USE_LOCAL_ARRAYS
  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n, &alpha, &beta,
	      POLYBENCH_ARRAY(A),
	      POLYBENCH_ARRAY(B),
	      POLYBENCH_ARRAY(x),
	      POLYBENCH_ARRAY(tmp),
	      POLYBENCH_ARRAY(y));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_gesummv (n, alpha, beta,
		  POLYBENCH_ARRAY(A),
		  POLYBENCH_ARRAY(B),
		  POLYBENCH_ARRAY(tmp),
		  POLYBENCH_ARRAY(x),
		  POLYBENCH_ARRAY(y));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
#else

  /* Run kernel. */
  kernel_gesummv (n, alpha, beta);

#endif
  return 0;
}
