/**
 * atax.c: This file is part of the PolyBench/C 3.2 test suite.
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
#include "atax.h"


/* Array initialization. */
static
void init_array (int nx, int ny,
		 DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(x,NY,ny),
		 DATA_TYPE POLYBENCH_1D(y,NY,ny),
		 DATA_TYPE POLYBENCH_1D(tmp,NX,nx))
{
  int i, j;

  for (i = 0; i < ny; i++)
      x[i] = i * M_PI;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      A[i][j] = ((DATA_TYPE) i*(j+1)) / nx;

  for (i = 0; i < ny; i++)
      y[i] = 0;

  for (i = 0; i < nx; i++)
      tmp[i] = 0;
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int nx,
		 DATA_TYPE POLYBENCH_1D(y,NX,nx))

{
  int i;

  for (i = 0; i < nx; i++) {
    fprintf (stderr, DATA_PRINTF_MODIFIER, y[i]);
    if (i % 20 == 0) fprintf (stderr, "\n");
  }
  fprintf (stderr, "\n");
}

#define __DECLARATION_OF_A DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny)
#define __DECLARATION_OF_x DATA_TYPE POLYBENCH_1D(x,NY,ny)
#define __DECLARATION_OF_y DATA_TYPE POLYBENCH_1D(y,NY,ny)
#define __DECLARATION_OF_tmp DATA_TYPE POLYBENCH_1D(tmp,NX,nx)

#define __PLACE_TO_INSERT_FORWARD_DECLARATIONS

/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_atax(int nx, int ny
#ifndef USE_LOCAL_ARRAYS
		 ,DATA_TYPE POLYBENCH_2D(A,NX,NY,nx,ny),
		 DATA_TYPE POLYBENCH_1D(x,NY,ny),
		 DATA_TYPE POLYBENCH_1D(y,NY,ny),
		 DATA_TYPE POLYBENCH_1D(tmp,NX,nx)
#endif
    )
{
  int i, j;

#pragma scop
  for (i = 0; i < _PB_NX; i++)
    {
      for (j = 0; j < _PB_NY; j++)
	tmp[i] = tmp[i] + A[i][j] * x[j];
      for (j = 0; j < _PB_NY; j++)
	y[j] = y[j] + A[i][j] * tmp[i];
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int nx = NX;
  int ny = NY;

#ifndef USE_LOCAL_ARRAYS
  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, NX, NY, nx, ny);
  POLYBENCH_1D_ARRAY_DECL(x, DATA_TYPE, NY, ny);
  POLYBENCH_1D_ARRAY_DECL(y, DATA_TYPE, NY, ny);
  POLYBENCH_1D_ARRAY_DECL(tmp, DATA_TYPE, NX, nx);

  /* Initialize array(s). */
  init_array (nx, ny, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(x),
	       POLYBENCH_ARRAY(y),
	       POLYBENCH_ARRAY(tmp));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_atax (nx, ny,
	       POLYBENCH_ARRAY(A),
	       POLYBENCH_ARRAY(x),
	       POLYBENCH_ARRAY(y),
	       POLYBENCH_ARRAY(tmp));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(nx, POLYBENCH_ARRAY(y)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(x);
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(tmp);

#else
      /* Run kernel. */
      kernel_atax (nx, ny);

#endif
  return 0;
}
