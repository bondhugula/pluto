/*
 * John Conway's Game of Life - 2D, Periodic - B3/S23
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

/*
 * N is the length of one side of the GoL world
 * T is the number of generations to evolve
 */
#ifdef HAS_DECLS
#include "decls.h"
#else
#define N 16000L
#define T 500L
#endif

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

double t_start, t_end;


#define NUM_FP_OPS 8

/* Define our arrays */
int life[2][N][N];

double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}



/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y) {
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }

  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;

    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
   * tv_usec is certainly positive.
   */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

/*
 * Calculate the status of a cell in the next evolution given
 * the neighbors it has
 */
int b2s23(int cell, int neighbors) {
  if((cell == 1 && ((neighbors < 2) || (neighbors > 3)))) {
    return 0;
  }

  if((cell == 1 && (neighbors == 2 || neighbors == 3))) {
    return 1;
  }

  if((cell == 0 && neighbors == 3)) {
    return 1;
  }

  return cell;
}

/*
 * Print the final array
 */
void print_points(int t) {
  int a,b;

  for(a = 0; a < N; a++) {
    for(b = 0; b < N; b++) {
      fprintf(stderr, "%c ", (life[t%2][a][b] ? '1' : '.'));
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
}

int main(int argc, char * argv[]) {
  long int t, i, j, k;
  long int ip, im, jp, jm;
  int a, b, neighbors, flag;

  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  double total = 0.0;

#ifdef DEBUG
  printf("Size of the world = %ldx%ld\nNumber of generations to evolve = %ld\n", N, N, T);
#endif

  /* Initialization */
  srand(42); // seed with a constant value to verify results

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      life[0][i][j] = (rand() & 0x1) ? 1 : 0;
      life[1][i][j] = 0;
    }
  }

  /* // preset patterns */
  /* for (i = 0; i < N; i++) { */
  /*   for (j = 0; j < N; j++) { */
  /*     life[0][i][j] = 0; */
  /*     life[1][i][j] = 0; */
  /*   } */
  /* } */

  /* // toad - oscillator - period 2  */
  /* life[0][2][2] = 1; */
  /* life[0][2][3] = 1; */
  /* life[0][2][4] = 1; */
  /* life[0][3][1] = 1; */
  /* life[0][3][2] = 1; */
  /* life[0][3][3] = 1; */


  /* // pulsar - oscillator - period 3 */
  /* life[0][2][4] = 1; */
  /* life[0][2][5] = 1; */
  /* life[0][2][6] = 1; */
  /* life[0][2][10] = 1; */
  /* life[0][2][11] = 1; */
  /* life[0][2][12] = 1; */
  /* life[0][4][2] = 1; */
  /* life[0][4][7] = 1; */
  /* life[0][4][9] = 1; */
  /* life[0][4][14] = 1; */
  /* life[0][5][2] = 1; */
  /* life[0][5][7] = 1; */
  /* life[0][5][9] = 1; */
  /* life[0][5][14] = 1; */
  /* life[0][6][2] = 1; */
  /* life[0][6][7] = 1; */
  /* life[0][6][9] = 1; */
  /* life[0][6][14] = 1; */
  /* life[0][7][4] = 1; */
  /* life[0][7][5] = 1; */
  /* life[0][7][6] = 1; */
  /* life[0][7][10] = 1; */
  /* life[0][7][11] = 1; */
  /* life[0][7][12] = 1; */
  /* life[0][9][4] = 1; */
  /* life[0][9][5] = 1; */
  /* life[0][9][6] = 1; */
  /* life[0][9][10] = 1; */
  /* life[0][9][11] = 1; */
  /* life[0][9][12] = 1; */
  /* life[0][10][2] = 1; */
  /* life[0][10][7] = 1; */
  /* life[0][10][9] = 1; */
  /* life[0][10][14] = 1; */
  /* life[0][11][2] = 1; */
  /* life[0][11][7] = 1; */
  /* life[0][11][9] = 1; */
  /* life[0][11][14] = 1; */
  /* life[0][12][2] = 1; */
  /* life[0][12][7] = 1; */
  /* life[0][12][9] = 1; */
  /* life[0][12][14] = 1; */
  /* life[0][14][4] = 1; */
  /* life[0][14][5] = 1; */
  /* life[0][14][6] = 1; */
  /* life[0][14][10] = 1; */
  /* life[0][14][11] = 1; */
  /* life[0][14][12] = 1; */

#ifdef DEBUG
  // display the initial setting
  print_points(0);
#endif

    IF_TIME(t_start = rtclock());

#pragma scop
  for (t = 0; t < T; t++) {
    for (i = 1; i < N-1; i++) {
      for (j = 1; j < N-1; j++) {
        life[(t+1)%2][i][j] = b2s23(life[t%2][i][j], life[t%2][i-1][j+1] + life[t%2][i-1][j] + life[t%2][i-1][j-1] 
                    + life[t%2][i][j+1] + life[t%2][i][j-1]                     
                    + life[t%2][i+1][j+1] + life[t%2][i+1][j] + life[t%2][i+1][j-1]);
      }
    }
    /* print_points(1-k); */
  }
#pragma endscop

    IF_TIME(t_end = rtclock());
    IF_TIME(fprintf(stdout, "%0.6lfs\n", t_end - t_start));


  if (fopen(".test", "r")) {
      print_points(T);
  }

  return 0;
}

// icc -O3 -DTIME -DDEBUG life.c -o op-life
