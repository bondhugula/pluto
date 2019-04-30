/*
 * Calculating the price of American Put Option
 * Adapted from Pochoir test bench
 *
 * Irshad Pananilath: irshad@csa.iisc.ernet.in
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define SIMPLE 0
#define N_RANK 1

/* apop_pochoir -S 100 -E 95 -r 10 -V 30 -T 1 -s 2000000 -t 10000 */

/* #define DEFAULT_S 100.00 */
/* #define DEFAULT_E 100.00 */
/* #define DEFAULT_r   0.10 */
/* #define DEFAULT_V   0.25 */
/* #define DEFAULT_T   1.00 */
/*  */
/* #define DEFAULT_s 100 */
/* #define DEFAULT_t 100 */

#define DEFAULT_S 100.00
#define DEFAULT_E 95.00
#define DEFAULT_r   0.10
#define DEFAULT_V   0.30
#define DEFAULT_T   1.00

#define DEFAULT_s 2000000
#define DEFAULT_t 10000

double C[3][DEFAULT_s + 1]; // TODO::change the space dimension
double F[2][DEFAULT_s + 1]; // TODO::change the space dimension

/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec)
  {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }

  if (x->tv_usec - y->tv_usec > 1000000)
  {
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

void print_usage( char *prog )
{
  printf( "Usage: %s [ options ]\n\n", prog );

  printf( "Options:\n" );

  printf( "\t-S value : spot price ( default: %0.2lf )\n", DEFAULT_S );
  printf( "\t-E value : exercise price ( default: %0.2lf )\n", DEFAULT_E );
  printf( "\t-r value : interest rate ( default: %0.2lf )\n", DEFAULT_r * 100 );
  printf( "\t-V value : volatility ( default: %0.2lf )\n", DEFAULT_V * 100 );
  printf( "\t-T value : time to mature in years ( default: %0.2lf )\n\n", DEFAULT_T );

  printf( "\t-s value : steps in space dimension ( default: %d )\n", DEFAULT_s );
  printf( "\t-t value : steps in time dimension ( default: %d )\n\n", DEFAULT_t );

  printf( "\t-i               : Run iterative stencil\n\n" );

  printf( "\t-h               : print this help screen\n\n" );
}

void computeCoeffs( double r, double V, double T, int ns, int nt) {
  double V2 = V * V;
  double dt = T / nt;
  double r1 = 1.0 / ( 1.0 + r * dt );
  double r2 = dt / ( 1.0 + r * dt );
  int x;

  for ( x = 0; x <= ns; ++x ) {
    C[0][x] = r2 * 0.5 * x * ( - r + V2 * x );
    C[1][x] = r1 * ( 1 - V2 * x * x * dt );
    C[2][x] = r2 * 0.5 * x * ( r + V2 * x );
  }
}

int main( int argc, char *argv[ ] ) {
  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  double S, E, r, V, T;
  int ns, nt;
  int i, t;
  double v, price1;

  S = DEFAULT_S;
  E = DEFAULT_E;
  r = DEFAULT_r;
  V = DEFAULT_V;
  T = DEFAULT_T;

  ns = DEFAULT_s;
  nt = DEFAULT_t;

#ifdef DEBUG
  printf( "\nStencil-based DP for the price of American put option ( Run with option -h for help ).\n\n" );
  printf( "Parameters:\n\n" );

  printf( "\t spot price = %0.2lf\n", S );
  printf( "\t exercise price = %0.2lf\n", E );
  printf( "\t interest rate = %0.2lf\%\n", r * 100 );
  printf( "\t volatility = %0.2lf\%\n", V * 100 );
  printf( "\t time to mature ( in years ) = %0.2lf\n\n", T );

  printf( "\t steps in space dimension = %d\n", ns );
  printf( "\t steps in time dimension = %d\n\n", nt );
#endif

  ns = ns + ( ns & 1 );
  double dS = 2.0 * S / ns;

  /* computeCoeffs( r, V, T, ns, nt ); */

  // initialize
  double V2 = V * V;
  double dt = T / nt;
  double r1 = 1.0 / ( 1.0 + r * dt );
  double r2 = dt / ( 1.0 + r * dt );
  int x;

  for ( x = 0; x <= ns; ++x ) {
    C[0][x] = r2 * 0.5 * x * ( - r + V2 * x );
    C[1][x] = r1 * ( 1 - V2 * x * x * dt );
    C[2][x] = r2 * 0.5 * x * ( r + V2 * x );
  }

  for ( i = 0; i <= ns; ++i )
    F[0][i] = max( 0.0, E - i * dS );

  F[1][0] = E;

  for (t = 0; t < nt; ++t ) {
    F[0][ns] = 0;
    F[1][ns] = 0;
  }

#ifdef TIME
  gettimeofday(&start, 0);
#endif

#pragma scop
  for ( t = 0; t < nt; ++t ) {
    for ( i = 1; i < ns; ++i ) {
      F[(t+1)%2][i] = max(C[0][i] * F[t%2][i-1] + C[1][i] * F[t%2][i] + C[2][i] * F[t%2][i+1], E - i * dS );
    }
  }
#pragma endscop

  price1 = F[nt%2][( ns >> 1 )];

#ifdef TIME
  gettimeofday(&end, 0);

  ts_return = timeval_subtract(&result, &end, &start);
  tdiff = (double) (result.tv_sec + result.tv_usec * 1.0e-6);

  printf( "\n\nIterative Stencil:\n" );
  printf( "\t option price = %.2lf\n", price1);

  printf("Time taken: %7.5lfms\n", tdiff * 1.0e3);
#endif

  return 0;
}
