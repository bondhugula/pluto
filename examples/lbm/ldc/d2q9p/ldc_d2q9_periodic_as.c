/************************************************************************/
/* Lid Driven Cavity(LDC) - 2D simulation using LBM with a D2Q9 lattice */
/*                                                                      */
/* Scheme: 2-Grid, Pull, Compute, Keep                                  */
/*                                                                      */
/* Lattice naming convention                                            */
/*   c6[NW]  c2[N]  c5[NE]                                              */
/*   c3[W]   c0[C]  c1[E]                                               */
/*   c7[SW]  c4[S]  c8[SE]                                              */
/*                                                                      */
/* Author: Irshad Pananilath <pmirshad+code@gmail.com>                  */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

#ifndef nX
#define nX 1024  // Size of domain along x-axis
#endif

#ifndef nY
#define nY 1024  // Size of domain along y-axis
#endif

#ifndef nK
#define nK 9  // D2Q9
#endif

#ifndef nTimesteps
#define nTimesteps 50000  // Number of timesteps for the simulation
#endif

#define C 0
#define N 2
#define S 4
#define E 1
#define W 3
#define NE 5
#define NW 6
#define SE 8
#define SW 7

/* Tunable parameters */
const double re = 100;
const double uTop = 0.01;

/* Calculated parameters: Set in main */
double omega = 0.0;
double one_minus_omega = 0.0;

const double w1 = 4.0 / 9.0;
const double w2 = 1.0 / 9.0;
const double w3 = 1.0 / 36.0;

double grid[2][nY][nX][nK] __attribute__((aligned(32)));

void lbm_kernel(double, double, double, double, double, double, double, double,
                double, double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict, int, int,
                int);
void dumpVelocities(int timestep);
int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y);

int main(void) {
  int t, y, x, k;
  double total_lattice_pts = (double)nY * (double)nX * (double)nTimesteps;

  /* For timekeeping */
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  /* Compute values for global parameters */
  omega = 2.0 / ((6.0 * uTop * (nX - 0.5) / re) + 1.0);
  one_minus_omega = 1.0 - omega;

  printf(
      "2D Lid Driven Cavity simulation with D2Q9 lattice:\n"
      "\tfilename   : %s\n"
      "\tscheme     : 2-Grid, Fused, Pull\n"
      "\tgrid size  : %d x %d = %.2lf * 10^3 Cells\n"
      "\tnTimeSteps : %d\n"
      "\tRe         : %.2lf\n"
      "\tuTop       : %.6lf\n"
      "\tomega      : %.6lf\n",
      __FILE__, nX, nY, nX * nY / 1.0e3, nTimesteps, re, uTop, omega);

  /* Initialize all 9 PDFs for each point in the domain to 1.0 */
  for (y = 0; y < nY; y++) {
    for (x = 0; x < nX; x++) {
      for (k = 0; k < nK; k++) {
        grid[0][y][x][k] = 1.0;
        grid[1][y][x][k] = 1.0;
      }
    }
  }

  /* To satisfy PET */
  short _nX = nX/2;
  short _nY = nY/2;
  int _nTimesteps = nTimesteps;

#ifdef TIME
  gettimeofday(&start, 0);
#endif

#pragma scop
  /* Kernel of the code */
  for (t = 0; t < _nTimesteps; t++) {
#pragma omp parallel for private(x)
    for (y = 0; y < 2*_nY; y++) {
#pragma ivdep
      for (x = 0; x < 2*_nX; x++) {
        lbm_kernel(grid[t % 2][y][x][C],
                   grid[t % 2][y==0?2*_nY-1:y - 1][x + 0][N],
                   grid[t % 2][y==2*_nY-1?0:y + 1][x + 0][S],
                   grid[t % 2][y + 0][x==0?2*_nX-1:x - 1][E],
                   grid[t % 2][y + 0][x==2*_nX-1?0:x + 1][W],
                   grid[t % 2][y==0?2*_nY-1:y - 1][x==0?2*_nX-1:x - 1][NE],
                   grid[t % 2][y==0?2*_nY-1:y-1][x==2*_nX-1?0:x+1][NW],
                   grid[t % 2][y==2*_nY-1?0:y+1][x==0?2*_nX-1:x-1][SE],
                   grid[t % 2][y==2*_nY-1?0:y+1][x==2*_nX-1?0:x+1][SW],
                   &grid[(t + 1) % 2][y][x][C],
                   &grid[(t + 1) % 2][y][x][N],
                   &grid[(t + 1) % 2][y][x][S],
                   &grid[(t + 1) % 2][y][x][E],
                   &grid[(t + 1) % 2][y][x][W],
                   &grid[(t + 1) % 2][y][x][NE],
                   &grid[(t + 1) % 2][y][x][NW],
                   &grid[(t + 1) % 2][y][x][SE],
                   &grid[(t + 1) % 2][y][x][SW],
                   t, y, x);
      }
    }
  }
#pragma endscop

#ifdef TIME
  gettimeofday(&end, 0);

  ts_return = timeval_subtract(&result, &end, &start);
  tdiff = (double)(result.tv_sec + result.tv_usec * 1.0e-6);

  printf("\tTime taken : %7.5lfm\n", tdiff / 60.0);
  printf("\tMLUPS      : %7.5lf\n", (total_lattice_pts / (1.0e6 * tdiff)));
#endif

#ifdef DEBUG
  /* Dump rho, uX, uY for the entire domain to verify results */
  dumpVelocities(t);
#endif

  return 0;
}

/* Performs LBM kernel for one lattice cell */
void lbm_kernel(double in_c, double in_n, double in_s, double in_e, double in_w,
                double in_ne, double in_nw, double in_se, double in_sw,
                double *restrict out_c, double *restrict out_n,
                double *restrict out_s, double *restrict out_e,
                double *restrict out_w, double *restrict out_ne,
                double *restrict out_nw, double *restrict out_se,
                double *restrict out_sw, int t, int y, int x) {
  /* %%INIT%% */
  double lbm_sum_e, lbm_sum_w, lbm_sum_n, lbm_sum_s, lbm_rho, lbm_inv_rho;
  double lbm_uX, lbm_uY, lbm_uX2, lbm_uY2, lbm_u2, lbm_one_minus_u2;
  double lbm_omega_w1_rho, lbm_omega_w2_rho, lbm_omega_w3_rho;
  double lbm_30_uY, lbm_45_uY2, lbm_30_uX, lbm_45_uX2;

  double lbm_out_c, lbm_out_n, lbm_out_s, lbm_out_e, lbm_out_w, lbm_out_ne,
         lbm_out_nw, lbm_out_se, lbm_out_sw;
  double lbm_uX_plus_uY, lbm_uX_minus_uY, lbm_30_uX_plus_uY, lbm_30_uX_minus_uY,
         lbm_45_uX_plus_uY, lbm_45_uX_minus_uY;
  /* %%INIT%% */

  /* #<{(| %%BOUNDARY%% |)}># */
  /* if (((y == 0 || y == nY - 1)) || */
  /*     ((x == 0 || x == nX - 1) ) { */
  /*   #<{(| Bounce back boundary condition at walls |)}># */
  /*   *out_c = in_c; */
  /*  */
  /*   *out_s = in_n; */
  /*   *out_n = in_s; */
  /*   *out_w = in_e; */
  /*   *out_e = in_w; */
  /*  */
  /*   *out_sw = in_ne; */
  /*   *out_se = in_nw; */
  /*   *out_nw = in_se; */
  /*   *out_ne = in_sw; */
  /*  */
  /*   return; */
  /* } */
  /* %%BOUNDARY%% */

  /* Algebraic simplification */

  lbm_sum_e = in_ne + in_e + in_se;
  lbm_sum_w = in_nw + in_w + in_sw;
  lbm_sum_n = in_nw + in_n + in_ne;
  lbm_sum_s = in_sw + in_s + in_se;

  /* Compute rho */
  lbm_rho = lbm_sum_e + (in_n + in_c + in_s) + lbm_sum_w;
  lbm_inv_rho = 1.0 / lbm_rho;

  /* Compute velocity along x-axis */
  lbm_uX = lbm_inv_rho * (lbm_sum_e - lbm_sum_w);

  /* Compute velocity along y-axis */
  lbm_uY = lbm_inv_rho * (lbm_sum_n - lbm_sum_s);

  /* %%BOUNDARY%% */
  /* Impose the constantly moving upper wall */
  if (y == nY -1){
    lbm_uX = uTop;
    lbm_uY = 0.00;
  }
  /* %%BOUNDARY%% */

  lbm_uX2 = lbm_uX * lbm_uX;
  lbm_uY2 = lbm_uY * lbm_uY;

  lbm_u2 = 1.5 * (lbm_uX2 + lbm_uY2);

  lbm_one_minus_u2 = 1.0 - lbm_u2;

  lbm_omega_w1_rho = omega * w1 * lbm_rho;
  lbm_omega_w2_rho = omega * w2 * lbm_rho;
  lbm_omega_w3_rho = omega * w3 * lbm_rho;

  lbm_30_uY = 3.0 * lbm_uY;
  lbm_45_uY2 = 4.5 * lbm_uY2;
  lbm_30_uX = 3.0 * lbm_uX;
  lbm_45_uX2 = 4.5 * lbm_uX2;

  lbm_uX_plus_uY = lbm_uX + lbm_uY;
  lbm_uX_minus_uY = lbm_uX - lbm_uY;

  lbm_30_uX_plus_uY = 3.0 * lbm_uX_plus_uY;
  lbm_30_uX_minus_uY = 3.0 * lbm_uX_minus_uY;
  lbm_45_uX_plus_uY = 4.5 * lbm_uX_plus_uY * lbm_uX_plus_uY;
  lbm_45_uX_minus_uY = 4.5 * lbm_uX_minus_uY * lbm_uX_minus_uY;

  /* Compute and keep new PDFs */
  lbm_out_c = one_minus_omega * in_c + lbm_omega_w1_rho * (lbm_one_minus_u2);

  lbm_out_n = one_minus_omega * in_n + lbm_omega_w2_rho * (lbm_one_minus_u2 + lbm_30_uY + lbm_45_uY2);
  lbm_out_s = one_minus_omega * in_s + lbm_omega_w2_rho * (lbm_one_minus_u2 - lbm_30_uY + lbm_45_uY2);
  lbm_out_e = one_minus_omega * in_e + lbm_omega_w2_rho * (lbm_one_minus_u2 + lbm_30_uX + lbm_45_uX2);
  lbm_out_w = one_minus_omega * in_w + lbm_omega_w2_rho * (lbm_one_minus_u2 - lbm_30_uX + lbm_45_uX2);

  lbm_out_ne = one_minus_omega * in_ne + lbm_omega_w3_rho * (lbm_one_minus_u2 + lbm_30_uX_plus_uY + lbm_45_uX_plus_uY);
  lbm_out_nw = one_minus_omega * in_nw + lbm_omega_w3_rho * (lbm_one_minus_u2 - lbm_30_uX_minus_uY + lbm_45_uX_minus_uY);
  lbm_out_se = one_minus_omega * in_se + lbm_omega_w3_rho * (lbm_one_minus_u2 + lbm_30_uX_minus_uY + lbm_45_uX_minus_uY);
  lbm_out_sw = one_minus_omega * in_sw + lbm_omega_w3_rho * (lbm_one_minus_u2 - lbm_30_uX_plus_uY + lbm_45_uX_plus_uY);

  *out_c = lbm_out_c;

  *out_n = lbm_out_n;
  *out_s = lbm_out_s;
  *out_e = lbm_out_e;
  *out_w = lbm_out_w;

  *out_ne = lbm_out_ne;
  *out_nw = lbm_out_nw;
  *out_se = lbm_out_se;
  *out_sw = lbm_out_sw;
}

void dumpVelocities(int t) {
  int y, x;
  double rho, uX, uY;

  for (y = 1; y < nY-1 ; y++) {
    for (x = 1; x < nX-1; x++) {
      /* Compute rho */
      rho = grid[t % 2][y - 1][x + 0][N] + grid[t % 2][y + 1][x + 0][S] +
            grid[t % 2][y + 0][x - 1][E] + grid[t % 2][y + 0][x + 1][W] +
            grid[t % 2][y - 1][x - 1][NE] + grid[t % 2][y - 1][x + 1][NW] +
            grid[t % 2][y + 1][x - 1][SE] + grid[t % 2][y + 1][x + 1][SW] +
            grid[t % 2][y + 0][x + 0][C];

      /* Compute velocity along x-axis */
      uX = (grid[t % 2][y + 1][x - 1][SE] + grid[t % 2][y + 0][x - 1][E] +
            grid[t % 2][y - 1][x - 1][NE]) -
           (grid[t % 2][y + 1][x + 1][SW] + grid[t % 2][y + 0][x + 1][W] +
            grid[t % 2][y - 1][x + 1][NW]);
      uX = uX / rho;

      /* Compute velocity along y-axis */
      uY = (grid[t % 2][y - 1][x + 1][NW] + grid[t % 2][y - 1][x + 0][N] +
            grid[t % 2][y - 1][x - 1][NE]) -
           (grid[t % 2][y + 1][x + 1][SW] + grid[t % 2][y + 1][x + 0][S] +
            grid[t % 2][y + 1][x - 1][SE]);
      uY = uY / rho;

      fprintf(stderr, "%0.6lf,%0.6lf,%0.6lf\n", rho, uX, uY);
    }
  }
}

int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y) {
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

  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  return x->tv_sec < y->tv_sec;
}

/* icc -O3 -fp-model precise -restrict -DDEBUG -DTIME ldc_d2q9.c -o ldc_d2q9.elf */
