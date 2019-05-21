/****************************************************************************/
/* Lid Driven Cavity(LDC) - 2D simulation using MRT-LBM with a D2Q9 lattice */
/*                                                                          */
/* Scheme: 2-Grid, Pull, Compute, Keep                                      */
/*                                                                          */
/* Lattice naming convention                                                */
/*   c6[NW]  c2[N]  c5[NE]                                                  */
/*   c3[W]   c0[C]  c1[E]                                                   */
/*   c7[SW]  c4[S]  c8[SE]                                                  */
/*                                                                          */
/* Author: Irshad Pananilath <pmirshad+code@gmail.com>                      */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <sys/time.h>

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
#define nTimesteps 20000  // Number of timesteps for the simulation
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
const double re = 1000;
const double uTop = 0.1;

/* Calculated parameters: Set in main */
double omega = 0.0;
double nu = 0.0;

/* Moment relaxation rates */
double s_0 = 0.0;
double s_1 = 0.0;
double s_2 = 0.0;
double s_3 = 0.0;
double s_4 = 0.0;
double s_5 = 0.0;
double s_6 = 0.0;
double s_7 = 0.0;
double s_8 = 0.0;

const double w1 = 4.0 / 9.0;
const double w2 = 1.0 / 9.0;
const double w3 = 1.0 / 36.0;

double grid[2][nY + 2 + 4][nX + 2 + 2][nK];

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
  double rho, uX, uY;
  double m_eq0, m_eq1, m_eq2, m_eq3, m_eq4, m_eq5, m_eq6, m_eq7, m_eq8;

  /* For timekeeping */
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  /* Compute values for global parameters */
  omega = 2.0 / ((6.0 * uTop * (nX - 0.5) / re) + 1.0);
  nu = (1.0 / omega - 0.5) / 3.0;

  s_0 = 1.0;
  s_1 = 1.4;
  s_2 = 1.4;
  s_3 = 1.0;
  s_4 = 1.2;
  s_5 = 1.0;
  s_6 = 1.2;
  /* s_7 = omega; */
  /* s_8 = omega; */
  s_7 = 2.0 / (1.0 + 6.0 * nu);
  s_8 = 2.0 / (1.0 + 6.0 * nu);

  printf(
      "2D Lid Driven Cavity simulation using MRT-LBM with D2Q9 lattice:\n"
      "\tscheme     : 2-Grid, Fused, Pull\n"
      "\tgrid size  : %d x %d = %.2lf * 10^3 Cells\n"
      "\tnTimeSteps : %d\n"
      "\tRe         : %.2lf\n"
      "\tuTop       : %.6lf\n"
      "\tomega      : %.6lf\n",
      nX, nY, nX * nY / 1.0e3, nTimesteps, re, uTop, omega);

  /* Initialize the PDFs to equilibrium values */
  /* Calculate initial equilibrium moments */
  rho = 1.0;
  uX = 0.0;
  uY = 0.0;

  m_eq0 = rho;
  m_eq1 = -2.0 * rho + 3.0 * (uX * uX + uY * uY);
  m_eq2 = rho - 3.0 * (uX * uX + uY * uY);
  m_eq3 = uX;
  m_eq4 = -uX;
  m_eq5 = uY;
  m_eq6 = -uY;
  m_eq7 = uX * uX - uY * uY;
  m_eq8 = uX * uY;

  for (y = 0; y < nY + 2 + 4; y++) {
    for (x = 0; x < nX + 2 + 2; x++) {
      grid[0][y][x][0] = w1;
      grid[1][y][x][0] = w1;

      for (k = 1; k < 5; k++) {
        grid[0][y][x][k] = w2;
        grid[1][y][x][k] = w2;
      }

      for (k = 5; k < nK; k++) {
        grid[0][y][x][k] = w3;
        grid[1][y][x][k] = w3;
      }
    }
  }

  /* To satisfy PET */
  short _nX = nX + 3;
  short _nY = nY + 4;
  int _nTimesteps = nTimesteps;

#ifdef TIME
  gettimeofday(&start, 0);
#endif

#pragma scop
  /* Kernel of the code */
  for (t = 0; t < _nTimesteps; t++) {
    for (y = 2; y < _nY; y++) {
      for (x = 1; x < _nX; x++) {
        lbm_kernel(grid[t % 2][y][x][C],
                   grid[t % 2][y - 1][x + 0][N],
                   grid[t % 2][y + 1][x + 0][S],
                   grid[t % 2][y + 0][x - 1][E],
                   grid[t % 2][y + 0][x + 1][W],
                   grid[t % 2][y - 1][x - 1][NE],
                   grid[t % 2][y - 1][x + 1][NW],
                   grid[t % 2][y + 1][x - 1][SE],
                   grid[t % 2][y + 1][x + 1][SW],
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
  double rho, uX, uY, u2;
  double m_eq0, m_eq1, m_eq2, m_eq3, m_eq4, m_eq5, m_eq6, m_eq7, m_eq8;
  double m_0, m_1, m_2, m_3, m_4, m_5, m_6, m_7, m_8;

  if (((y == 2 || y == nY + 3) && (x > 1 && x < nX + 2)) ||
      ((x == 1 || x == nX + 2) && (y > 2 && y < nY + 3))) {
    /* Bounce back boundary condition at walls */
    *out_c = in_c;

    *out_s = in_n;
    *out_n = in_s;
    *out_w = in_e;
    *out_e = in_w;

    *out_sw = in_ne;
    *out_se = in_nw;
    *out_nw = in_se;
    *out_ne = in_sw;

    return;
  }

  /* Compute rho */
  rho = in_n + in_s + in_e + in_w + in_ne + in_nw + in_se + in_sw + in_c;

  /* Compute velocity along x-axis */
  uX = (in_se + in_e + in_ne) - (in_sw + in_w + in_nw);
  /* uX = uX / rho; */

  /* Compute velocity along y-axis */
  uY = (in_nw + in_n + in_ne) - (in_sw + in_s + in_se);
  /* uY = uY / rho; */

  /* Impose the constantly moving upper wall */
  if (((y == nY + 2) && (x > 1 && x < nX + 2))) {
    uX = uTop * rho;
    uY = 0.00 * rho;
  }

  /* Convert PDFs to moment space */
  m_0 = +1 * in_c + 1 * in_e + 1 * in_n + 1 * in_w + 1 * in_s + 1 * in_ne + 1 * in_nw + 1 * in_sw + 1 * in_se;
  m_1 = -1 * in_c - 1 * in_e - 1 * in_n - 1 * in_w + 2 * in_s + 2 * in_ne + 2 * in_nw + 2 * in_sw - 4 * in_se;
  m_2 = -2 * in_c - 2 * in_e - 2 * in_n - 2 * in_w + 1 * in_s + 1 * in_ne + 1 * in_nw + 1 * in_sw + 4 * in_se;
  m_3 = +1 * in_c + 0 * in_e - 1 * in_n + 0 * in_w + 1 * in_s - 1 * in_ne - 1 * in_nw + 1 * in_sw + 0 * in_se;
  m_4 = -2 * in_c + 0 * in_e + 2 * in_n + 0 * in_w + 1 * in_s - 1 * in_ne - 1 * in_nw + 1 * in_sw + 0 * in_se;
  m_5 = +0 * in_c + 1 * in_e + 0 * in_n - 1 * in_w + 1 * in_s + 1 * in_ne - 1 * in_nw - 1 * in_sw + 0 * in_se;
  m_6 = +0 * in_c - 2 * in_e + 0 * in_n + 2 * in_w + 1 * in_s + 1 * in_ne - 1 * in_nw - 1 * in_sw + 0 * in_se;
  m_7 = +1 * in_c - 1 * in_e + 1 * in_n - 1 * in_w + 0 * in_s + 0 * in_ne + 0 * in_nw + 0 * in_sw + 0 * in_se;
  m_8 = +0 * in_c + 0 * in_e + 0 * in_n + 0 * in_w + 1 * in_s - 1 * in_ne + 1 * in_nw - 1 * in_sw + 0 * in_se;

  /* Calculate equilibrium moments */
  m_eq0 = rho;
  m_eq1 = -2.0 * rho + 3.0 * (uX * uX + uY * uY);
  m_eq2 = rho - 3.0 * (uX * uX + uY * uY);
  m_eq3 = uX;
  m_eq4 = -uX;
  m_eq5 = uY;
  m_eq6 = -uY;
  m_eq7 = uX * uX - uY * uY;
  m_eq8 = uX * uY;

  /* Calculate difference, relax in momentum space and collide
   * (m - s * (m - m_eq)) */
  m_eq0 = m_0 - s_0 * (m_0 - m_eq0);
  m_eq1 = m_1 - s_1 * (m_1 - m_eq1);
  m_eq2 = m_2 - s_2 * (m_2 - m_eq2);
  m_eq3 = m_3 - s_3 * (m_3 - m_eq3);
  m_eq4 = m_4 - s_4 * (m_4 - m_eq4);
  m_eq5 = m_5 - s_5 * (m_5 - m_eq5);
  m_eq6 = m_6 - s_6 * (m_6 - m_eq6);
  m_eq7 = m_7 - s_7 * (m_7 - m_eq7);
  m_eq8 = m_8 - s_8 * (m_8 - m_eq8);

  /* Inverse transform from moment space */
  *out_c = (1.0 / 36.0) * (+4 * m_eq0 - 4 * m_eq1 + 4 * m_eq2 + 0 * m_eq3 + 0 * m_eq4 + 0 * m_eq5 + 0 * m_eq6 + 0 * m_eq7 + 0 * m_eq8);
  *out_e = (1.0 / 36.0) * (+4 * m_eq0 - 1 * m_eq1 - 2 * m_eq2 + 6 * m_eq3 - 6 * m_eq4 + 0 * m_eq5 + 0 * m_eq6 + 9 * m_eq7 + 0 * m_eq8);
  *out_n = (1.0 / 36.0) * (+4 * m_eq0 - 1 * m_eq1 - 2 * m_eq2 + 0 * m_eq3 - 0 * m_eq4 + 6 * m_eq5 - 6 * m_eq6 - 9 * m_eq7 + 0 * m_eq8);
  *out_w = (1.0 / 36.0) * (+4 * m_eq0 - 1 * m_eq1 - 2 * m_eq2 - 6 * m_eq3 + 6 * m_eq4 + 0 * m_eq5 + 0 * m_eq6 + 9 * m_eq7 + 0 * m_eq8);
  *out_s = (1.0 / 36.0) * (+4 * m_eq0 - 1 * m_eq1 - 2 * m_eq2 + 0 * m_eq3 + 0 * m_eq4 - 6 * m_eq5 + 6 * m_eq6 - 9 * m_eq7 + 0 * m_eq8);
  *out_ne = (1.0 / 36.0) * (+4 * m_eq0 + 2 * m_eq1 + 1 * m_eq2 + 6 * m_eq3 + 3 * m_eq4 + 6 * m_eq5 + 3 * m_eq6 + 0 * m_eq7 + 9 * m_eq8);
  *out_nw = (1.0 / 36.0) * (+4 * m_eq0 + 2 * m_eq1 + 1 * m_eq2 - 6 * m_eq3 - 3 * m_eq4 + 6 * m_eq5 + 3 * m_eq6 + 0 * m_eq7 - 9 * m_eq8);
  *out_sw = (1.0 / 36.0) * (+4 * m_eq0 + 2 * m_eq1 + 1 * m_eq2 - 6 * m_eq3 - 3 * m_eq4 - 6 * m_eq5 - 3 * m_eq6 + 0 * m_eq7 + 9 * m_eq8);
  *out_se = (1.0 / 36.0) * (+4 * m_eq0 + 2 * m_eq1 + 1 * m_eq2 + 6 * m_eq3 + 3 * m_eq4 - 6 * m_eq5 - 3 * m_eq6 + 0 * m_eq7 - 9 * m_eq8);
}

void dumpVelocities(int t) {
  int y, x;
  double rho, uX, uY;

  for (y = 3; y < nY + 3; y++) {
    for (x = 2; x < nX + 2; x++) {
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

/* icc -O3 -fp-model precise -restrict -DDEBUG -DTIME mrt_ldc_d2q9.c -o mrt_ldc_d2q9.elf */
