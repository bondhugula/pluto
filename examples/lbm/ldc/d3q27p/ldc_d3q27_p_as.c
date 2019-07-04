/*************************************************************************/
/* Lid Driven Cavity(LDC) - 3D simulation using LBM with a D3Q27 lattice */
/*                                                                       */
/* Scheme: 3-Grid, Pull, Compute, Keep                                   */
/*                                                                       */
/* Lattice naming convention                                             */
/*   http://arxiv.org/pdf/1202.6081.pdf                                  */
/*                                                                       */
/* Author: Irshad Pananilath <pmirshad+code@gmail.com>                   */
/*                                                                       */
/*************************************************************************/

#include <math.h>
#include <stdio.h>
#include <sys/time.h>

#ifndef nX
#define nX 200  // Size of domain along x-axis
#endif

#ifndef nY
#define nY 200  // Size of domain along y-axis
#endif

#ifndef nZ
#define nZ 200  // Size of domain along z-axis
#endif

#ifndef nK
#define nK 27  // D3Q27
#endif

#ifndef nTimesteps
#define nTimesteps 1000  // Number of timesteps for the simulation
#endif

#define C 0
#define E 1
#define W 2
#define N 3
#define S 4
#define T 5
#define B 6
#define NE 7
#define NW 8
#define SE 9
#define SW 10
#define TE 11
#define TW 12
#define BE 13
#define BW 14
#define TN 15
#define TS 16
#define BN 17
#define BS 18
#define TNE 19
#define TNW 20
#define TSE 21
#define TSW 22
#define BNE 23
#define BNW 24
#define BSE 25
#define BSW 26

/* Tunable parameters */
const double re = 100;
const double uTopX = 0.05;
const double uTopY = 0.02;

/* Calculated parameters: Set in main */
double omega = 0.0;
double one_minus_omega = 0.0;

const double w1 = 8.0 / 27.0;
const double w2 = 2.0 / 27.0;
const double w3 = 1.0 / 54.0;
const double w4 = 1.0 / 216.0;

double grid[2][nZ][nY][nX][nK] __attribute__((aligned(32)));

void lbm_kernel(double, double, double, double, double, double, double, double,
                double, double, double, double, double, double, double, double,
                double, double, double, double, double, double, double, double,
                double, double, double, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, double *restrict, double *restrict,
                double *restrict, int, int, int, int);

void dumpVelocities(int timestep);
int timeval_subtract(struct timeval *result, struct timeval *x,
                     struct timeval *y);

int main(void) {
  int t = 0, z, y, x, k;
  double total_lattice_pts =
      (double)nZ * (double)nY * (double)nX * (double)nTimesteps;

  /* For timekeeping */
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0;

  /* Compute values for global parameters */
  omega = 2.0 /
          ((6.0 * sqrt(uTopX * uTopX + uTopY * uTopY) * (nX - 0.5) / re) + 1.0);
  one_minus_omega = 1.0 - omega;

  printf(
      "3D Lid Driven Cavity simulation with D3Q27 lattice:\n"
      "\tfilename   : %s\n"
      "\tscheme     : 3-Grid, Fused, Pull\n"
      "\tgrid size  : %d x %d x %d = %.2lf * 10^3 Cells\n"
      "\tnTimeSteps : %d\n"
      "\tRe         : %.2lf\n"
      "\tuTopX      : %.6lf\n"
      "\tuTopY      : %.6lf\n"
      "\tomega      : %.6lf\n",
      __FILE__, nX, nY, nZ, nX * nY * nZ / 1.0e3, nTimesteps, re, uTopX, uTopY,
      omega);

  /* Initialize all 27 PDFs for each point in the domain to w1, w2 or w3
   * accordingly */
  for (z = 0; z < nZ; z++) {
    for (y = 0; y < nY; y++) {
      for (x = 0; x < nX; x++) {
        grid[0][z][y][x][0] = w1;
        grid[1][z][y][x][0] = w1;

        for (k = 1; k < 7; k++) {
          grid[0][z][y][x][k] = w2;
          grid[1][z][y][x][k] = w2;
        }

        for (k = 7; k < 19; k++) {
          grid[0][z][y][x][k] = w3;
          grid[1][z][y][x][k] = w3;
        }

        for (k = 19; k < nK; k++) {
          grid[0][z][y][x][k] = w4;
          grid[1][z][y][x][k] = w4;
        }
      }
    }
  }

  /* To satisfy PET */
  short _nX = nX;
  short _nY = nY;
  short _nZ = nZ;
  short _nTimesteps = nTimesteps;

#ifdef TIME
  gettimeofday(&start, 0);
#endif

#pragma scop
  /* Kernel of the code */
  for (t = 0; t < _nTimesteps; t++) {
    for (z = 0; z < _nZ; z++) {
      for (y = 0; y < _nY; y++) {
        for (x = 0; x < _nX; x++) {
          lbm_kernel(
              grid[t % 2][z + 0][y + 0][x + 0][C],
              grid[t % 2][z + 0][y + 0][x==0?_nX-1:x-1][E],
              grid[t % 2][z + 0][y + 0][x==_nX-1?0:x+1][W],
              grid[t % 2][z + 0][y==0?_nY-1:y-1][x + 0][N],
              grid[t % 2][z + 0][y==_nY-1?0:y+1][x + 0][S],
              grid[t % 2][z==0?_nZ-1:z-1][y + 0][x + 0][T],
              grid[t % 2][z==_nZ-1?0:z+1][y + 0][x + 0][B],
              grid[t % 2][z + 0][y==0?_nY-1:y-1][x==0?_nX-1:x-1][NE],
              grid[t % 2][z + 0][y==0?_nY-1:y-1][x==_nX-1?0:x+1][NW],
              grid[t % 2][z + 0][y==_nY-1?0:y+1][x==0?_nX-1:x-1][SE],
              grid[t % 2][z + 0][y==_nY-1?0:y+1][x==_nX-1?0:x+1][SW],
              grid[t % 2][z==0?_nZ-1:z-1][y + 0][x==0?_nX-1:x-1][TE],
              grid[t % 2][z==0?_nZ-1:z-1][y + 0][x==_nX-1?0:x+1][TW],
              grid[t % 2][z==_nZ-1?0:z+1][y + 0][x==0?_nX-1:x-1][BE],
              grid[t % 2][z==_nZ-1?0:z+1][y + 0][x==_nX-1?0:x+1][BW],
              grid[t % 2][z==0?_nZ-1:z-1][y==0?_nY-1:y-1][x + 0][TN],
              grid[t % 2][z==0?_nZ-1:z-1][y==_nY-1?0:y+1][x + 0][TS],
              grid[t % 2][z==_nZ-1?0:z+1][y==0?_nY-1:y-1][x + 0][BN],
              grid[t % 2][z==_nZ-1?0:z+1][y==_nY-1?0:y+1][x + 0][BS],
              grid[t % 2][z==0?_nZ-1:z-1][y==0?_nY-1:y-1][x==0?_nX-1:x-1][TNE],
              grid[t % 2][z==0?_nZ-1:z-1][y==0?_nY-1:y-1][x==_nX-1?0:x+1][TNW],
              grid[t % 2][z==0?_nZ-1:z-1][y==_nY-1?0:y+1][x==0?_nX-1:x-1][TSE],
              grid[t % 2][z==0?_nZ-1:z-1][y==_nY-1?0:y+1][x==_nX-1?0:x+1][TSW],
              grid[t % 2][z==_nZ-1?0:z+1][y==0?_nY-1:y-1][x==0?_nX-1:x-1][BNE],
              grid[t % 2][z==_nZ-1?0:z+1][y==0?_nY-1:y-1][x==_nX-1?0:x+1][BNW],
              grid[t % 2][z==_nZ-1?0:z+1][y==_nY-1?0:y+1][x==0?_nX-1:x-1][BSE],
              grid[t % 2][z==_nZ-1?0:z+1][y==_nY-1?0:y+1][x==_nX-1?0:x+1][BSW],
              &grid[(t + 1) % 2][z][y][x][C],
              &grid[(t + 1) % 2][z][y][x][E],
              &grid[(t + 1) % 2][z][y][x][W],
              &grid[(t + 1) % 2][z][y][x][N],
              &grid[(t + 1) % 2][z][y][x][S],
              &grid[(t + 1) % 2][z][y][x][T],
              &grid[(t + 1) % 2][z][y][x][B],
              &grid[(t + 1) % 2][z][y][x][NE],
              &grid[(t + 1) % 2][z][y][x][NW],
              &grid[(t + 1) % 2][z][y][x][SE],
              &grid[(t + 1) % 2][z][y][x][SW],
              &grid[(t + 1) % 2][z][y][x][TE],
              &grid[(t + 1) % 2][z][y][x][TW],
              &grid[(t + 1) % 2][z][y][x][BE],
              &grid[(t + 1) % 2][z][y][x][BW],
              &grid[(t + 1) % 2][z][y][x][TN],
              &grid[(t + 1) % 2][z][y][x][TS],
              &grid[(t + 1) % 2][z][y][x][BN],
              &grid[(t + 1) % 2][z][y][x][BS],
              &grid[(t + 1) % 2][z][y][x][TNE],
              &grid[(t + 1) % 2][z][y][x][TNW],
              &grid[(t + 1) % 2][z][y][x][TSE],
              &grid[(t + 1) % 2][z][y][x][TSW],
              &grid[(t + 1) % 2][z][y][x][BNE],
              &grid[(t + 1) % 2][z][y][x][BNW],
              &grid[(t + 1) % 2][z][y][x][BSE],
              &grid[(t + 1) % 2][z][y][x][BSW],
              t, z, y, x);
        }
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

void lbm_kernel(
    double in_c, double in_e, double in_w, double in_n, double in_s,
    double in_t, double in_b, double in_ne, double in_nw, double in_se,
    double in_sw, double in_te, double in_tw, double in_be, double in_bw,
    double in_tn, double in_ts, double in_bn, double in_bs, double in_tne,
    double in_tnw, double in_tse, double in_tsw, double in_bne, double in_bnw,
    double in_bse, double in_bsw, double *restrict out_c,
    double *restrict out_e, double *restrict out_w, double *restrict out_n,
    double *restrict out_s, double *restrict out_t, double *restrict out_b,
    double *restrict out_ne, double *restrict out_nw, double *restrict out_se,
    double *restrict out_sw, double *restrict out_te, double *restrict out_tw,
    double *restrict out_be, double *restrict out_bw, double *restrict out_tn,
    double *restrict out_ts, double *restrict out_bn, double *restrict out_bs,
    double *restrict out_tne, double *restrict out_tnw,
    double *restrict out_tse, double *restrict out_tsw,
    double *restrict out_bne, double *restrict out_bnw,
    double *restrict out_bse, double *restrict out_bsw, int t, int z, int y,
    int x) {
  /* %%INIT%% */
  double lbm_sum_e, lbm_sum_w, lbm_sum_n, lbm_sum_s, lbm_sum_t, lbm_sum_b,
      lbm_tne_plus_ne_plus_bne, lbm_tse_plus_se_plus_bse,
      lbm_tnw_plus_nw_plus_bnw, lbm_tsw_plus_sw_plus_bsw, lbm_te_plus_e_plus_be,
      lbm_tw_plus_w_plus_bw, lbm_tn_plus_n_plus_bn, lbm_ts_plus_s_plus_bs;

  double lbm_rho, lbm_inv_rho, lbm_omega_rho;
  double lbm_uX, lbm_uY, lbm_uZ;
  double lbm_uX2, lbm_uY2, lbm_uZ2;
  double lbm_u2, lbm_1_minus_u2;
  double lbm_omega_w1_rho, lbm_omega_w2_rho, lbm_omega_w3_rho, lbm_omega_w4_rho;
  double lbm_30_uX, lbm_45_uX2, lbm_30_uY, lbm_45_uY2, lbm_30_uZ, lbm_45_uZ2;
  double lbm_uX_plus_uY, lbm_uX_minus_uY, lbm_uY_plus_uZ, lbm_uY_minus_uZ,
      lbm_uX_plus_uZ, lbm_uX_minus_uZ, lbm_uX_plus_uY_plus_uZ,
      lbm_uX_minus_uY_plus_uZ, lbm_uX_plus_uY_minus_uZ,
      lbm_uX_minus_uY_minus_uZ;

  double lbm_30_uX_plus_uY, lbm_30_uX_minus_uY, lbm_30_uY_plus_uZ,
      lbm_30_uY_minus_uZ, lbm_30_uX_plus_uZ, lbm_30_uX_minus_uZ,
      lbm_45_uX_plus_uY, lbm_45_uX_minus_uY, lbm_45_uY_plus_uZ,
      lbm_45_uY_minus_uZ, lbm_45_uX_plus_uZ, lbm_45_uX_minus_uZ,
      lbm_30_uX_plus_uY_plus_uZ, lbm_30_uX_minus_uY_plus_uZ,
      lbm_30_uX_plus_uY_minus_uZ, lbm_30_uX_minus_uY_minus_uZ,
      lbm_45_uX_plus_uY_plus_uZ, lbm_45_uX_minus_uY_plus_uZ,
      lbm_45_uX_plus_uY_minus_uZ, lbm_45_uX_minus_uY_minus_uZ;

  double lbm_out_c, lbm_out_e, lbm_out_w, lbm_out_n, lbm_out_s, lbm_out_t,
      lbm_out_b, lbm_out_ne, lbm_out_nw, lbm_out_se, lbm_out_sw, lbm_out_tn,
      lbm_out_ts, lbm_out_te, lbm_out_tw, lbm_out_bn, lbm_out_bs, lbm_out_be,
      lbm_out_bw, lbm_out_tsw, lbm_out_tse, lbm_out_tnw, lbm_out_tne,
      lbm_out_bsw, lbm_out_bse, lbm_out_bnw, lbm_out_bne;

  double rho, uX, uY, uZ, u2;
  double temp;
  /* %%INIT%% */

  /* #<{(| %%BOUNDARY%% |)}># */
  /* if (((z == 2 || z == nZ + 3) && (x > 0 && x < nX + 3) && */
  /*      (y > 0 && y < nY + 3)) || */
  /*     ((x == 1 || x == nX + 2) && (y > 0 && y < nY + 3) && */
  /*      (z > 1 && z < nZ + 4)) || */
  /*     ((y == 1 || y == nY + 2) && (x > 0 && x < nX + 3) && */
  /*      (z > 1 && z < nZ + 4))) { */
  /*   #<{(| Bounce back boundary condition at walls |)}># */
  /*   *out_c = in_c; */
  /*  */
  /*   *out_s = in_n; */
  /*   *out_n = in_s; */
  /*   *out_w = in_e; */
  /*   *out_e = in_w; */
  /*   *out_b = in_t; */
  /*   *out_t = in_b; */
  /*  */
  /*   *out_sw = in_ne; */
  /*   *out_se = in_nw; */
  /*   *out_nw = in_se; */
  /*   *out_ne = in_sw; */
  /*  */
  /*   *out_bs = in_tn; */
  /*   *out_bn = in_ts; */
  /*   *out_bw = in_te; */
  /*   *out_be = in_tw; */
  /*   *out_ts = in_bn; */
  /*   *out_tn = in_bs; */
  /*   *out_tw = in_be; */
  /*   *out_te = in_bw; */
  /*  */
  /*   *out_bsw = in_tne; */
  /*   *out_bse = in_tnw; */
  /*   *out_bnw = in_tse; */
  /*   *out_bne = in_tsw; */
  /*   *out_tsw = in_bne; */
  /*   *out_tse = in_bnw; */
  /*   *out_tnw = in_bse; */
  /*   *out_tne = in_bsw; */
  /*  */
  /*   return; */
  /* } */
  /* #<{(| %%BOUNDARY%% |)}># */

  /* Algebraic Simplification Begin */

  lbm_tne_plus_ne_plus_bne =
      in_tne + in_ne + in_bne; /* partial sum for planar NE */
  lbm_tse_plus_se_plus_bse =
      in_tse + in_se + in_bse;                  /* partial sum for planar SE */
  lbm_te_plus_e_plus_be = in_te + in_e + in_be; /* partial sum for planar E */

  /* sum for volumetric E */
  lbm_sum_e = lbm_te_plus_e_plus_be + lbm_tne_plus_ne_plus_bne +
              lbm_tse_plus_se_plus_bse;

  lbm_tnw_plus_nw_plus_bnw =
      in_tnw + in_nw + in_bnw; /* partial sum for planar NW */
  lbm_tsw_plus_sw_plus_bsw =
      in_tsw + in_sw + in_bsw;                  /* partial sum for planar SW */
  lbm_tw_plus_w_plus_bw = in_tw + in_w + in_bw; /* partial sum for planar W */

  /* sum for volumetric W */
  lbm_sum_w = lbm_tw_plus_w_plus_bw + lbm_tnw_plus_nw_plus_bnw +
              lbm_tsw_plus_sw_plus_bsw;

  lbm_tn_plus_n_plus_bn = in_tn + in_n + in_bn; /* partial sum for planar N */
  lbm_ts_plus_s_plus_bs = in_ts + in_s + in_bs; /* partial sum for planar S */

  /* sum for volumetric N */
  lbm_sum_n = lbm_tn_plus_n_plus_bn + lbm_tne_plus_ne_plus_bne +
              lbm_tnw_plus_nw_plus_bnw;

  /* sum for volumetric S */
  lbm_sum_s = lbm_ts_plus_s_plus_bs + lbm_tse_plus_se_plus_bse +
              lbm_tsw_plus_sw_plus_bsw;

  /* 24 adds */

  /* Compute rho */
  lbm_rho = lbm_sum_e + lbm_sum_w + lbm_ts_plus_s_plus_bs +
            lbm_tn_plus_n_plus_bn + in_t + in_c + in_b;

  lbm_inv_rho = 1.0 / lbm_rho;

  /* Compute velocity along x-axis */
  lbm_uX = lbm_inv_rho * (lbm_sum_e - lbm_sum_w);

  /* Compute velocity along y-axis */
  lbm_uY = lbm_inv_rho * (lbm_sum_n - lbm_sum_s);

  /* Compute velocity along z-axis */
  lbm_uZ = lbm_inv_rho * ((in_tn + in_ts + in_te + in_tw + in_tne + in_tnw +
                           in_tse + in_tsw + in_t) -
                          (in_bn + in_bs + in_be + in_bw + in_bne + in_bnw +
                           in_bse + in_bsw + in_b));

  /* 49 add/sub in total + 3 mul */
  /* 28/77 add/sub (36.36%) eliminated from naive code */
  /* TODO: Should check whether Quadrant-wise addition is better */

  /* %%BOUNDARY%% */
  /* Impose the constantly moving upper wall */
  if (z == nZ -1) {
    lbm_uX = uTopX;
    lbm_uY = uTopY;
    lbm_uZ = 0.00;
  }
  /* %%BOUNDARY%% */

  lbm_uX2 = lbm_uX * lbm_uX;
  lbm_uY2 = lbm_uY * lbm_uY;
  lbm_uZ2 = lbm_uZ * lbm_uZ;

  lbm_u2 = 1.5 * (lbm_uX2 + lbm_uY2 + lbm_uZ2);
  lbm_1_minus_u2 = 1.0 - lbm_u2;

  lbm_omega_rho = omega * lbm_rho;

  lbm_omega_w1_rho = w1 * lbm_omega_rho;
  lbm_omega_w2_rho = w2 * lbm_omega_rho;
  lbm_omega_w3_rho = w3 * lbm_omega_rho;
  lbm_omega_w4_rho = w4 * lbm_omega_rho;

  lbm_30_uX = 3.0 * lbm_uX;
  lbm_45_uX2 = 4.5 * lbm_uX2;
  lbm_30_uY = 3.0 * lbm_uY;
  lbm_45_uY2 = 4.5 * lbm_uY2;
  lbm_30_uZ = 3.0 * lbm_uZ;
  lbm_45_uZ2 = 4.5 * lbm_uZ2;

  lbm_uX_plus_uY = lbm_uX + lbm_uY;
  lbm_uX_minus_uY = lbm_uX - lbm_uY;
  lbm_uY_plus_uZ = lbm_uY + lbm_uZ;
  lbm_uY_minus_uZ = lbm_uY - lbm_uZ;
  lbm_uX_plus_uZ = lbm_uX + lbm_uZ;
  lbm_uX_minus_uZ = lbm_uX - lbm_uZ;

  lbm_30_uX_plus_uY = 3.0 * lbm_uX_plus_uY;
  lbm_30_uX_minus_uY = 3.0 * lbm_uX_minus_uY;
  lbm_30_uY_plus_uZ = 3.0 * lbm_uY_plus_uZ;
  lbm_30_uY_minus_uZ = 3.0 * lbm_uY_minus_uZ;
  lbm_30_uX_plus_uZ = 3.0 * lbm_uX_plus_uZ;
  lbm_30_uX_minus_uZ = 3.0 * lbm_uX_minus_uZ;

  lbm_45_uX_plus_uY = 4.5 * lbm_uX_plus_uY * lbm_uX_plus_uY;
  lbm_45_uX_minus_uY = 4.5 * lbm_uX_minus_uY * lbm_uX_minus_uY;
  lbm_45_uY_plus_uZ = 4.5 * lbm_uY_plus_uZ * lbm_uY_plus_uZ;
  lbm_45_uY_minus_uZ = 4.5 * lbm_uY_minus_uZ * lbm_uY_minus_uZ;
  lbm_45_uX_plus_uZ = 4.5 * lbm_uX_plus_uZ * lbm_uX_plus_uZ;
  lbm_45_uX_minus_uZ = 4.5 * lbm_uX_minus_uZ * lbm_uX_minus_uZ;

  lbm_uX_plus_uY_plus_uZ = lbm_uX_plus_uY + lbm_uZ;
  lbm_uX_minus_uY_plus_uZ = lbm_uX_minus_uY + lbm_uZ;
  lbm_uX_plus_uY_minus_uZ = lbm_uX_plus_uY - lbm_uZ;
  lbm_uX_minus_uY_minus_uZ = lbm_uX_minus_uY - lbm_uZ;

  lbm_30_uX_plus_uY_plus_uZ = 3.0 * lbm_uX_plus_uY_plus_uZ;
  lbm_30_uX_minus_uY_plus_uZ = 3.0 * lbm_uX_minus_uY_plus_uZ;
  lbm_30_uX_plus_uY_minus_uZ = 3.0 * lbm_uX_plus_uY_minus_uZ;
  lbm_30_uX_minus_uY_minus_uZ = 3.0 * lbm_uX_minus_uY_minus_uZ;

  lbm_45_uX_plus_uY_plus_uZ =
      4.5 * lbm_uX_plus_uY_plus_uZ * lbm_uX_plus_uY_plus_uZ;
  lbm_45_uX_minus_uY_plus_uZ =
      4.5 * lbm_uX_minus_uY_plus_uZ * lbm_uX_minus_uY_plus_uZ;
  lbm_45_uX_plus_uY_minus_uZ =
      4.5 * lbm_uX_plus_uY_minus_uZ * lbm_uX_plus_uY_minus_uZ;
  lbm_45_uX_minus_uY_minus_uZ =
      4.5 * lbm_uX_minus_uY_minus_uZ * lbm_uX_minus_uY_minus_uZ;

  /* Compute and keep new PDFs */

  lbm_out_c = one_minus_omega * in_c + lbm_omega_w1_rho * (lbm_1_minus_u2);

  lbm_out_n = one_minus_omega * in_n +
              lbm_omega_w2_rho * (lbm_1_minus_u2 + lbm_30_uY + lbm_45_uY2);
  lbm_out_s = one_minus_omega * in_s +
              lbm_omega_w2_rho * (lbm_1_minus_u2 - lbm_30_uY + lbm_45_uY2);
  lbm_out_e = one_minus_omega * in_e +
              lbm_omega_w2_rho * (lbm_1_minus_u2 + lbm_30_uX + lbm_45_uX2);
  lbm_out_w = one_minus_omega * in_w +
              lbm_omega_w2_rho * (lbm_1_minus_u2 - lbm_30_uX + lbm_45_uX2);
  lbm_out_t = one_minus_omega * in_t +
              lbm_omega_w2_rho * (lbm_1_minus_u2 + lbm_30_uZ + lbm_45_uZ2);
  lbm_out_b = one_minus_omega * in_b +
              lbm_omega_w2_rho * (lbm_1_minus_u2 - lbm_30_uZ + lbm_45_uZ2);

  lbm_out_ne = one_minus_omega * in_ne +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 + lbm_30_uX_plus_uY + lbm_45_uX_plus_uY);
  lbm_out_nw = one_minus_omega * in_nw +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 - lbm_30_uX_minus_uY + lbm_45_uX_minus_uY);
  lbm_out_se = one_minus_omega * in_se +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 + lbm_30_uX_minus_uY + lbm_45_uX_minus_uY);
  lbm_out_sw = one_minus_omega * in_sw +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 - lbm_30_uX_plus_uY + lbm_45_uX_plus_uY);

  lbm_out_tn = one_minus_omega * in_tn +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 + lbm_30_uY_plus_uZ + lbm_45_uY_plus_uZ);
  lbm_out_ts = one_minus_omega * in_ts +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 - lbm_30_uY_minus_uZ + lbm_45_uY_minus_uZ);
  lbm_out_te = one_minus_omega * in_te +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 + lbm_30_uX_plus_uZ + lbm_45_uX_plus_uZ);
  lbm_out_tw = one_minus_omega * in_tw +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 - lbm_30_uX_minus_uZ + lbm_45_uX_minus_uZ);

  lbm_out_bn = one_minus_omega * in_bn +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 + lbm_30_uY_minus_uZ + lbm_45_uY_minus_uZ);
  lbm_out_bs = one_minus_omega * in_bs +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 - lbm_30_uY_plus_uZ + lbm_45_uY_plus_uZ);
  lbm_out_be = one_minus_omega * in_be +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 + lbm_30_uX_minus_uZ + lbm_45_uX_minus_uZ);
  lbm_out_bw = one_minus_omega * in_bw +
               lbm_omega_w3_rho *
                   (lbm_1_minus_u2 - lbm_30_uX_plus_uZ + lbm_45_uX_plus_uZ);

  lbm_out_tne = one_minus_omega * in_tne +
                lbm_omega_w4_rho * (lbm_1_minus_u2 + lbm_30_uX_plus_uY_plus_uZ +
                                    lbm_45_uX_plus_uY_plus_uZ);

  lbm_out_tnw =
      one_minus_omega * in_tnw +
      lbm_omega_w4_rho * (lbm_1_minus_u2 - lbm_30_uX_minus_uY_minus_uZ +
                          lbm_45_uX_minus_uY_minus_uZ);

  lbm_out_tse =
      one_minus_omega * in_tse +
      lbm_omega_w4_rho * (lbm_1_minus_u2 + lbm_30_uX_minus_uY_plus_uZ +
                          lbm_45_uX_minus_uY_plus_uZ);

  lbm_out_tsw =
      one_minus_omega * in_tsw +
      lbm_omega_w4_rho * (lbm_1_minus_u2 - lbm_30_uX_plus_uY_minus_uZ +
                          lbm_45_uX_plus_uY_minus_uZ);

  lbm_out_bne =
      one_minus_omega * in_bne +
      lbm_omega_w4_rho * (lbm_1_minus_u2 + lbm_30_uX_plus_uY_minus_uZ +
                          lbm_45_uX_plus_uY_minus_uZ);

  lbm_out_bnw =
      one_minus_omega * in_bnw +
      lbm_omega_w4_rho * (lbm_1_minus_u2 - lbm_30_uX_minus_uY_plus_uZ +
                          lbm_45_uX_minus_uY_plus_uZ);

  lbm_out_bse =
      one_minus_omega * in_bse +
      lbm_omega_w4_rho * (lbm_1_minus_u2 + lbm_30_uX_minus_uY_minus_uZ +
                          lbm_45_uX_minus_uY_minus_uZ);

  lbm_out_bsw = one_minus_omega * in_bsw +
                lbm_omega_w4_rho * (lbm_1_minus_u2 - lbm_30_uX_plus_uY_plus_uZ +
                                    lbm_45_uX_plus_uY_plus_uZ);

  *out_c = lbm_out_c;
  *out_n = lbm_out_n;
  *out_s = lbm_out_s;
  *out_e = lbm_out_e;
  *out_w = lbm_out_w;
  *out_t = lbm_out_t;
  *out_b = lbm_out_b;

  *out_ne = lbm_out_ne;
  *out_nw = lbm_out_nw;
  *out_se = lbm_out_se;
  *out_sw = lbm_out_sw;

  *out_tn = lbm_out_tn;
  *out_ts = lbm_out_ts;
  *out_te = lbm_out_te;
  *out_tw = lbm_out_tw;

  *out_bn = lbm_out_bn;
  *out_bs = lbm_out_bs;
  *out_be = lbm_out_be;
  *out_bw = lbm_out_bw;

  *out_tne = lbm_out_tne;
  *out_tnw = lbm_out_tnw;
  *out_tse = lbm_out_tse;
  *out_tsw = lbm_out_tsw;
  *out_bne = lbm_out_bne;
  *out_bnw = lbm_out_bnw;
  *out_bse = lbm_out_bse;
  *out_bsw = lbm_out_bsw;

  /* Algebraic Simplification End */
}

void dumpVelocities(int t) {
  int z, y, x;
  double rho, uX, uY, uZ;

  for (z = 1; z < nZ - 1; z++) {
    for (y = 1; y < nY - 1; y++) {
      for (x = 1; x < nX - 1; x++) {
        /* Compute rho */
        rho = grid[t % 2][z + 0][y - 1][x + 0][N] +
              grid[t % 2][z + 0][y + 1][x + 0][S] +
              grid[t % 2][z + 0][y + 0][x - 1][E] +
              grid[t % 2][z + 0][y + 0][x + 1][W] +
              grid[t % 2][z - 1][y + 0][x + 0][T] +
              grid[t % 2][z + 1][y + 0][x + 0][B] +
              grid[t % 2][z + 0][y - 1][x - 1][NE] +
              grid[t % 2][z + 0][y - 1][x + 1][NW] +
              grid[t % 2][z + 0][y + 1][x - 1][SE] +
              grid[t % 2][z + 0][y + 1][x + 1][SW] +
              grid[t % 2][z - 1][y - 1][x + 0][TN] +
              grid[t % 2][z - 1][y + 1][x + 0][TS] +
              grid[t % 2][z - 1][y + 0][x - 1][TE] +
              grid[t % 2][z - 1][y + 0][x + 1][TW] +
              grid[t % 2][z + 1][y - 1][x + 0][BN] +
              grid[t % 2][z + 1][y + 1][x + 0][BS] +
              grid[t % 2][z + 1][y + 0][x - 1][BE] +
              grid[t % 2][z + 1][y + 0][x + 1][BW] +
              grid[t % 2][z - 1][y - 1][x - 1][TNE] +
              grid[t % 2][z - 1][y - 1][x + 1][TNW] +
              grid[t % 2][z - 1][y + 1][x - 1][TSE] +
              grid[t % 2][z - 1][y + 1][x + 1][TSW] +
              grid[t % 2][z + 1][y - 1][x - 1][BNE] +
              grid[t % 2][z + 1][y - 1][x + 1][BNW] +
              grid[t % 2][z + 1][y + 1][x - 1][BSE] +
              grid[t % 2][z + 1][y + 1][x + 1][BSW] +
              grid[t % 2][z + 0][y + 0][x + 0][C];

        /* Compute velocity along x-axis */
        uX = (grid[t % 2][z + 0][y + 1][x - 1][SE] +
              grid[t % 2][z + 0][y - 1][x - 1][NE] +
              grid[t % 2][z - 1][y + 0][x - 1][TE] +
              grid[t % 2][z + 1][y + 0][x - 1][BE] +
              grid[t % 2][z - 1][y - 1][x - 1][TNE] +
              grid[t % 2][z - 1][y + 1][x - 1][TSE] +
              grid[t % 2][z + 1][y - 1][x - 1][BNE] +
              grid[t % 2][z + 1][y + 1][x - 1][BSE] +
              grid[t % 2][z + 0][y + 0][x - 1][E]) -
             (grid[t % 2][z + 0][y + 1][x + 1][SW] +
              grid[t % 2][z + 0][y - 1][x + 1][NW] +
              grid[t % 2][z - 1][y + 0][x + 1][TW] +
              grid[t % 2][z + 1][y + 0][x + 1][BW] +
              grid[t % 2][z - 1][y - 1][x + 1][TNW] +
              grid[t % 2][z - 1][y + 1][x + 1][TSW] +
              grid[t % 2][z + 1][y - 1][x + 1][BNW] +
              grid[t % 2][z + 1][y + 1][x + 1][BSW] +
              grid[t % 2][z + 0][y + 0][x + 1][W]);
        uX = uX / rho;

        /* Compute velocity along y-axis */
        uY = (grid[t % 2][z + 0][y - 1][x + 1][NW] +
              grid[t % 2][z + 0][y - 1][x - 1][NE] +
              grid[t % 2][z - 1][y - 1][x + 0][TN] +
              grid[t % 2][z + 1][y - 1][x + 0][BN] +
              grid[t % 2][z - 1][y - 1][x - 1][TNE] +
              grid[t % 2][z - 1][y - 1][x + 1][TNW] +
              grid[t % 2][z + 1][y - 1][x - 1][BNE] +
              grid[t % 2][z + 1][y - 1][x + 1][BNW] +
              grid[t % 2][z + 0][y - 1][x + 0][N]) -
             (grid[t % 2][z + 0][y + 1][x + 1][SW] +
              grid[t % 2][z + 0][y + 1][x - 1][SE] +
              grid[t % 2][z - 1][y + 1][x + 0][TS] +
              grid[t % 2][z + 1][y + 1][x + 0][BS] +
              grid[t % 2][z - 1][y + 1][x - 1][TSE] +
              grid[t % 2][z - 1][y + 1][x + 1][TSW] +
              grid[t % 2][z + 1][y + 1][x - 1][BSE] +
              grid[t % 2][z + 1][y + 1][x + 1][BSW] +
              grid[t % 2][z + 0][y + 1][x + 0][S]);
        uY = uY / rho;

        /* Compute velocity along z-axis */
        uZ = (grid[t % 2][z - 1][y - 1][x + 0][TN] +
              grid[t % 2][z - 1][y + 1][x + 0][TS] +
              grid[t % 2][z - 1][y + 0][x - 1][TE] +
              grid[t % 2][z - 1][y + 0][x + 1][TW] +
              grid[t % 2][z - 1][y - 1][x - 1][TNE] +
              grid[t % 2][z - 1][y - 1][x + 1][TNW] +
              grid[t % 2][z - 1][y + 1][x - 1][TSE] +
              grid[t % 2][z - 1][y + 1][x + 1][TSW] +
              grid[t % 2][z - 1][y + 0][x + 0][T]) -
             (grid[t % 2][z + 1][y - 1][x + 0][BN] +
              grid[t % 2][z + 1][y + 1][x + 0][BS] +
              grid[t % 2][z + 1][y + 0][x - 1][BE] +
              grid[t % 2][z + 1][y + 0][x + 1][BW] +
              grid[t % 2][z + 1][y - 1][x - 1][BNE] +
              grid[t % 2][z + 1][y - 1][x + 1][BNW] +
              grid[t % 2][z + 1][y + 1][x - 1][BSE] +
              grid[t % 2][z + 1][y + 1][x + 1][BSW] +
              grid[t % 2][z + 1][y + 0][x + 0][B]);
        uZ = uZ / rho;

        fprintf(stderr, "%0.6lf,%0.6lf,%0.6lf,%0.6lf\n", rho, uX, uY, uZ);
      }
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

/* icc -O3 -fp-model precise -restrict -DDEBUG -DTIME ldc_d3q27_as.c -o
 * ldc_d3q27_as.elf */
