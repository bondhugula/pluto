//*********************************************************************
//
// sparse matrix based fdtd code
//
// on Tech-X machines, compile with:
//   g++  -O3 -I/home/research/messmer/soft/aztec/lib fdtd.c
//              /home/research/messmer/soft/aztec/lib/libaztec.a -lm -lg2c -lblas
// or
//    mpicxx  -DAZ_MPI -O3 -I/usr/local/aztecmpi/lib fdtd.c
//             -L/usr/local/aztecmpi/lib -laztec -lm -lg2c -lblas
//
//*********************************************************************

#include <stdio.h>
#include <stdlib.h>
#include "az_aztec.h"
#include "math.h"

#ifdef AZ_MPI
#include "mpi.h"
#endif

#define perror(str) { fprintf(stderr,"%s\n",str);   exit(-1); }

const int ncomp = 3;            // number of components in each field
const int compX = 0;            // index of the X component
const int compY = 1;            // index of the Y component
const int compZ = 2;            // index of the Z component
/*
const int nx = 128;             // domain size
const int ny = 128;
const int nz = 128;
const int nsteps = 100;         // number of time steps
*/
const int nx = 20;             // domain size
const int ny = 40;
const int nz = 1;
const int nsteps = 1000;         // number of time steps


const double dx = 0.5;          // grid spacing
const double dy = 0.5;
const double dz = 0.5;
const double dt  = 1e-10;       // timestep

const double ox = nx/2 * dx;
const double oy = ny/2 * dy;
const double oz = nz/2 * dz;

const double R = nx/2*dx*0.45;

const double c   = 2.9979e8;    // speed of light
const double c2   = c * c;
const double omega = 5.e8;      // frequency of launched pulse

#define MAX_NZ_ROW 5 /*  Max number of nonzero elements in any matrix row    */

void create_curl_matrix_row_edge(int row,int i,double val[],int bindx[]);
void create_curl_matrix_row_face(int row,int i,double val[],int bindx[]);
void setup_bc(int *update, int *update_index, int N_update,
              int **bc_indx, int *n_bc);
int pos_to_row(int x, int y, int z);
int row_to_pos(int row, int *comp, int *x, int *y, int *z);
void write_file(char *filename, double *tmp_vec, int N_update);
long currentTimeMillis();
	

int main(int argc, char *argv[])
{

  /* See Aztec User's Guide for the variables that follow:         */
  int    proc_config[AZ_PROC_SIZE];/* Processor information.                */
  int    N_update;                 /* # of unknowns updated on this node    */
  int    *update;                  /* vector elements updated on this node  */

  int    *data_orgA;               /* Array to specify data layout          */
  int    *externalA;               /* vector elements needed by this node.  */
  int    *update_indexA;           /* ordering of update[] and external[]   */
  int    *extern_indexA;           /* locally on this processor.            */
  int    *bindxA;                  /* Sparse matrix to be solved is stored  */
  double *valA;                    /* in these MSR arrays.                  */
  AZ_MATRIX *mat_curl_edge;        /* curl operator matrix                  */

  int    *data_orgB;                /* Array to specify data layout          */
  int    *externalB;                /* vector elements needed by this node.  */
  int    *update_indexB;            /* ordering of update[] and external[]   */
  int    *extern_indexB;            /* locally on this processor.            */
  int    *bindxB;                   /* Sparse matrix to be solved is stored  */
  double *valB;                     /* in these MSR arrays.                  */
  AZ_MATRIX *mat_curl_face;         /* curl operator matrix                  */

  int *bc_indx;
  int n_bc;

  double *efield;
  double *bfield;
  double *epsilon;
  double *tmp_vec;
  double *tmp_vec2;

  int    i, nrow, x, y, z;
  int k, t;
  long startTime, endTime;
  int myrank;
  int vec_len;


  /* get number of processors and the name of this processor */
#ifdef AZ_MPI
  MPI_Init(&argc,&argv);
  AZ_set_proc_config(proc_config, MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
  myrank = 0;
  AZ_set_proc_config(proc_config, AZ_NOT_MPI);
#endif

  nrow = ncomp * nx * ny * nz;  /* overll number of matrix rows  */

  // Define partitioning:  matrix rows (ascending order) owned by this node
  // Here it is done automatically, but it can also be specified by hand
  AZ_read_update(&N_update, &update, proc_config, nrow, 1, AZ_linear);


  // In the following we set up the matrix for the edge centered curl operator
  // All the steps are described in detail in the AZTEC manual.
  // first: allocate space for the first matrix.
  bindxA = (int    *) malloc((N_update*MAX_NZ_ROW+1)*sizeof(int));
  valA   = (double *) malloc((N_update*MAX_NZ_ROW+1)*sizeof(double));
  if (valA == NULL) perror("Error: Not enough space to create matrix");
  // Initialize the index for the first off diagonal element
  bindxA[0] = N_update+1;

  // Create the matrix row by row. Each processor creates only rows appearing
  // in update[] (using global col. numbers).
  for (i = 0; i < N_update; i++)
    create_curl_matrix_row_edge(update[i], i, valA, bindxA);

  // convert matrix to a local distributed matrix
  AZ_transform(proc_config, &externalA, bindxA, valA, update, &update_indexA,
               &extern_indexA, &data_orgA, N_update, NULL, NULL, NULL, NULL,
               AZ_MSR_MATRIX);

  // convert the matrix arrays into a matrix structure, used in the
  // matrix vector multiplication
  mat_curl_edge = AZ_matrix_create(data_orgA[AZ_N_internal] + data_orgA[AZ_N_border]);
  AZ_set_MSR(mat_curl_edge, bindxA, valA, data_orgA, 0, NULL, AZ_LOCAL);

  // at this point the edge centered curl matrix is completed.

  // In the following we set up the matrix for the face centered curl operator
  // All the steps are described in detail in the AZTEC manual.
  // first: allocate space for the first matrix.
  bindxB = (int    *) malloc((N_update*MAX_NZ_ROW+1)*sizeof(int));
  valB   = (double *) malloc((N_update*MAX_NZ_ROW+1)*sizeof(double));
  if (valB == NULL) perror("Error: Not enough space to create matrix");

  // Initialize the index for the first off diagonal element
  bindxB[0] = N_update+1;

  // Create the matrix row by row. Each processor creates only rows appearing
  // in update[] (using global col. numbers).
  for (i = 0; i < N_update; i++)
      create_curl_matrix_row_face(update[i], i, valB, bindxB);

  // convert matrix to a local distributed matrix
  AZ_transform(proc_config, &externalB, bindxB, valB, update, &update_indexB,
			                   &extern_indexB, &data_orgB,
					   N_update, NULL, NULL, NULL, NULL,
					                  AZ_MSR_MATRIX);
  // convert the matrix arrays into a matrix structure, used in the
  // matrix vector multiplication
  mat_curl_face = AZ_matrix_create(data_orgB[AZ_N_internal] + data_orgB[AZ_N_border]);
  AZ_set_MSR(mat_curl_face, bindxB, valB, data_orgB, 0, NULL, AZ_LOCAL);

  // at this point the face centered curl matrix is completed.


  //  allocate memory for the fields and a temporary vector
  vec_len = N_update + data_orgA[AZ_N_external];
  efield = (double *) malloc(vec_len*sizeof(double));
  bfield = (double *) malloc(vec_len*sizeof(double));
  epsilon = (double *) malloc(vec_len*sizeof(double));
  tmp_vec = (double *) malloc(vec_len*sizeof(double));
  tmp_vec2 = (double *) malloc(vec_len*sizeof(double));

  // setup the boundary condition. We will get an arry that tells us
  // which positions need to be updated and where the results needs
  // to be stored in the E field.
  setup_bc(update, update_indexB, N_update, &bc_indx, &n_bc);

  // initialize the field vectors
  for(k = 0; k < vec_len; k++){
     efield[k] = 0.;
     bfield[k] = 0.;
     epsilon[k] = 1.;
     tmp_vec[k] = 0.;
  }

  // initialize the dielectric structure. Ugly hard-coded stuff,
  // needs to be cleaned out...
  for(y=45; y<55; y++){
   for(x = y; x<100; x++)
      epsilon[compZ + pos_to_row(x, y, 0)] = 0.95;
  }	
  // reorder the dielectric vector in order to align with the B field
  AZ_reorder_vec(epsilon, data_orgA, update_indexA, NULL);


  printf("Begin iteration \n");

  // just some timing ...
  startTime = currentTimeMillis();

  // *******************
  // begin of the time stepping loop
  // *******************

  for( t = 0; t < nsteps; t++){

    // first we do the e field update

    // convert the B field to the H field

    for(k = 0 ; k < vec_len; k++)
      bfield[k] *= epsilon[k];

    // setup the initial condition
    for( k = 0; k < n_bc; k++){
      x = bc_indx[4*k];
      y = bc_indx[4*k+1];
      z = bc_indx[4*k+2];
      efield[bc_indx[4*k+3]] =
           sin((double) y * 5. * 3.14159 / (double) ny) *
           sin(omega * dt * (double) (t + 1));
    }

    //  E field update:
    //  tmp_vec = Curl_Op * bfield
    //  efield = efield +  c^2 * dt * tmp_vec
    AZ_MSR_matvec_mult( bfield, tmp_vec, mat_curl_edge, proc_config);

    // reorder the result in tmp_vec so that it aligns with the
    // decomposition of the E field
    AZ_invorder_vec(tmp_vec, data_orgA, update_indexA, NULL, tmp_vec2);
    AZ_reorder_vec(tmp_vec2, data_orgB, update_indexB, NULL);

    // update the efield
    for(k = 0 ; k < N_update; k++)
      efield[k] = efield[k] + c2 * tmp_vec2[k] * dt;

    // bfield update :
    // tmp_vec = DualCurl_Op * efield
    // bfield = bfield - tmp_vec * dt
    AZ_MSR_matvec_mult( efield, tmp_vec, mat_curl_face, proc_config);

    // reorder the result so that it fits the decomposition of the bfield
    AZ_invorder_vec(tmp_vec, data_orgB, update_indexB, NULL, tmp_vec2);
    AZ_reorder_vec(tmp_vec2, data_orgA, update_indexA, NULL);

    // update the b field
    for(k = 0;  k < N_update; k++)
	  bfield[k] = bfield[k] - tmp_vec2[k] * dt;

    if(myrank == 0)
      printf("Taking step %d at time %g\n", t,
		       (double) (currentTimeMillis() - startTime) / 1000.);
  }
  // ******************
  // end of timestepping loop
  // *****************

  endTime = currentTimeMillis();
  printf("After iteration: %g\n", (double)(endTime - startTime) / 1000. );

#if 1
  system("rm efield.txt bfield.txt");

  // dump filed data: efield
  AZ_invorder_vec(efield, data_orgB, update_indexB, NULL, tmp_vec);
  write_file("efield.txt", tmp_vec, N_update);

  // dump filed data: bfield
  AZ_invorder_vec(bfield, data_orgA, update_indexA, NULL, tmp_vec);
  write_file("bfield.txt", tmp_vec, N_update);
#endif

  /* Free allocated memory */
  AZ_matrix_destroy( &mat_curl_edge);
  free((void *) update);   free((void *) update_indexA);
  free((void *) externalA); free((void *) extern_indexA);
  free((void *) bindxA);    free((void *) valA);  free((void *) data_orgA);

  AZ_matrix_destroy( &mat_curl_face);
  free((void *) externalB); free((void *) extern_indexB);
  free((void *) bindxB);    free((void *) valB);  free((void *) data_orgB);

  free((void *) efield);    free((void *) bfield);
  free((void *) tmp_vec);   free((void *) tmp_vec2);


#ifdef AZ_MPI
  MPI_Finalize();
#endif
  return(1);

}

//***************************************************************************
//***************************************************************************
//***************************************************************************

//
// determine the global index for a given position
//
int pos_to_row(int x, int y, int z){
 const int Lx =  ncomp * nx;
 const int Ly =  ncomp * nx * ny;
 return (((x + nx) % nx) * ncomp + ( (y + ny) % ny) * Lx + ((z+nz) % nz) * Ly);
}

//
// determine the position for a given matrix row
//
int row_to_pos(int row, int *comp, int *x, int *y, int *z){
  const int Lx = ncomp * nx;
  const int Ly = ncomp * nx * ny;
  const int Lz = ncomp * nx * ny * nz;

  *comp = row % ncomp;
  *z = row / Ly;
  *y = (row - (*z) * Ly) / Lx;
  *x = (row - (*z) * Ly - (*y) * Lx) / ncomp;
}

//
// setup the boundary condition.
//
void setup_bc(int *update, int *update_index,
    int N_update, int **bc_indx, int *n_bc){
 int x, y, z;
 int indx;
 int ymin = 4 * ny / 10;
 int ymax = 6 * ny / 10;
 int entry;

 x = 0;
 z = 0;

 *n_bc = 0;
#define MAX_BC ny

 *bc_indx = (int *) malloc( 4 * MAX_BC * sizeof(int));
 for(y = ymin; y < ymax;  y++){
  indx = pos_to_row(x, y, z) + 1;
  if((entry = AZ_find_index(indx, update, N_update)) != -1){
    (*bc_indx)[4 * (*n_bc)+0] = x;
    (*bc_indx)[4 * (*n_bc)+1] = y - ymin;
    (*bc_indx)[4 * (*n_bc)+2] = z;
    (*bc_indx)[4 * (*n_bc)+3] = update_index[entry];
    (*n_bc)++;
  }
 }

}

void getPhysCoord(int x, int y, int z, double* px, double* py, double* pz){
  *px = x * dx - ox;
  *py = y * dy - oy;
  *pz = z * dz - oz;
}

bool isInside(int x, int y, int z){
 double px, py, pz;
 getPhysCoord(x, y, z, &px, &py, &pz);
 return ((px*px + py*py + pz*pz) <= R*R );
}

// determine the intersection length for a pair of points where
// one point is inside and the other point is outside the geometry.
// this function is specific to the spherical cavity.
double getIntersect(int x1, int y1, int z1, int x2, int y2, int z2){
 double px1, py1, pz1, px2, py2, pz2;
 double dpx, dpy, dpz;
 double a, b, c, n1, n2, n;

 getPhysCoord(x1, y1, z1, &px1, &py1, &pz1);
 getPhysCoord(x2, y2, z2, &px2, &py2, &pz2);

 dpx = px2 - px1;
 dpy = py2 - py1;
 dpz = pz2 - pz1;

 a = dpx * dpx + dpy * dpy + dpz * dpz;
 b = 2 * px1 * dpx + 2 * py1 * dpy + 2 * pz1 * dpz;
 c = px1 * px1 + py1 * py1 + pz1 * pz1 - R * R;

 n1 = (-b + sqrt(b*b - 4. * a * c)) / (2. * a);
 n2 = (-b - sqrt(b*b - 4. * a * c)) / (2. * a);

// printf("pt1 = %g/%g/%g pt2=%g/%g/%g n1=%g n2=%g\n",
//	px1, py1, pz1, px2, py2, pz2, n1, n2);

 n = (n1 < 1) ? n1 : n2;

 return n;
}

double  getEdgeLen(int x1, int y1, int z1, int x2, int y2, int z2){
 bool inside1, inside2;
 double dl;
 inside1 = isInside(x1, y1, z1);
 inside2 = isInside(x2, y2, z2);
 if (!inside1 & !inside2) return 0;
 if ( inside1 &  inside2) return 1.;
 dl = getIntersect(x1, y1, z1, x2, y2, z2);
// printf("dl = %g inside1 = %d\n", dl, inside1);
 return (inside1) ? dl : 1. - dl;
}

double getArea(double dl1, double dl2, double dl3, double dl4,
		double l1, double l2 ){
 // center of cell in the cell.  	
 if((dl1 + dl2 > l1) || (dl3 + dl4 > l2) ||
   (dl1 == 0. && dl2 == 0. && dl3 == 0. && dl4 == 0.)){
	  return l1 * l2;
 }

 // trapeze
 if(dl1 == l1 || dl2 == l1){
   return l1 * ( ((dl3 < dl4) ? dl3 : dl4) + fabs(dl3 - dl4) / 2);
 }

 if(dl3 == l2 || dl4 == l2){
   return  l2 *( ((dl1 < dl2) ? dl1 : dl2) + fabs(dl1 - dl2) / 2);
 }

 // triangle
 return  (dl1 + dl2) * (dl3 + dl4) / 2;

}


//
// create the matrix for the curl operator. The edge flag indicates if
// the curl will result in a surface centered or edge centered result.
//
void create_curl_matrix_row_edge(int row, int location, double val[],
              int bindx[]){
  int k;
  int comp, x, y, z;
  int upper, lower;

  row_to_pos(row, &comp, &x, &y, &z);

  //  for edge centered curl (result is on the edge, e.g. for curl B in e update)
  //  we finite difference via f(x) - f(x-1). For the face centered curl (result
  //  is on the face, e.g. for curl E in b update), we finite difference via
  //  f(x+1) - f(x);

  upper =  0;
  lower = -1;

  k = bindx[location];

  switch(comp) {

  case 0: {
    /* dBz / dy */
    bindx[k] = compZ + pos_to_row(x, y+upper, z); val[k++] = +1. / dy;
    bindx[k] = compZ + pos_to_row(x, y+lower, z); val[k++] = -1. / dy;

    /* - dBy / dz */
    bindx[k] = compY + pos_to_row(x, y, z+upper); val[k++] = -1. / dz;
    bindx[k] = compY + pos_to_row(x, y, z+lower); val[k++] = +1. / dz;
    break;
  }

  case 1: {
    /* dBx / dz  */
    bindx[k] = compX + pos_to_row(x, y, z+upper); val[k++] = +1. / dz;
    bindx[k] = compX + pos_to_row(x, y, z+lower); val[k++] = -1. / dz;
    /* -dBz / dx  */
     bindx[k] = compZ + pos_to_row(x+upper, y, z); val[k++] = -1. / dx;
     bindx[k] = compZ + pos_to_row(x+lower, y, z); val[k++] = +1. / dx;
     break;
  }

  case 2: {
    /* dBy/dx  */
    bindx[k] = compY + pos_to_row(x+upper, y, z); val[k++] = +1. / dx;
    bindx[k] = compY + pos_to_row(x+lower, y, z); val[k++] = -1. / dx;

    /* -dBx/dy  */
    bindx[k] = compX + pos_to_row(x, y+upper, z); val[k++] = -1. / dy;
    bindx[k] = compX + pos_to_row(x, y+lower, z); val[k++] = +1. / dy;
    break;
  }
  }

  bindx[location+1] = k;  val[location]     = 0.; /* matrix diagonal */
}



//
// create the matrix for the curl operator. The edge flag indicates if
// the curl will result in a surface centered or edge centered result.
//
void create_curl_matrix_row_face(int row, int location, double val[],
              int bindx[]){
  int k;
  int comp, x, y, z;
  int upper, lower;
  double dl1, dl2, dl3, dl4, dA;

  row_to_pos(row, &comp, &x, &y, &z);

  //  for edge centered curl (result is on the edge, e.g. for curl B in e update)
  //  we finite difference via f(x) - f(x-1). For the face centered curl (result
  //  is on the face, e.g. for curl E in b update), we finite difference via
  //  f(x+1) - f(x);

  upper = +1;
  lower =  0;

  k = bindx[location];

  switch(comp) {

  case 0: {
    dl1 = dz * getEdgeLen(x, y+upper, z, x, y+upper, z+1);
    dl2 = dz * getEdgeLen(x, y+lower, z, x, y+lower, z+1);
    dl3 = dy * getEdgeLen(x, y, z+upper, x, y+1, z+upper);
    dl4 = dy * getEdgeLen(x, y, z+lower, x, y+1, z+lower);

    dA = getArea(dl1, dl2, dl3, dl4, dz, dy);
//    if(dA == 0) printf("dA(0) = %g dl1=%g dl2=%g dl3=%g dl4=%g\n", dA, dl1, dl2, dl3, dl4);

    /* dEz / dy */
    bindx[k] = compZ + pos_to_row(x, y+upper, z); val[k++] = +dl1 / dA;
    bindx[k] = compZ + pos_to_row(x, y+lower, z); val[k++] = -dl2 / dA;

    /* - dEy / dz */
    bindx[k] = compY + pos_to_row(x, y, z+upper); val[k++] = -dl3 / dA;
    bindx[k] = compY + pos_to_row(x, y, z+lower); val[k++] = +dl4 / dA;
    break;
  }

  case 1: {
    dl1 = dx * getEdgeLen(x, y, z+upper, x+1, y, z+upper);
    dl2 = dx * getEdgeLen(x, y, z+lower, x+1, y, z+lower);
    dl3 = dz * getEdgeLen(x+upper, y, z, x+upper, y, z+1);
    dl4 = dz * getEdgeLen(x+lower, y, z, x+lower, y, z+1);

    dA = getArea(dl1, dl2, dl3, dl4, dx, dz);
 //   if(dA == 0) printf("dA(1) = %g dl1=%g dl2=%g dl3=%g dl4=%g\n", dA, dl1, dl2, dl3, dl4);

    /*dEx / dz  */
    bindx[k] = compX + pos_to_row(x, y, z+upper); val[k++] = +dl1 / dA;
    bindx[k] = compX + pos_to_row(x, y, z+lower); val[k++] = -dl2 / dA;
    /* -dEz / dx  */
     bindx[k] = compZ + pos_to_row(x+upper, y, z); val[k++] = -dl3 / dA;
     bindx[k] = compZ + pos_to_row(x+lower, y, z); val[k++] = +dl4 / dA;
     break;
  }

  case 2: {
    dl1 = dy * getEdgeLen(x+upper, y, z, x+upper, y+1, z);
    dl2 = dy * getEdgeLen(x+lower, y, z, x+lower, y+1, z);
    dl3 = dx * getEdgeLen(x, y+upper, z, x+1, y+upper, z);
    dl4 = dx * getEdgeLen(x, y+lower, z, x+1, y+lower, z);

    dA = getArea(dl1, dl2, dl3, dl4, dy, dx);
//    if(dA == 0) printf("dA(2) = %g dl1=%g dl2=%g dl3=%g dl4=%g\n", dA, dl1, dl2, dl3, dl4);

    /* dEy/dx  */
    bindx[k] = compY + pos_to_row(x+upper, y, z); val[k++] = +dl1 / dA;
    bindx[k] = compY + pos_to_row(x+lower, y, z); val[k++] = -dl2 / dA;

    /* -dEx/dy  */
    bindx[k] = compX + pos_to_row(x, y+upper, z); val[k++] = -dl3 / dA;
    bindx[k] = compX + pos_to_row(x, y+lower, z); val[k++] = +dl4 / dA;
    break;
  }
  }

  bindx[location+1] = k;  val[location]     = 0.; /* matrix diagonal */
}


//
// write a field to the file. Some primitive serialization is done in
// case of a parallel run.
//
void write_file(char *filename, double *tmp_vec, int N_update){
  FILE* outfile;

#ifdef AZ_MPI
  {
    int commSize;
    int myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    for(int i = 0; i < commSize; i++){
      if(myRank == i){
        if(myRank == 0)
          outfile = fopen(filename, "w");
        else
          outfile = fopen(filename, "a");
        for(int k = 0; k < N_update; k++)
          fprintf(outfile, "%g\n", tmp_vec[k]);

        fclose(outfile);
      }

     MPI_Barrier(MPI_COMM_WORLD);
    }
  }
#else
  outfile = fopen(filename, "w");
  for(int k = 0; k < N_update; k++){
    fprintf(outfile, "%g\n", tmp_vec[k]);
  }
  fclose(outfile);
#endif
}



//
// just some timer for performance menasurement
//
#include <sys/time.h>
#include <unistd.h>

long currentTimeMillis() {

    struct timeval tv ;
    struct timezone tz ;

    gettimeofday(&tv, &tz) ;

    return tv.tv_sec * 1000 + tv.tv_usec / 1000 ;
}
