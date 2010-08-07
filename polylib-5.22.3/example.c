
/* 
  main.c
    This file along with Zpolyhedron.c, polyhedron.c, Lattice.c,  
    Matop.c SolveDio.c, matrix.c and vector.c does the following :

     - Intersection of two Z-Domains.
     - Difference of two Z-domains. 
     - Image of a Z-domain by a invertible, 
        affine rational function. 
*/

#include <polylib/polylib.h>

int main () {
  
  Matrix *a, *b; 
  Polyhedron *P;
  
  a = Matrix_Read ();
  P = Constraints2Polyhedron (a, 200);

  Polyhedron_Print(stdout, P_VALUE_FMT, P);
  
  Matrix_Free (a);

  Domain_Free (P);
  
  return 0;
} /* main */
