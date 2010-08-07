
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
  ZPolyhedron *Z1, *Z2, *Z3, *Z4;
  
  a = Matrix_Read ();
  b = Matrix_Read ();
  P = Constraints2Polyhedron (b, 200);
  Z1 = ZPolyhedron_Alloc (a, P);
  
  Matrix_Free (a);
  Matrix_Free (b);
  Domain_Free (P);
  
  a = Matrix_Read ();
  b = Matrix_Read ();
  P = Constraints2Polyhedron (b, 200);
  Z2 = ZPolyhedron_Alloc (a, P);
  
  Matrix_Free (a); 
  Matrix_Free (b); 
  Domain_Free (P);
  
  Z3 = ZDomainIntersection (Z1, Z2);
  printf ("\nZ3 = Z1 and Z2");
  ZDomainPrint(stdout,P_VALUE_FMT, Z3);

  a = Matrix_Read ();
  Z4 = ZDomainImage (Z1, a);
  printf ("\nZ4 = image (Z1 by a)");
  ZDomainPrint (stdout,P_VALUE_FMT, Z4);

  Matrix_Free (a);
  ZDomain_Free (Z1);
  ZDomain_Free (Z2);
  ZDomain_Free (Z3);
  ZDomain_Free (Z4);
  
  return 0;
} /* main */
