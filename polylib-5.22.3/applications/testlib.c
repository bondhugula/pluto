/* main.c
     COPYRIGHT
          Both this software and its documentation are

              Copyright 1993 by IRISA /Universite de Rennes I - France,
              Copyright 1995,1996 by BYU
                         all rights reserved.

          Permission is granted to copy, use, and distribute
          for any commercial or noncommercial purpose under the terms
          of the GNU General Public license, version 2, June 1991
          (see file : LICENSING).

   This file along with polyhedron.c and vector.c do the following functions:
    - Extraction of a minimal set of constraints from some set of constraints
    - Intersection of two convexes
    - Application of a linear function to some convex
    - Verification that a convex is included in some other convex

   They are compiled together into an executable file called "test".
   The file test.in contains sample input data for the program.
   The file test.out contains the output produced by the program.

   This directory also contains a makefile to build and run the test.

   This file is a tutorial on how to use the library.  The comments
   explain whats going on.  You can use this file as a pattern for your
   own program.  Good Luck !

   --Doran
*/

#include <stdio.h>

#include <polylib/polylib.h>


int main() { 
  
  Matrix *a, *b, *t;
  Polyhedron *A, *B, *C, *D;
  
  printf("Polyhedral Library Test\n\n");
  
  /* read in a matrix containing your equations */
  /* for example, run this program and type in these five  lines:
     4 4
     1 0 1 -1
     1 -1 0 6
     1 0 -1 7
     1 1 0 -2
     This is a matrix for the following inequalities
     1 = inequality,  0x +  1y -1 >=0  -->	y >= 1
     1 = inequality, -1x +  0y +6 >=0  -->	x <= 6
     1 = inequality,  0x + -1y +7 >=0  -->	y <= 7
     1 = inequality,  1x +  0y -2 >=0  -->	x >= 2
     If the first number is a 0 instead of a 1, then that constraint
     is an 'equality' instead of an 'inequality'.
  */
  a = Matrix_Read();
  
  /* read in a second matrix containing a second set of constraints:
     for example :
     4 4
     1 1 0 -1
     1 -1 0 3
     1 0 -1 5
     1 0 1 -2
  */
  b = Matrix_Read();

  /* Convert the constraints to a Polyhedron.
     This operation 1. Computes the dual ray/vertice form of the
     system, and 2. Eliminates redundant constraints and reduces
     them to a minimal form.
  */
  A = Constraints2Polyhedron(a, 200);
  B = Constraints2Polyhedron(b, 200);

  /* the 200 is the size of the working space (in terms of number
     of rays) that is allocated temporarily
     -- you can enlarge or reduce it as needed */
  
  /* There is likewise a rays to polyhedron procedure */
  
  /* Since you are done with the matrices a and b, be a good citizen
     and clean up your garbage */
  Matrix_Free(a);
  Matrix_Free(b);
  
  /* If you want the the reduced set of equations back, you can
     either learn to read the polyhedron structure (not hard,
     look in "types.h"...
     or you can get the matrix back in the same format it started
     in... */
  a = Polyhedron2Constraints(A);
  b = Polyhedron2Constraints(B);

  /* Take a look at them if you want */
  printf("\na =");
  Matrix_Print(stdout,P_VALUE_FMT,a);
  printf("\nb =");
  Matrix_Print(stdout,P_VALUE_FMT,b);

  Matrix_Free(a);
  Matrix_Free(b);
  
  /* To intersect the two systems, use the polyhedron formats... */
  C = DomainIntersection(A, B, 200);
  
  /* Again, the 200 is the size of the working space */
  
  /* This time, lets look a the polyhedron itself... */
  printf("\nC = A and B =");
  Polyhedron_Print(stdout,P_VALUE_FMT,C);
  
  /* 
   * The operations DomainUnion, DomainDifference, and DomainConvex
   * are also available 
   */
  
  /* 
   * Read in a third matrix containing a transformation matrix,
   * this one swaps the indices (x,y --> y,x):
   * 3 3
   * 0 1 0
   * 1 0 0
   * 0 0 1
   */
  
  
  t = Matrix_Read();
  
  /* Take the preimage (transform the equations) of the domain C to
     get D.  Are you catching on to this 200 thing yet ??? */
  
  D = Polyhedron_Preimage(C, t, 200);
  
  /* cleanup please */
  Matrix_Free(t);
  
  printf("\nD = transformed C =");
  Polyhedron_Print(stdout,P_VALUE_FMT,D);
  Domain_Free(D);
  
  /* in a similar way, Polyhedron_Image(dom, mat, 200), takes the image
     of dom under matrix mat  (transforms the vertices/rays) */
  
  /* The function PolyhedronIncludes(Pol1, Pol2) returns 1 if Pol1
     includes (covers) Pol2, and a 0 otherwise */
  
  if (PolyhedronIncludes(A,C))
    printf("\nWe expected A to cover C since C = A intersect B\n");
  if (!PolyhedronIncludes(C,B))
    printf("and C does not cover B...\n");
  
  /* Final note:  some functions are defined for Domains, others
   * for Polyhedrons.  A domain is simply a list of polyhedrons.
   * Every polyhedron structure has a "next" pointer used to 
   * make a list of polyhedrons...  For instance, the union of
   * two disjoint polyhedra is a domain consisting of two polyhedra.
   * If you want the convex domain... you have to call 
   * DomainConvex(Pol1, 200) explicity.   
   * Note that includes does not work on domains, only on simple
   * polyhedrons...
   * Thats the end of the demo...  write me if you have questions.
   * And remember to clean up... 
   */
  
  Domain_Free(A);
  Domain_Free(B);
  Domain_Free(C);
  
  return 0;
}


