#include <stdio.h>
#include <stdlib.h>

#include <polylib/polylib.h>

int main() {

  Matrix *a, *b;
  Polyhedron *A, *B;
  Param_Polyhedron *PA;
  char **param_name;
  
  a = Matrix_Read();
  A = Constraints2Polyhedron(a,200);
  Matrix_Free(a);
  
  b = Matrix_Read();
  B = Constraints2Polyhedron(b,200);
  Matrix_Free(b);
  
  /* Read the name of the parameters */
  param_name = Read_ParamNames(stdin,B->Dimension);  
  PA = Polyhedron2Param_Vertices(A,B,500);  
  Param_Vertices_Print(stdout,PA->V,param_name);
  Domain_Free(A);
  Domain_Free(B);
  Param_Polyhedron_Free( PA );
  free(param_name);
  return 0;
} /* main */



