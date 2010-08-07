#include <stdio.h>
#include <stdlib.h>

#include <polylib/polylib.h>

#define WS 0

int main() {
	
  Matrix *a, *b;
  Polyhedron *A, *B;
  Param_Polyhedron *PA;
  Param_Domain *P;
  Param_Vertices *V;
  int nbPV;
  char **param_name;
  
  a = Matrix_Read();
  if(!a || a->NbColumns == 0) {
    fprintf(stderr,"Input error: empty matrix\n");
    exit(0);
  }
  A = Constraints2Polyhedron(a, WS);
  Matrix_Free(a);
  b = Matrix_Read();
  
  if(!b || b->NbColumns == 0) {
    fprintf(stderr, "Input error: empty matrix\n");
    exit(0);
  }
  B = Constraints2Polyhedron(b, WS);
  Matrix_Free(b);
  
  /* Read the name of the parameters */
  param_name = Read_ParamNames(stdin, B->Dimension);
  PA = Polyhedron2Param_Domain(A,B,WS);
  if(!PA || PA->D==NULL) {
    printf("---------------------------------------\n");
    printf("Empty polyhedron\n");
    return 0;
  }
  nbPV = PA->nbV;
  Domain_Free(A);
  Domain_Free(B);

  /*****************************/
  /* Scan the validity domains */
  for(P=PA->D;P;P=P->next) {
    
    /* prints current val. dom. */
    printf( "---------------------------------------\n" );
    printf( "Domain :\n");
    Print_Domain( stdout, P->Domain, param_name );
    
    /* scan the vertices */
    printf( "Vertices :\n");
    FORALL_PVertex_in_ParamPolyhedron(V,P,PA) {
	
      /* prints each vertex */
      Print_Vertex( stdout, V->Vertex, param_name );
      printf( "\n" );
    }
    END_FORALL_PVertex_in_ParamPolyhedron;
  }
  /*****************************/
  
  Param_Polyhedron_Free( PA );
  free(param_name);
  
  return 0;
} /* main */ 

