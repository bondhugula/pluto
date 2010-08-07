#define WS 0

#include <stdlib.h>
#include <polylib/polylib.h>

int main() {
	
  Matrix *A;
  Polyhedron *P;
  
  A = Matrix_Read();
  if(A->NbColumns < 2) {
    printf("Wrong input: %d columns\n", A->NbColumns );
    Matrix_Free(A);
    exit(1);
  }
  Matrix_Print(stdout,P_VALUE_FMT,A);
  P = Constraints2Polyhedron(A,WS);
  Matrix_Free(A);  
  Polyhedron_Print(stdout,P_VALUE_FMT,P);
  Domain_Free(P);
  return 0;  
}
