/*       Polyhedron disjoint intersections
 */

/*
 disjoint_union_sep computes the disjoint union of the given list of domains.
 input:
    (integer) # of polyhedra
    list of polyhedra in the usual matrix (constraints) format

 output:
    list of polyhedra (constraint matrices) having no integer point in common
*/

#include <stdio.h>
#include <stdlib.h>

#include <polylib/polylib.h>

#define WS 0


/* Procedure to print constraints of a domain */
void AffContraintes(Polyhedron *p)
{
	for( ;p;p=p->next)
	{
		Polyhedron_PrintConstraints(stdout, P_VALUE_FMT, p );
		printf("\n");
	}
}


int main() {
  
	int np, i;

	Matrix *a;
	Polyhedron *A, *tmp, *DD;

	scanf( "%d", &np );

	A = NULL;
	for( i=0 ; i<np ; i++ )
	{
		a = Matrix_Read();
		tmp = Constraints2Polyhedron(a,WS);
		Matrix_Free(a);
		tmp ->next = A;
		A = tmp;
	}


	DD = Disjoint_Domain( A, 0, WS );

	AffContraintes(DD);

	Domain_Free( DD );
	Domain_Free( A );

	return 0;
}



