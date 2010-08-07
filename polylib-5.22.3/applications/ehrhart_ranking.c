/* Program for testing the ranking function. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <polylib/polylib.h>
#include <polylib/ranking.h>

int main( int argc, char **argv)
{
  int i;
  char ** param_name = NULL;
  Matrix *M;
  Polyhedron *P, *D, *C;
  Enumeration *e, *en;

  int nb_parms;
  
#ifdef EP_EVALUATION
  Value *p, *tmp;
  int k;
#endif

  M = Matrix_Read();
  P = Constraints2Polyhedron(M, POL_NO_DUAL);
  Matrix_Free(M);
  M = Matrix_Read();
  D = Constraints2Polyhedron(M, POL_NO_DUAL);
  Matrix_Free(M);
  M = Matrix_Read();
  C = Constraints2Polyhedron(M, POL_NO_DUAL);
  Matrix_Free(M);

  nb_parms = D->Dimension;

   /* Read the name of the parameters */
  param_name = Read_ParamNames(stdin,nb_parms);

  /* compute a polynomial approximation of the Ehrhart polynomial */
  printf("============ Ranking function ============\n");
  e = Polyhedron_LexSmallerEnumerate(P, D, D->Dimension-C->Dimension, 
				     C, POL_NO_DUAL);
  
  Polyhedron_Free(P);
  Polyhedron_Free(D);
  Polyhedron_Free(C);

  for (en=e; en; en=en->next) {
    Print_Domain(stdout,en->ValidityDomain, param_name);
    print_evalue(stdout,&en->EP, param_name);
    printf( "\n-----------------------------------\n" );
  }

 
#ifdef EP_EVALUATION
  if( isatty(0) && nb_parms != 0)
  {  /* no tty input or no polyhedron -> no evaluation. */
    printf("Evaluation of the Ehrhart polynomial :\n");
    p = (Value *)malloc(sizeof(Value) * (nb_parms));
    for(i=0;i<nb_parms;i++) 
      value_init(p[i]);
    FOREVER {
      fflush(stdin);
      printf("Enter %d parameters : ",nb_parms);
      for(k=0;k<nb_parms;++k) {
	scanf("%s",str);
	value_read(p[k],str);
      }
      fprintf(stdout,"EP( ");
      value_print(stdout,VALUE_FMT,p[0]);
      for(k=1;k<nb_parms;++k) {
	fprintf(stdout,",");
	value_print(stdout,VALUE_FMT,p[k]);
      }
      fprintf(stdout," ) = ");
      value_print(stdout,VALUE_FMT,*(tmp=compute_poly(en,p)));
      free(tmp);
      fprintf(stdout,"\n");  
    }
  }
#endif /* EP_EVALUATION */
  
  while( e )
    {
      free_evalue_refs( &(e->EP) );
      Polyhedron_Free( e->ValidityDomain );
      en = e ->next;
      free( e );
      e = en;
    }
  for( i=0 ; i<nb_parms ; i++ )
    free( param_name[i] );
  free(param_name);
  return 0;
}
