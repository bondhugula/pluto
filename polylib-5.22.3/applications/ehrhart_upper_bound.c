#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <polylib/polylib.h>


int main( int argc, char **argv)
{
  int i;
  char ** param_name;
  Matrix *C1, *P1;
  Polyhedron *P, *C;
  Enumeration *e, *en;

  Matrix * Validity_Lattice;
  int nb_parms;
  
#ifdef EP_EVALUATION
  Value *p, *tmp;
  int k;
#endif

  P1 = Matrix_Read();
  C1 = Matrix_Read();
  nb_parms = C1->NbColumns-2;
  if(nb_parms < 0) {
    fprintf( stderr, "Not enough parameters !\n" );
    exit(0);
  }

   /* Read the name of the parameters */
  param_name = Read_ParamNames(stdin,nb_parms);

  /* inflate the polyhedron, so that the inflated EP approximation will be an
     upper bound for the EP's polyhedron. */
  mpolyhedron_inflate(P1,nb_parms);

  /* compute a polynomial approximation of the Ehrhart polynomial */
  e = Ehrhart_Quick_Apx(P1, C1, &Validity_Lattice, 1024);
  
  Matrix_Free(C1);
  Matrix_Free(P1);

  printf("============ Ehrhart polynomial quick polynomial upper bound ============\n");
  show_matrix(Validity_Lattice);
  for( en=e ; en ; en=en->next ) {    
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
