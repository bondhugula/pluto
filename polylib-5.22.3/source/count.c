/*************************************************/
/*     count.c                                   */
/* program to count the number of points         */
/* in a parameterized polytope                   */
/*                                               */
/* input : polytope                              */
/*         context                               */
/*                                               */
/* written by Vincent Loechner, aug. 2000.       */
/*  loechner@icps.u-strasbg.fr                   */
/*************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <polylib/polylib.h>
#define MAXRAYS 1024

/****************************************************************/
int main(int argc,char *argv[]) {
	
  Matrix *C1, *P1;
  Polyhedron *C, *P, *S;
  Polyhedron *CC, *PP;
  Enumeration *en;
  Value *p;
  int i,j,k;
  int m,M;
  char str[1024];
  Value c;

  /******* Read the input *********/
  P1 = Matrix_Read();
  C1 = Matrix_Read();
  if(C1->NbColumns < 2) {
    fprintf(stderr,"Not enough parameters !\n");
    exit(0);
  }
  P = Constraints2Polyhedron(P1, MAXRAYS);
  C = Constraints2Polyhedron(C1, MAXRAYS);
  Matrix_Free(C1);
  Matrix_Free(P1);
  
  
  /******* Compute the true context *******/
  CC = align_context(C,P->Dimension,MAXRAYS);
  PP = DomainIntersection(P,CC,MAXRAYS);
  Domain_Free(CC);
  C1 = Matrix_Alloc(C->Dimension+1,P->Dimension+1);
  for(i=0;i<C1->NbRows;i++)
    for(j=0;j<C1->NbColumns;j++)
      if(i==j-P->Dimension+C->Dimension)
	value_set_si(C1->p[i][j],1);
      else
	value_set_si(C1->p[i][j],0);
  CC = Polyhedron_Image(PP,C1,MAXRAYS);
  Domain_Free(C);
  Domain_Free(PP);
  Matrix_Free(C1);
  C = CC;
  
  /******* Initialize parameters *********/
  p = (Value *)malloc(sizeof(Value) * (P->Dimension+2));
  for(i=0;i<=P->Dimension;i++) {
    value_init(p[i]);
    value_set_si(p[i],0);
  }
  value_init(p[i]);
  value_set_si(p[i],1);
  
  /*** S = scanning list of polyhedra ***/
  S = Polyhedron_Scan(P,C,MAXRAYS);

  value_init(c);
  
  /******* Count now *********/
  FOREVER {
    fflush(stdin);
    printf("Enter %d parameters : ",C->Dimension);
    for(k=S->Dimension-C->Dimension+1;k<=S->Dimension;++k) {
      scanf(" %s", str);
      value_read(p[k],str);
    }      
    printf("EP( ");
    value_print(stdout,VALUE_FMT,p[S->Dimension-C->Dimension+1]);
    for(k=S->Dimension-C->Dimension+2;k<=S->Dimension;++k) {
      printf(", ");
      value_print(stdout,VALUE_FMT,p[k]);
    }  
    printf(" ) = "); 
    count_points(1,S,p,&c);
    value_print(stdout,VALUE_FMT,c);
    printf("\n"); 
  }
  for(i=0;i<=(P->Dimension+1);i++)
    value_clear(p[i]);
  value_clear(c);
  return(0);
} /* main */



