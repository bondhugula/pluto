/* polytest.c */
#include <stdio.h>
#include <polylib/polylib.h>

#define WS 0

char s[128];

int main() { 
  
  Matrix *a=NULL, *b=NULL, *c, *d, *e, *f;
  Polyhedron *A, *B, *C, *D, *last, *tmp;
  int i, nbPol, nbMat, func;
  
  fgets(s, 128, stdin);
  nbPol = nbMat = 0;
  while ((*s=='#') ||
	  ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
    fgets(s, 128, stdin);

  for (i=0, A=last=(Polyhedron *)0; i<nbPol; i++) {
    a = Matrix_Read();
    tmp = Constraints2Polyhedron(a,WS);
    Matrix_Free(a);
    if (!last) A = last = tmp;
    else {
      last->next = tmp;
      last = tmp;
    }
    }

    if (nbMat)
    {  a = Matrix_Read(); }

    fgets(s,128,stdin);
    nbPol = nbMat = 0;
    while ( (*s=='#') ||
        ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
      fgets(s, 128, stdin);

    for (i=0, B=last=(Polyhedron *)0; i<nbPol; i++) {
      b = Matrix_Read();
      tmp = Constraints2Polyhedron(b,WS);
      Matrix_Free(b);
      if (!last) B = last = tmp;
      else {
	last->next = tmp;
	last = tmp;
      }
    }
    
    if (nbMat)
      {  b = Matrix_Read(); }
    
    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "F %d", &func)<1))
      fgets(s, 128, stdin);
    
    switch (func) {
    case 1:
      C = DomainUnion(A, B, WS);
      D = DomainConvex(C, WS);
      d = Polyhedron2Constraints(D);
      Matrix_Print(stdout,P_VALUE_FMT,d);
      Matrix_Free(d);
      Domain_Free(C);
      Domain_Free(D);
      break;
    case 2:
      D = DomainSimplify(A, B, WS);
      d = Polyhedron2Constraints(D);
      Matrix_Print(stdout,P_VALUE_FMT,d);
      Matrix_Free(d);
      Domain_Free(D);
      break;
    case 3:
      a = Polyhedron2Constraints(A);
      Matrix_Print(stdout,P_VALUE_FMT,a);
      b = Polyhedron2Constraints(B);
      Matrix_Print(stdout,P_VALUE_FMT,b);
      break;
    case 4:
      a = Polyhedron2Rays(A);
      Matrix_Print(stdout,P_VALUE_FMT,a);
      break;
    case 5:
      
      /* a = ec , da = c , ed = 1 */
      right_hermite(a,&c,&d,&e);
      Matrix_Print(stdout,P_VALUE_FMT,c);
      Matrix_Print(stdout,P_VALUE_FMT,d);
      Matrix_Print(stdout,P_VALUE_FMT,e);
      f = Matrix_Alloc(e->NbRows,c->NbColumns);
      Matrix_Product(e,c,f);
      Matrix_Print(stdout,P_VALUE_FMT,f);
      Matrix_Free(f);
      f = Matrix_Alloc(d->NbRows,a->NbColumns);
      Matrix_Product(d,a,f);
      Matrix_Print(stdout,P_VALUE_FMT,f);
      Matrix_Free(f);
      f = Matrix_Alloc(e->NbRows, d->NbColumns);
      Matrix_Product(e,d,f);
      Matrix_Print(stdout,P_VALUE_FMT,f);
      break;
    case 6:
      
      /* a = ce , ad = c , de = 1 */
      left_hermite(a,&c,&d,&e);
      Matrix_Print(stdout,P_VALUE_FMT,c);
      Matrix_Print(stdout,P_VALUE_FMT,d);
      Matrix_Print(stdout,P_VALUE_FMT,e);
      f = Matrix_Alloc(c->NbRows, e->NbColumns);
      Matrix_Product(c,e,f);
      Matrix_Print(stdout,P_VALUE_FMT,f);
      Matrix_Free(f);
      f = Matrix_Alloc(a->NbRows, d->NbColumns);
      Matrix_Product(a,d,f);
      Matrix_Print(stdout,P_VALUE_FMT,f);
      Matrix_Free(f);
      f = Matrix_Alloc(d->NbRows, e->NbColumns);
      Matrix_Product(d,e,f);
      Matrix_Print(stdout,P_VALUE_FMT,f);
      break;
    case 7:	          
     
      /* Polyhedron_Print(stdout,"%5d", A); */
      /* Matrix_Print(stdout,"%4d", b);     */
      
      C = Polyhedron_Image(A, b, WS);
      Polyhedron_Print(stdout,P_VALUE_FMT,C);
      break;
    case 8:
      
      printf("%s\n",
	     Polyhedron_Not_Empty(A,B,WS) ? "Not Empty" : "Empty");
      break;
    case 9:
      
      i = PolyhedronLTQ(A,B,1,0,WS);
      printf("%s\n",
	     i==-1 ? "A<B" : i==1 ? "A>B" : i==0 ? "A><B" : "error");
      i = PolyhedronLTQ(B,A,1,0,WS);
      printf("%s\n",
	     i==-1 ? "A<B" : i==1 ? "A>B" : i==0 ? "A><B" : "error");
      break;
    case 10:
      i = GaussSimplify(a,b);
      Matrix_Print(stdout,P_VALUE_FMT,b);
      break;
     
    default:
      printf("? unknown function\n");
    }
    
    Domain_Free(A);
    Domain_Free(B);
    
    return 0;
}


