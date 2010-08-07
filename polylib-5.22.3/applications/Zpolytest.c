/* zpolytest.c
This is a testbench for the Zpolylib (part of polylib manipulating 
Z-polyhedra. */

#include <stdio.h>
#include <polylib/polylib.h>

#define WS 0

char s[128];

int main() {
  
  Matrix *a=NULL, *b=NULL, *c=NULL, *d, *e, *g;
  LatticeUnion *l1,*l2,*l3,*l4,*temp;
  Polyhedron *A=NULL, *B=NULL, *C=NULL, *D;
  ZPolyhedron *ZA, *ZB, *ZC, *ZD, *Zlast;
  int  nbPol, nbMat, func, rank ;
  Vector *v=NULL;
    
  /* The structure of the input file to this program  is the following: 
       First a line containing        
           M nbMat
       Where nbMat is an integer indicating how many Matrices will 
       be described in the following. temporary debugging. Next 
       the matrice are described. For each matrix, the first row is two
       integers:
           nbRows nbColumns 
       Then the matrix is written row by row. a line starting with 
      a `#' is considered as a comment 

      Then a line containing 
           D nbDomain
      where nbDomain is an integer indicating how many domain will 
       be described in the following. Domains are describled as for 
       polylib,  the first row is two integers:
          nbConstraints dimension
      then the constraints are described in the Polylib format.
      The last line of the input file contains :
          F numTest
          which indicates which test will be performed on the data.
      Warning, currently no more than 3 matrice of Polyhedra can be read*/
  
  fgets(s, 128, stdin);
  nbPol = nbMat = 0;
  while ( (*s=='#') ||
	  ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
    fgets(s, 128, stdin);
  
  
  /* debug */
  /*     fprintf(stdout,"nbMat=%d",nbMat);fflush(stdout); */
  
  switch (nbMat) {
    
  case 1: 
    a = Matrix_Read();
    break;
  
  case 2: 
    a = Matrix_Read();
    b = Matrix_Read();
    break;
  
  case 3: a = Matrix_Read();
    b = Matrix_Read();
    c = Matrix_Read();
    break;
  }
  
  fgets(s, 128, stdin);
  while ((*s=='#') ||
	 ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
    fgets(s, 128, stdin);
  
  
  /* debug */
  /*  fprintf(stdout,"nbPol=%d",nbPol);fflush(stdout);  */
  
  switch (nbPol) { 
  
  case 1:  
    g = Matrix_Read();
    A = Constraints2Polyhedron(g,WS);
    Matrix_Free(g);
    break;
  
  case 2:         
    g = Matrix_Read();
    A = Constraints2Polyhedron(g,WS);
    Matrix_Free(g);
    g = Matrix_Read();
    B = Constraints2Polyhedron(g,WS);
    Matrix_Free(g);
    break;
  
  case 3:
    g = Matrix_Read();
    A = Constraints2Polyhedron(g,WS);
    Matrix_Free(g);
    g = Matrix_Read();
    B = Constraints2Polyhedron(g,WS);
    Matrix_Free(g);
    g = Matrix_Read();
    C = Constraints2Polyhedron(g,WS);
    Matrix_Free(g);
    break;
  }
  
  fgets(s, 128, stdin);
  while ((*s=='#') || (sscanf(s, "F %d", &func)<1) ) fgets(s, 128, stdin);
  

  switch (func) {
    
  case 1:
    
    /* just a  test of polylib functions */
    C = DomainUnion(A, B, 200);
    D = DomainConvex(C, 200);
    d = Polyhedron2Constraints(D);
    Matrix_Print(stdout,P_VALUE_FMT, d);
    Matrix_Free(d);
    Domain_Free(D);
    break;
    
  case 2: /* AffineHermite */
    
    AffineHermite(a,&b,&c);
    Matrix_Print(stdout,P_VALUE_FMT, b);
    Matrix_Print(stdout,P_VALUE_FMT, c);
    break;
    
  case 3: /* LatticeIntersection */
    
    c = LatticeIntersection(a,b);
    Matrix_Print(stdout,P_VALUE_FMT, c);
    break;
    
  case 4: /* LatticeDifference */
        
    fprintf(stdout," 2 in 1 : %d\n",LatticeIncludes(b,a));
    fprintf(stdout," 1 in 3 : %d\n",LatticeIncludes(c,a));
    fprintf(stdout," 1 in 2 : %d\n",LatticeIncludes(a,b));
    break;
  
  case 5: /* LatticeDifference */
    
    l1=LatticeDifference(a,b);
    l2=LatticeDifference(b,a);
    l3=LatticeDifference(c,a);
    l4=LatticeDifference(b,c);
    fprintf(stdout,"L1 - L2 :\n");
    temp=l1;
    while (temp!=NULL) {
      
      Matrix_Print(stdout,P_VALUE_FMT,temp->M);
      temp=temp->next; 
    };
    fprintf(stdout,"Diff2:\n");
    temp=l2;
    while (temp!=NULL) {
      Matrix_Print(stdout,P_VALUE_FMT, temp->M);
      temp=temp->next; 
    };
    fprintf(stdout,"Diff3:\n");
    temp=l3;
    while (temp!=NULL) {
      Matrix_Print(stdout,P_VALUE_FMT, temp->M);
      temp=temp->next; 
    };
    fprintf(stdout,"Diff4:\n");
    temp=l4;
    while (temp!=NULL) {
      Matrix_Print(stdout,P_VALUE_FMT, temp->M);
      temp=temp->next; 
    };
    break;
    
  case 6: /* isEmptyZPolyhedron */
    
    ZA=ZPolyhedron_Alloc(a,A);
    fprintf(stdout,"is Empty? :%d \n", isEmptyZPolyhedron(ZA));
    ZDomain_Free(ZA);
    break;
    
  case 7: /* ZDomainIntersection */
        
    ZA=ZPolyhedron_Alloc(a,A);
    ZB=ZPolyhedron_Alloc(b,B);
    ZC = ZDomainIntersection(ZA,ZB);
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    ZDomain_Free(ZA);
    ZDomain_Free(ZB);
    ZDomain_Free(ZC);
    break;
    
  case 8: /* ZDomainUnion */
    
    ZA=ZPolyhedron_Alloc(a,A);
    ZB=ZPolyhedron_Alloc(b,B);
    ZC = ZDomainUnion(ZA,ZB);
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    break;
    
  case 9: /* ZDomainDifference */
    
    ZA=ZPolyhedron_Alloc(a,A);
    ZB=ZPolyhedron_Alloc(b,B);
    ZC = ZDomainDifference(ZA,ZB);
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    break;
    
  case 10: /* ZDomainImage */
    
    ZA=ZPolyhedron_Alloc(a,A);
    ZC = ZDomainImage(ZA,b); 
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    break;
    
  case 11: /* ZDomainPreimage */
    
    ZA=ZPolyhedron_Alloc(a,A);
    ZC = ZDomainPreimage(ZA,b); 
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    break;
    
  case 12: /* ZDomainDifference */
    ZA=ZPolyhedron_Alloc(a,A);
    ZC = ZDomainPreimage(ZA,b); 
    ZD = ZDomainImage(ZC,b); 
    Zlast=ZDomainDifference(ZD,ZC);
    fprintf(stdout,"the Two zpol are equal? :%d\n",
	    isEmptyZPolyhedron(Zlast));
    break;
  
  case 13:  /* ZDomainSimplify */
    
    ZA=ZPolyhedron_Alloc(a,A);
    ZA->next = ZPolyhedron_Alloc(b,B);
    ZDomainPrint(stdout,P_VALUE_FMT, ZA);
    ZD = ZDomainSimplify(ZA);
    ZDomainPrint(stdout,P_VALUE_FMT, ZD);
    break;
    
  case 14:  /* EmptyZpolyhedron */
        
    ZA=EmptyZPolyhedron(3);
    fprintf(stdout,"is Empty? :%d \n", isEmptyZPolyhedron(ZA));
    ZDomain_Free(ZA);
    break;
    
  case 15:  /* ZDomainInclude */
  
    ZA=ZPolyhedron_Alloc(a,A);
    ZB=ZPolyhedron_Alloc(b,B);
    fprintf(stdout,"A in B  :%d \nB in A  :%d \n", 
	    ZPolyhedronIncludes(ZA,ZB),
	    ZPolyhedronIncludes(ZB,ZA));
    break;
  
  case 16: /* LatticePreimage */
        
    c = LatticePreimage(a,b);
    Matrix_Print(stdout,P_VALUE_FMT, c);
    AffineHermite(c,&d,&e);
    Matrix_Print(stdout,P_VALUE_FMT, d);
    break;
    
  case 17: /* LatticeImage */
    
    c = LatticeImage(a,b);
    Matrix_Print(stdout,P_VALUE_FMT, c);
    AffineHermite(c,&d,&e);
    Matrix_Print(stdout,P_VALUE_FMT, d);
    break;
      
  case 18:  /* EmptyLattice */
        
    fprintf(stdout,"is Empty? :%d \n", isEmptyLattice(a));
    fprintf(stdout,"is Empty? :%d \n", isEmptyLattice(EmptyLattice(3)));
    break;
  
  case 19:  /* CanonicalForm */
     
    ZA=ZPolyhedron_Alloc(a,A);
    ZB=ZPolyhedron_Alloc(a,B);
    CanonicalForm(ZA,&ZC,&c);
    CanonicalForm(ZB,&ZD,&d);
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    ZDomainPrint(stdout,P_VALUE_FMT, ZD);
    break;
    
  case 20: /* LatticeSimplify */
    
    l1=LatticeUnion_Alloc();
    l2=LatticeUnion_Alloc();
    l1->M=Matrix_Copy(a);
    l1->next=l2;
    l2->M=Matrix_Copy(b);
    l1=LatticeSimplify(l1);
    PrintLatticeUnion(stdout,P_VALUE_FMT,l1); 
    LatticeUnion_Free(l1);
    break;
    
  case 21: /* AffineSmith */
  
    AffineSmith(a,&b,&c, &d);
    Matrix_Print(stdout,P_VALUE_FMT, b); 
    Matrix_Print(stdout,P_VALUE_FMT, c);
    Matrix_Print(stdout,P_VALUE_FMT, d);
    Matrix_Free(d);
    break;
  
  case 22: /* SolveDiophantine */
        
    rank=SolveDiophantine(a,&d,&v); 
    Matrix_Print(stdout,P_VALUE_FMT, a);
    fprintf(stdout," rank: %d \n ",rank);
    Matrix_Print(stdout,P_VALUE_FMT, d); 
    Vector_Print(stdout,P_VALUE_FMT, v);
    rank=SolveDiophantine(b,&d,&v); 
    Matrix_Print(stdout,P_VALUE_FMT, b);
    fprintf(stdout," rank: %d \n ",rank);
    Matrix_Print(stdout,P_VALUE_FMT, d); 
    Vector_Print(stdout,P_VALUE_FMT, v);
    rank=SolveDiophantine(c,&d,&v); 
    Matrix_Print(stdout,P_VALUE_FMT, c);
    fprintf(stdout," rank: %d \n ",rank);
    Matrix_Print(stdout,P_VALUE_FMT, d); 
    Vector_Print(stdout,P_VALUE_FMT, v);
    Vector_Free(v);
    break;

  case 23: /* SplitZPolyhedron */
        
    ZA=ZPolyhedron_Alloc(a,A);
    ZC = SplitZpolyhedron(ZA,b);
    ZDomainPrint(stdout,P_VALUE_FMT, ZC);
    break;


  case 100: /* debug */
    
    ZA=ZPolyhedron_Alloc(a,A);
    ZDomainPrint(stdout,P_VALUE_FMT, ZA);
    ZDomain_Free(ZA);
    break;
    
  default:
    printf("? unknown function\n");
  }

  /*    Polyhedron_Free(A); */
  if (a)
    Matrix_Free(a);
  if (b)
    Matrix_Free(b);
  if (c)
    Matrix_Free(c);

  if (A)
    Domain_Free(A);
  if (B)
    Domain_Free(B);
  if (C)
    Domain_Free(C);
  
  return 0;
} /* main */
