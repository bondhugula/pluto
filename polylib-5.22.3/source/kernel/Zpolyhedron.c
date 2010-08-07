#include <polylib/polylib.h> 
#include <stdlib.h>

static ZPolyhedron * ZPolyhedronIntersection(ZPolyhedron *, ZPolyhedron *);
static ZPolyhedron *ZPolyhedron_Copy(ZPolyhedron *A);
static void ZPolyhedron_Free(ZPolyhedron *Zpol);
static ZPolyhedron * ZPolyhedronDifference(ZPolyhedron *, ZPolyhedron *);
static ZPolyhedron * ZPolyhedronImage(ZPolyhedron *, Matrix *);
static ZPolyhedron * ZPolyhedronPreimage(ZPolyhedron *, Matrix *);
static ZPolyhedron *AddZPolytoZDomain(ZPolyhedron *A, ZPolyhedron *Head);
static void ZPolyhedronPrint(FILE *fp, char *format, ZPolyhedron *A);

typedef struct forsimplify {
  Polyhedron *Pol;
  LatticeUnion *LatUni;
  struct forsimplify *next;
} ForSimplify;


/* 
 * Returns True if 'Zpol' is empty, otherwise returns False 
 */
Bool isEmptyZPolyhedron (ZPolyhedron *Zpol) {
  
  if(Zpol == NULL)
    return True;
  if((isEmptyLattice (Zpol->Lat)) || (emptyQ(Zpol->P)))
    return True;
  return False;
} /* isEmptyZPolyhedron */

/*
 * Given Lattice 'Lat' and a Polyhderon 'Poly', allocate space, and return
 * the Z-polyhderon corresponding to the image of the polyhderon 'Poly' by the 
 * lattice 'Lat'. If the input lattice 'Lat' is not integeral, it integralises 
 * it, i.e. the lattice of the Z-polyhderon returned is integeral. 
 */ 
ZPolyhedron *ZPolyhedron_Alloc(Lattice *Lat, Polyhedron *Poly) {
  
  ZPolyhedron *A;
  
  POL_ENSURE_FACETS(Poly);
  POL_ENSURE_VERTICES(Poly);

  if(Lat->NbRows != Poly->Dimension+1) {
    fprintf(stderr,"\nInZPolyAlloc - The Lattice  and the Polyhedron");
    fprintf(stderr," are not compatible to form a ZPolyhedra\n");
    return NULL;
  }  
  if((!(isEmptyLattice(Lat))) && (!isfulldim (Lat))) {
    fprintf(stderr,"\nZPolAlloc: Lattice not Full Dimensional\n");
    return NULL;
  }
  A = (ZPolyhedron *)malloc(sizeof(ZPolyhedron));
  if (!A) {
    fprintf(stderr,"ZPolAlloc : Out of Memory\n");
    return NULL;
  } 
  A->next = NULL;
  A->P = Domain_Copy(Poly);
  A->Lat = Matrix_Copy(Lat);
  
  if(IsLattice(Lat) == False) {
    ZPolyhedron *Res;
    
    Res = IntegraliseLattice (A);
    ZPolyhedron_Free (A);
    return Res;
  }
  return A;
} /* ZPolyhedron_Alloc */

/*
 * Free the memory used by the Z-domain 'Head' 
 */
void ZDomain_Free (ZPolyhedron *Head) {
  
  if (Head == NULL)
    return;
  if (Head->next != NULL)
    ZDomain_Free(Head->next);  
  ZPolyhedron_Free(Head);
} /* ZDomain_Free */

/*
 * Free the memory used by the Z-polyhderon 'Zpol' 
 */
static void ZPolyhedron_Free (ZPolyhedron *Zpol) {
  
  if (Zpol == NULL)
    return;
  Matrix_Free((Matrix *) Zpol->Lat);
  Domain_Free(Zpol->P);
  free(Zpol);
  return;
} /* ZPolyhderon_Free */

/*
 * Return a copy of the Z-domain 'Head' 
 */ 
ZPolyhedron *ZDomain_Copy(ZPolyhedron *Head) {
  
  ZPolyhedron *Zpol;
  Zpol = ZPolyhedron_Copy(Head);
  
  if (Head->next != NULL)
    Zpol->next = ZDomain_Copy(Head->next);
  return Zpol;
} /* ZDomain_Copy */

/*
 * Return a copy of the Z-polyhderon 'A' 
 */
static ZPolyhedron *ZPolyhedron_Copy(ZPolyhedron *A) {

  ZPolyhedron *Zpol;
  
  Zpol = ZPolyhedron_Alloc(A->Lat, A->P);
  return Zpol;
} /* ZPolyhderon_Copy */

/* 
 * Add the ZPolyhedron 'Zpol' to the Z-domain 'Result' and return a pointer 
 * to the new Z-domain. 
 */
static ZPolyhedron *AddZPoly2ZDomain(ZPolyhedron *Zpol, ZPolyhedron *Result) {
 
  ZPolyhedron *A;
  
  if (isEmptyZPolyhedron(Zpol))
    return Result;
  A = ZPolyhedron_Copy(Zpol);
  A->next = NULL;
  
  if (isEmptyZPolyhedron (Result)) {
    ZDomain_Free (Result);
    return A;
  }
  A->next = Result;  
  return A;
} /* AddZPoly2ZDomain */
  
/*
 * Given a Z-polyhderon 'A' and a Z-domain 'Head', return a new Z-domain with 
 * 'A' added to it. If the new Z-polyhedron 'A', is already included in the 
 * Z-domain 'Head', it is not added in the list. Othewise, the function checks 
 * if the new Z-polyhedron 'A' to be added to the Z-domain 'Head' has a common
 * lattice with some other Z-polyhderon already present in the Z-domain. If it 
 * is so, it takes the union of the underlying polyhdera; domains and returns. 
 * The function tries to make sure that the added Z-polyhedron 'A' is in the 
 * canonical form.
 */
static ZPolyhedron *AddZPolytoZDomain(ZPolyhedron *A, ZPolyhedron *Head) {
  
  ZPolyhedron *Zpol, *temp, *temp1;
  Polyhedron *i;
  Bool Added;  
  
  if ((A == NULL) || (isEmptyZPolyhedron(A)))
    return Head;
  
  /* For each "underlying" Pol, find the Cnf and add Zpol in Cnf*/  
  for(i=A->P; i!= NULL; i=i->next) {
    ZPolyhedron *Z, *Z1;
    Polyhedron *Image;
    Matrix *H, *U;
    Lattice *Lat ;
    
    Added = False;    
    Image = Domain_Copy(i);
    Domain_Free(Image->next);
    Image->next = NULL;
    Z1 = ZPolyhedron_Alloc(A->Lat,Image);
    Domain_Free(Image);
    CanonicalForm(Z1,&Z,&H); 
    ZDomain_Free(Z1);
    Lat = (Lattice *)Matrix_Alloc(H->NbRows,Z->Lat->NbColumns);
    Matrix_Product(H,Z->Lat,(Matrix *)Lat);
    Matrix_Free(H);    
    AffineHermite(Lat,(Lattice **)&H,&U);
    Image = DomainImage(Z->P,U,MAXNOOFRAYS);
    ZDomain_Free(Z);
      
    Zpol=ZPolyhedron_Alloc((Lattice *)H,Image);     
    Domain_Free(Image);
    Matrix_Free((Matrix *)Lat);
    Matrix_Free(H);
    Matrix_Free(U);
    
    if ((Head == NULL) || (isEmptyZPolyhedron (Head))) {
      Head = Zpol;
      continue;
    }     
    temp1 = temp = Head;
    
    /* Check if the curr pol is included in the zpol or vice versa. */    
    for(; temp != NULL; temp = temp->next) {
      if (ZPolyhedronIncludes(Zpol, temp) == True) {
	ZPolyhedron_Free (Zpol);
	Added = True; 
	break;
      }
      else if (ZPolyhedronIncludes(temp, Zpol) == True) {
	if (temp == Head) {
	  Zpol->next = temp->next;
	  Head = Zpol;
	  ZPolyhedron_Free (temp);
	  Added = True;
	  break;
	}	
	temp1->next = Zpol;
	Zpol->next = temp->next;
	ZPolyhedron_Free (temp);
	Added = True;
	break ;
      }
      temp1 = temp ;
    }
    if(Added == True)
      continue ; 
    for(temp = Head; temp != NULL; temp = temp->next) {
      if(sameLattice(temp->Lat, Zpol->Lat) == True) {
	Polyhedron *Union;
	
	Union = DomainUnion (temp->P,Zpol->P,MAXNOOFRAYS);
	if (!Union)
	  fprintf (stderr,"\n In AddZPolytoZDomain: Out of memory\n");
	else {
	  Domain_Free(temp->P);
	  temp->P = Union;
	  Added = True;
	  ZPolyhedron_Free(Zpol);
	}
	break ;	
      }
      temp1 = temp;
    }
    if (Added == False) 
      temp1->next = Zpol;
  }
  return Head ;  
} /* AddZPolytoZDomain */

/*
 * Return the empty Z-polyhedron of dimension 'dimension' 
 */
ZPolyhedron *EmptyZPolyhedron(int dimension) {

  ZPolyhedron *Zpol;
  Lattice *E ;
  Polyhedron *P; 
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen ("_debug", "a");
  fprintf (fp, "\nEntered EMPTYZPOLYHEDRON\n");
  fclose (fp);
#endif

  E = EmptyLattice(dimension+1);
  P = Empty_Polyhedron(dimension);
  Zpol = ZPolyhedron_Alloc(E,P);
  Matrix_Free((Matrix *) E);
  Domain_Free(P);
  return Zpol;
} /* EmptyZPolyhedron */

/* 
 * Given Z-domains 'A' and 'B', return True if A is included in 'B', otherwise
 * return False.
 */
Bool ZDomainIncludes(ZPolyhedron *A, ZPolyhedron *B) {
  
  ZPolyhedron *Diff;
  Bool ret = False;
  
  Diff = ZDomainDifference(A,B);
  if (isEmptyZPolyhedron(Diff))
    ret = True;
  
  ZDomain_Free(Diff);
  return ret;
} /* ZDomainIncludes */

/* 
 * Given Z-polyhedra 'A' and 'B', return True if 'A' is included in 'B', 
 * otherwise return False 
 */
Bool ZPolyhedronIncludes(ZPolyhedron *A, ZPolyhedron *B) {

  Polyhedron *Diff = NULL ;
  Bool retval = False;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug","a");
  fprintf(fp,"\nEntered ZPOLYHEDRONINCLUDES\n");
  fclose(fp);
#endif
  
  if (LatticeIncludes(A->Lat, B->Lat) == True) {
    Polyhedron *ImageA, *ImageB ;
    
    ImageA = DomainImage(A->P,A->Lat,MAXNOOFRAYS);
    ImageB = DomainImage(B->P,B->Lat,MAXNOOFRAYS);
      
    Diff = DomainDifference(ImageA, ImageB, MAXNOOFRAYS);
    if(emptyQ (Diff))
      retval = True ;
    
    Domain_Free (ImageA);
    Domain_Free (ImageB);
    Domain_Free (Diff);
  }  
  return retval;
} /* ZPolyhedronIncludes */ 

/*
 * Print the contents of a Z-domain 'A' 
 */
void ZDomainPrint(FILE *fp,char *format,ZPolyhedron *A) {
  
#ifdef DOMDEBUG
  FILE *fp1;
  fp1 = fopen("_debug", "a");
  fprintf(fp1,"\nEntered ZDOMAINPRINT\n");
  fclose(fp1);
#endif
  
  ZPolyhedronPrint(fp,format,A);
  if (A->next != NULL) {
    fprintf(fp,"\nUNIONED with\n");
    ZDomainPrint(fp,format,A->next);
  }
  return;
} /* ZDomainPrint */  

/*
 * Print the contents of a ZPolyhderon 'A'
 */
static void ZPolyhedronPrint (FILE *fp, char *format, ZPolyhedron *A) {
  
  if (A == NULL)
    return ;
  fprintf(fp,"\nZPOLYHEDRON: Dimension %d \n",A->Lat->NbRows-1);
  fprintf(fp, "\nLATTICE: \n");
  Matrix_Print(fp,format,(Matrix *)A->Lat);
  Polyhedron_Print(fp,format,A->P);  
  return;
} /* ZPolyhedronPrint */

/*
 * Return the Z-domain union of the Z-domain 'A' and 'B'. The dimensions of the
 * Z-domains 'A' and 'B' must be equal. All the Z-polyhedra of the resulting 
 * union are expected to be in Canonical forms.
 */
ZPolyhedron *ZDomainUnion (ZPolyhedron *A, ZPolyhedron *B) {
  
  ZPolyhedron *Result = NULL, *temp;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ZDOMAINUNION\n");
  fclose(fp);
#endif

  for(temp = A; temp != NULL; temp = temp->next)
    Result = AddZPolytoZDomain(temp, Result); 
  for(temp = B; temp != NULL; temp = temp->next )
    Result = AddZPolytoZDomain(temp, Result);
  return Result;
} /* ZDomainUnion */

/* 
 * Return the Z-domain intersection of the Z-domains 'A' and 'B'.The dimensions
 * of domains 'A' and 'B' must be equal. 
 */
ZPolyhedron *ZDomainIntersection (ZPolyhedron *A, ZPolyhedron *B) {
  
  ZPolyhedron *Result = NULL, *tempA = NULL, *tempB = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ZDOMAININTERSECTION\n");
  fclose(fp);
#endif  
  
  for(tempA = A; tempA != NULL; tempA = tempA->next)
    for(tempB = B; tempB != NULL; tempB = tempB->next) {
      ZPolyhedron *Zpol;
      Zpol = ZPolyhedronIntersection(tempA, tempB);
      Result = AddZPolytoZDomain(Zpol, Result );
      ZPolyhedron_Free (Zpol);
    }
  if (Result == NULL)
    return EmptyZPolyhedron (A->Lat->NbRows-1);
  return Result;
} /* ZDomainIntersection */

/*
 * Return the Z-domain difference of the domains 'A' and 'B'. The dimensions of
 * the Z-domains 'A' and 'B' must be equal. Note that the difference of two 
 * Z-polyhedra is a Union of Z-polyhedra. The algorithms is as given below :-
 * Algorithm: (Given Z-domains A and B)
 *           Result <-- NULL
 *           for every Z-polyhderon Zpoly of A {
 *               temp <-- Zpoly;
 *               for every Z-polyhderon Z1 of B
 *                  temp = temp - Z1;
 *               }
 *           Add temp to Result;
 *           return;
 */
ZPolyhedron *ZDomainDifference(ZPolyhedron  *A, ZPolyhedron *B) {
  
  ZPolyhedron *Result = NULL, *tempA = NULL, *tempB = NULL;
  ZPolyhedron *templist, *res, *i, *j;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ZDOMAINDIFFERENCE\n");
  fclose(fp);
#endif
  
  if (A->Lat->NbRows != B->Lat->NbRows) {
    fprintf(stderr, "In ZDomainDifference : The Input ZDomains");
    fprintf(stderr, "Do not have the compatible dimensions\n");
    fprintf(stderr, "ZDomainDiffernce not performed\n");
    return NULL;
  }
  
  for(tempA = A; tempA != NULL; tempA = tempA->next) {
    ZPolyhedron *temp = NULL;
    temp = ZPolyhedron_Copy(tempA);
    
    for(tempB = B; tempB != NULL; tempB = tempB->next) {
      templist = NULL; res = NULL;
      for(i = temp; i != NULL; i = i->next) {
	i=temp;
	res = ZPolyhedronDifference(i,tempB);
	for (j = res; j != NULL; j = j->next )
	  templist = AddZPoly2ZDomain(j,templist);
	ZDomain_Free(res);
      }
      ZDomain_Free (temp);
      temp = NULL; 
      for(i = templist; i != NULL; i = i->next)
	temp = AddZPoly2ZDomain(i, temp);
      ZDomain_Free (templist);
    }
    for(i = temp; i != NULL; i = i->next)
      Result = AddZPolytoZDomain(i, Result);
    ZDomain_Free(temp);    
  }
  if (Result==NULL)
    return (EmptyZPolyhedron(A->Lat->NbRows-1));
  return Result;
} /* ZDomainDifference */

/*
 * Return the image of the Z-domain 'A' under the invertible, affine, rational
 * transformation function 'Func'. The matrix representing the function 'Func'
 * must be non-singular and the number of rows of the function must be equal
 * to the number of rows in the matrix representing the lattice of 'A'. 
 * Note:: Image((Z1 U Z2),F) = Image(Z1,F) U Image(Z2 U F).
 */
ZPolyhedron *ZDomainImage (ZPolyhedron *A, Matrix *Func) {

  ZPolyhedron *Result = NULL, *temp;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen ("_debug", "a");
  fprintf (fp, "\nEntered ZDOMAINIMAGE\n");
  fclose (fp);
#endif  
  
  for(temp = A; temp != NULL; temp = temp->next) {
    ZPolyhedron *Zpol;
    Zpol =  ZPolyhedronImage (temp, Func);
    Result = AddZPolytoZDomain (Zpol, Result);
    ZPolyhedron_Free (Zpol);
  }  
  if(Result == NULL)
    return EmptyZPolyhedron(A->Lat->NbRows-1);  
  return Result;
} /* ZDomainImage */ 

/*
 * Return the preimage of the Z-domain 'A' under the invertible, affine, ratio-
 * nal transformation 'Func'. The number of rows of the matrix representing 
 * the function 'Func' must be equal to the number of rows of the matrix repr-
 * senting the lattice of 'A'.  
 */
ZPolyhedron *ZDomainPreimage (ZPolyhedron *A, Matrix *Func) {

  ZPolyhedron *Result = NULL, *temp ;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ZDOMAINPREIMAGE\n");
  fclose(fp);
#endif
  
  if (A->Lat->NbRows != Func->NbRows) {
    fprintf(stderr,"\nError : In ZDomainPreimage, ");
    fprintf(stderr,"Incompatible dimensions of ZPolyhedron ");
    fprintf(stderr,"and the Function \n");
    return(EmptyZPolyhedron(Func->NbColumns-1)); 
  }  
  for(temp = A; temp != NULL; temp = temp->next) {	 
    ZPolyhedron *Zpol;
    Zpol = ZPolyhedronPreimage(temp, Func);
    Result = AddZPolytoZDomain(Zpol, Result);
    ZPolyhedron_Free(Zpol);
  }  
  if (Result == NULL)
    return(EmptyZPolyhedron(Func->NbColumns-1));
  return Result;   
} /* ZDomainPreimage */ 

/*
 * Return the Z-polyhedron intersection of the Z-polyhedra 'A' and 'B'. 
 * Note: If Z1 = L1 (intersect) P1 and Z2 = L2 (intersect) P2, then 
 *     Z1 (intersect) Z2 = (L1 (intersect) L2) (intersect) (P1 (intersect) P2) 
 */
static ZPolyhedron *ZPolyhedronIntersection(ZPolyhedron *A, ZPolyhedron *B) {
  
  ZPolyhedron *Result = NULL;
  Lattice *LInter;
  Polyhedron *PInter, *ImageA, *ImageB, *PreImage;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug","a");
  fprintf(fp,"\nEntered ZPOLYHEDRONINTERSECTION\n");
  fclose(fp);
#endif  
  
  LInter = LatticeIntersection(A->Lat,B->Lat);
  if(isEmptyLattice(LInter) == True) {
    ZPolyhedron_Free (Result);
    Matrix_Free(LInter);
    return (EmptyZPolyhedron(A->Lat->NbRows-1));
  }  
  ImageA = DomainImage(A->P,A->Lat,MAXNOOFRAYS);
  ImageB = DomainImage(B->P,B->Lat,MAXNOOFRAYS);
  PInter = DomainIntersection(ImageA,ImageB,MAXNOOFRAYS); 
  if (emptyQ(PInter))
    Result = EmptyZPolyhedron(LInter->NbRows-1);
  else { 
    PreImage = DomainPreimage(PInter,(Matrix *)LInter,MAXNOOFRAYS);
    Result = ZPolyhedron_Alloc(LInter, PreImage);    
    Domain_Free (PreImage);
  }  
  Matrix_Free(LInter);
  Domain_Free(PInter); 
  Domain_Free(ImageA);
  Domain_Free(ImageB);  
  return Result ;
} /* ZPolyhedronIntersection */ 

/* 
 * Return the difference of the two Z-polyhedra 'A' and 'B'. Below is the 
 * procedure to find the difference of 'A' and 'B' :-
 * Procedure: 
 *     Let A = L1 (intersect) P1' and B = L2 (intersect) P2' where 
 *     (P1' = DomImage(P1,L1) and P2' = DomImage(P2,L2)). Then 
 *     A-B = L1 (intersect) (P1'-P2') Union 
 *           (L1-L2) (intersect) (P1' (intersect) P2')
 */ 
static ZPolyhedron *ZPolyhedronDifference(ZPolyhedron *A, ZPolyhedron *B) {
  
  ZPolyhedron *Result = NULL ;
  LatticeUnion *LatDiff, *temp;
  Polyhedron *DomDiff, *DomInter, *PreImage, *ImageA, *ImageB;
  Bool flag = False;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ZPOLYHEDRONDIFFERENCE\n");
  fclose(fp);
#endif

  if(isEmptyZPolyhedron (A))
    return NULL;
  if(isEmptyZPolyhedron (B)) {
    Result = ZDomain_Copy (A);
    return Result;
  }  
  ImageA = DomainImage(A->P,(Matrix *)A->Lat,MAXNOOFRAYS);
  ImageB = DomainImage(B->P,(Matrix *)B->Lat,MAXNOOFRAYS);  
  DomDiff = DomainDifference(ImageA,ImageB,MAXNOOFRAYS);
  if (emptyQ (DomDiff))
    flag = True;
  else {
    ZPolyhedron *Z;
    PreImage = DomainPreimage(DomDiff,A->Lat,MAXNOOFRAYS);
    Z = ZPolyhedron_Alloc(A->Lat,PreImage);    
    Result = AddZPolytoZDomain(Z,Result);
  }  
  if (flag == True)  /* DomDiff = NULL; DomInter = A */
    DomInter = Domain_Copy(ImageA);
  else {
    DomInter = DomainIntersection(ImageA,ImageB,MAXNOOFRAYS);
    if (emptyQ(DomInter)) {
      if (flag == True)
	return (EmptyZPolyhedron(A->Lat->NbRows-1));
      else
	return Result;
    }
  }  
  LatDiff = LatticeDifference(A->Lat, B->Lat);
  if(LatDiff == NULL)
    if(flag == True )
      return(EmptyZPolyhedron (A->Lat->NbRows-1));
  
  while (LatDiff != NULL) {
    ZPolyhedron *tempZ = NULL;
    
    PreImage = DomainPreimage(DomInter, LatDiff->M, MAXNOOFRAYS);    
    tempZ = ZPolyhedron_Alloc(LatDiff->M, PreImage);
    Domain_Free(PreImage);
    Result = AddZPoly2ZDomain(tempZ,Result);
    ZPolyhedron_Free(tempZ);
    temp = LatDiff;
    LatDiff = LatDiff->next;  
    Matrix_Free ((Matrix *) temp->M);
    free (temp);
  }  
  Domain_Free (DomInter);
  Domain_Free (DomDiff);
  return Result;
} /* ZPolyhedronDifference */ 

/*
 * Return the image of the Z-polyhedron 'ZPol' under the invertible, affine, 
 * rational transformation function 'Func'. The matrix representing the funct-
 * ion must be non-singular and the number of rows of the function must be 
 * equal to the number of rows in the matrix representing the lattice of 'ZPol'
 * Algorithm: 
 *         1)  Let ZPol = L (intersect) Q
 *         2)  L1 = LatticeImage(L,F)
 *         3)  Q1 = DomainImage(Q,F)
 *         4)  Z1 = L1(Inverse(L1)*Q1)
 *         5)  Return Z1
 */
static ZPolyhedron *ZPolyhedronImage(ZPolyhedron *ZPol,Matrix *Func) {
 
  ZPolyhedron *Result = NULL ;
  Matrix *LatIm ;
  Polyhedron *Pol, *PolImage ; 
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered ZPOLYHEDRONIMAGE\n");
  fclose(fp);
#endif
  
  if ((Func->NbRows != ZPol->Lat->NbRows) || (Func->NbColumns != ZPol->Lat->NbColumns)) {
    fprintf (stderr, "In ZPolImage - The Function, is not compatible with the ZPolyhedron\n");
    return NULL;
  }
  LatIm = LatticeImage(ZPol->Lat,Func);
  if (isEmptyLattice(LatIm)) {
    Matrix_Free(LatIm);
    return NULL;
  }
  Pol = DomainImage(ZPol->P,ZPol->Lat,MAXNOOFRAYS);
  PolImage = DomainImage(Pol,Func,MAXNOOFRAYS);
  Domain_Free(Pol);
  if(emptyQ(PolImage)) {
    Matrix_Free (LatIm);
    Domain_Free (PolImage);
    return NULL;
  } 
  Pol = DomainPreimage(PolImage,LatIm,MAXNOOFRAYS);
  Result = ZPolyhedron_Alloc(LatIm,Pol);
  Domain_Free(Pol);
  Domain_Free(PolImage);
  Matrix_Free(LatIm);
  return Result;  
} /* ZPolyhedronImage */ 

/*
 * Return the preimage of the Z-polyhedron 'Zpol' under an affine transformati-
 * on function 'G'. The number of rows of matrix representing the function 'G',
 * must be equal to the number of rows of the matrix representing the lattice 
 * of Z1. 
 * Algorithm: 
 *            1) Let Zpol = L (intersect) Q
 *            2) L1 =LatticePreimage(L,F);
 *            3) Q1 = DomainPreimage(Q,F);
 *            4) Z1 = L1(Inverse(L1)*Q1);
 *            5) Return Z1
 */
static ZPolyhedron *ZPolyhedronPreimage(ZPolyhedron *Zpol, Matrix *G) {
  
  Lattice *Latpreim;
  Polyhedron *Qprime, *Q, *Polpreim;
  ZPolyhedron *Result;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug","a");
  fprintf(fp,"\nEntered ZPOLYHEDRONPREIMAGE\n");
  fclose(fp);
#endif
  
  if(G->NbRows != Zpol->Lat->NbRows) {
    fprintf(stderr,"\nIn ZPolyhedronPreimage: Error, The dimensions of the ");
    fprintf(stderr,"function are not compatible with that of the Zpolyhedron");
    return EmptyZPolyhedron(G->NbColumns-1);
  }
  Q = DomainImage(Zpol->P,Zpol->Lat,MAXNOOFRAYS);
  Polpreim = DomainPreimage(Q,G,MAXNOOFRAYS);
  if (emptyQ(Polpreim))
    Result = NULL;
  else {
    Latpreim = LatticePreimage(Zpol->Lat,G);
    if(isEmptyLattice(Latpreim))
      Result = NULL;
    else {
      Qprime = DomainPreimage(Polpreim, Latpreim, MAXNOOFRAYS);
      Result = ZPolyhedron_Alloc(Latpreim, Qprime);
      Domain_Free(Qprime);
    }
    Matrix_Free(Latpreim);
  }  
  Domain_Free(Q);
  return Result; 
} /* ZPolyhedronPreimage */ 

/* 
 * Return the Z-polyhderon 'Zpol' in canonical form: 'Result' (for the Z-poly-
 * hedron in canonical form) and Basis 'Basis' (for the basis with respect to 
 * which 'Result' is in canonical form.   
 */
void CanonicalForm(ZPolyhedron *Zpol,ZPolyhedron **Result,Matrix **Basis) {

  Matrix *B1 = NULL, *B2=NULL, *T1 , *B2inv;
  int i, l1, l2;
  Value tmp;
  Polyhedron *Image, *ImageP;
  Matrix *H, *U, *temp, *Hprime, *Uprime, *T2;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered CANONICALFORM\n");
  fclose(fp);
#endif
  
  if(isEmptyZPolyhedron (Zpol)) {
    Basis[0] = Identity(Zpol->Lat->NbRows);
    Result[0] = ZDomain_Copy (Zpol);
    return ;
  }
  value_init(tmp);
  l1 = FindHermiteBasisofDomain(Zpol->P,&B1);
  Image = DomainImage (Zpol->P,(Matrix *)Zpol->Lat,MAXNOOFRAYS);
  l2 = FindHermiteBasisofDomain(Image,&B2);
    
  if (l1 != l2)
    fprintf(stderr,"In CNF : Something wrong with the Input Zpolyhedra \n"); 
  
  B2inv = Matrix_Alloc(B2->NbRows, B2->NbColumns);
  temp = Matrix_Copy(B2);
  Matrix_Inverse(temp,B2inv);
  Matrix_Free(temp);
  
  temp = Matrix_Alloc(B2inv->NbRows,Zpol->Lat->NbColumns);
  T1 = Matrix_Alloc(temp->NbRows,B1->NbColumns);
  Matrix_Product(B2inv,(Matrix *)Zpol->Lat,temp);
  Matrix_Product(temp,B1,T1);  
  Matrix_Free(temp);
  
  T2 = ChangeLatticeDimension(T1,l1);
  temp = ChangeLatticeDimension(T2,T2->NbRows+1);

  /* Adding the affine part */
  for(i = 0; i < l1; i ++)
    value_assign(temp->p[i][temp->NbColumns-1],T1->p[i][T1->NbColumns-1]);
  
  AffineHermite(temp,&H,&U);
  Hprime = ChangeLatticeDimension(H,Zpol->Lat->NbRows);
  
  /* Exchanging the Affine part */
  for(i = 0; i < l1; i ++) {
    value_assign(tmp,Hprime->p[i][Hprime->NbColumns-1]);
    value_assign(Hprime->p[i][Hprime->NbColumns-1],Hprime->p[i][H->NbColumns-1]);
    value_assign(Hprime->p[i][H->NbColumns-1],tmp);
  }
  Uprime = ChangeLatticeDimension(U,Zpol->Lat->NbRows);
  
  /* Exchanging the Affine part */
  for (i = 0;i < l1; i++) {
    value_assign(tmp,Uprime->p[i][Uprime->NbColumns-1]);
    value_assign(Uprime->p[i][Uprime->NbColumns-1],Uprime->p[i][U->NbColumns-1]);
    value_assign(Uprime->p[i][U->NbColumns-1],tmp);
  }    
  Polyhedron_Free (Image);
  Matrix_Free (B2inv);
  B2inv = Matrix_Alloc(B1->NbRows, B1->NbColumns);
  Matrix_Inverse(B1,B2inv);
  ImageP = DomainImage(Zpol->P, B2inv, MAXNOOFRAYS);
  Matrix_Free(B2inv);
  Image = DomainImage(ImageP, Uprime, MAXNOOFRAYS);
  Domain_Free(ImageP);
  Result[0] = ZPolyhedron_Alloc(Hprime, Image);
  Basis[0] = Matrix_Copy(B2); 
  
  /* Free the variables */
  Polyhedron_Free (Image);
  Matrix_Free (B1);
  Matrix_Free (B2);
  Matrix_Free (temp);
  Matrix_Free (T1);
  Matrix_Free (T2);
  Matrix_Free (H);
  Matrix_Free (U);
  Matrix_Free (Hprime);
  Matrix_Free (Uprime);
  value_clear(tmp);
  return;
} /* CanonicalForm */ 

/*
 * Given a Z-polyhedron 'A' in which the Lattice is not integral, return the
 * Z-polyhedron which contains all the integral points in the input lattice.
 */
ZPolyhedron *IntegraliseLattice(ZPolyhedron *A) {
 
  ZPolyhedron *Result;
  Lattice *M = NULL, *Id;
  Polyhedron *Im = NULL, *Preim = NULL;

#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug","a");
  fprintf(fp,"\nEntered INTEGRALISELATTICE\n");
  fclose(fp);
#endif
  
  Im = DomainImage(A->P,A->Lat,MAXNOOFRAYS);
  Id = Identity(A->Lat->NbRows);
  M = LatticeImage(Id, A->Lat);
  if (isEmptyLattice(M))
    Result = EmptyZPolyhedron(A->Lat->NbRows-1);
  else {
    Preim = DomainPreimage(Im,M,MAXNOOFRAYS);
    Result = ZPolyhedron_Alloc(M,Preim);
  }  
  Matrix_Free(M);
  Domain_Free(Im);
  Domain_Free(Preim);  
  return Result;
} /* IntegraliseLattice */ 

/* 
 * Return the simplified representation of the Z-domain 'ZDom'. It attempts to 
 * convexize unions of polyhedra when they correspond to the same lattices and
 * to simplify union of lattices when they correspond to the same polyhdera. 
 */
ZPolyhedron *ZDomainSimplify(ZPolyhedron *ZDom) {
  
  ZPolyhedron *Ztmp, *Result;
  ForSimplify *Head, *Prev, *Curr;
  ZPolyhedron *ZDomHead, *Emp;
  
  if (ZDom == NULL) {
    fprintf(stderr,"\nError in ZDomainSimplify - ZDomHead = NULL\n");
    return NULL;
  }  
  if (ZDom->next == NULL)
    return (ZPolyhedron_Copy (ZDom));  
  Emp = EmptyZPolyhedron(ZDom->Lat->NbRows-1);
  ZDomHead = ZDomainUnion(ZDom, Emp);
  ZPolyhedron_Free(Emp);  
  Head = NULL;
  Ztmp = ZDomHead;
  do {
    Polyhedron *Img;
    Img = DomainImage(Ztmp->P,Ztmp->Lat,MAXNOOFRAYS);
    for(Curr = Head; Curr != NULL; Curr = Curr->next) {
      Polyhedron *Diff1;
      Bool flag = False;
      
      Diff1 = DomainDifference(Img,Curr->Pol,MAXNOOFRAYS);
      if (emptyQ(Diff1)) {
	Polyhedron *Diff2;

	Diff2 = DomainDifference(Curr->Pol,Img,MAXNOOFRAYS); 
	if (emptyQ(Diff2))
	  flag = True;	
	Domain_Free(Diff2);
      }      
      Domain_Free (Diff1);
      if (flag == True) {
	LatticeUnion *temp;	
	
	temp = (LatticeUnion *)malloc(sizeof(LatticeUnion));
	temp->M = (Lattice *)Matrix_Copy((Matrix *)Ztmp->Lat); 
	temp->next = Curr->LatUni;
	Curr->LatUni = temp;
	break;
      }
    }
    if(Curr == NULL) {
      Curr = (ForSimplify *)malloc(sizeof(ForSimplify));      
      Curr->Pol = Domain_Copy(Img);
      Curr->LatUni = (LatticeUnion *)malloc(sizeof(LatticeUnion));
      Curr->LatUni->M = (Lattice *)Matrix_Copy((Matrix *)Ztmp->Lat); 
      Curr->LatUni->next = NULL;
      Curr->next = Head;
      Head = Curr;
    }   
    Domain_Free (Img);
    Ztmp = Ztmp->next;
  } while(Ztmp != NULL);
  
  for (Curr = Head; Curr != NULL; Curr = Curr->next)
    Curr->LatUni = LatticeSimplify(Curr->LatUni);  
  Result = NULL;
  for(Curr = Head; Curr != NULL; Curr = Curr->next) {
    LatticeUnion *L;    
    for(L = Curr->LatUni; L != NULL; L = L->next) {
      Polyhedron *Preim;
      ZPolyhedron *Zpol;
      
      Preim = DomainPreimage(Curr->Pol,L->M,MAXNOOFRAYS);
      Zpol = ZPolyhedron_Alloc(L->M, Preim);
      Zpol->next = Result;
      Result = Zpol;
      Domain_Free(Preim);
    }
  }  
  Curr = Head;
  while (Curr != NULL) {
    Prev = Curr;
    Curr = Curr->next;     
    LatticeUnion_Free(Prev->LatUni);
    Domain_Free(Prev->Pol);
    free(Prev);
  }
  return Result;
} /* ZDomainSimplify */ 

ZPolyhedron *SplitZpolyhedron(ZPolyhedron *ZPol,Lattice *B) {
 
  Lattice *Intersection = NULL;
  Lattice *B1 = NULL, *B2 = NULL, *newB1 = NULL, *newB2 = NULL;
  Matrix *U = NULL,*M1 = NULL, *M2 = NULL, *M1Inverse = NULL,*MtProduct = NULL;
  Matrix *Vinv, *V , *temp, *DiagMatrix ;
  Matrix *H , *U1 , *X, *Y ;
  ZPolyhedron *zpnew, *Result;
  LatticeUnion *Head = NULL, *tempHead = NULL;
  int i;
  Value k;
  
#ifdef DOMDEBUG
  FILE *fp;
  fp = fopen("_debug", "a");
  fprintf(fp,"\nEntered SplitZpolyhedron \n"); 
  fclose(fp);
#endif

  
  if (B->NbRows != B->NbColumns) { 
    fprintf(stderr,"\n SplitZpolyhedron : The Input Matrix B is not a proper Lattice \n");
    return NULL;
  }
  
  if (ZPol->Lat->NbRows != B->NbRows) {
    fprintf(stderr,"\nSplitZpolyhedron : The Lattice in Zpolyhedron and B have ");
    fprintf(stderr,"incompatible dimensions \n");
    return NULL;
  }
  
  if (isinHnf (ZPol->Lat) != True) {
    AffineHermite(ZPol->Lat,&H,&U1);
    X = Matrix_Copy(H);    
    Matrix_Free(U1);
    Matrix_Free(H);
  }
  else
    X = Matrix_Copy(ZPol->Lat);
  
  if (isinHnf(B) != True) {
    AffineHermite(B,&H,&U1);
    Y = Matrix_Copy(H);   
    Matrix_Free(H);
    Matrix_Free(U1);
  }
  else
    Y = Matrix_Copy(B);  
  if (isEmptyLattice(X)) {
    return NULL;  
  }

  Head=Lattice2LatticeUnion(X,Y);

/* If the spliting operation can't be done the result is the original Zplyhedron. */

  if (Head == NULL) {
    Matrix_Free(X);
    Matrix_Free(Y);
    return ZPol;
  }  


  Result=NULL;

  if (Head)
   while(Head)
    {
      tempHead = Head;
      Head = Head->next;  
      zpnew=ZPolyhedron_Alloc(tempHead->M,ZPol->P);
      Result=AddZPoly2ZDomain(zpnew,Result);
      ZPolyhedron_Free(zpnew);
      tempHead->next = NULL; 
      free(tempHead);  
    }

  return Result;
}



