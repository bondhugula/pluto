/** 
 * $Id: testCompressParms.c,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 * 
 * Test routines for kernel/compress_parms.c functions
 * @author B. Meister, 3/2006
 * 
 */

#include <polylib/polylib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef dbg
#undef dbg
#endif
#define dbg 1

#define TEST(a) if (isOk = a) { \
                  printf(#a" tested ok.\n"); \
                } \
                else { \
                  printf(#a" NOT OK\n"); \
                } 

#define maxRays 200

char * origNames[] = {"n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"};

int main(int argc, char ** argv) {
  int isOk = 0;
  Matrix * A, * B;
  if (argc>1) {
    printf("Warning: No arguments taken into account: testing"
	   "remove_parm_eqs().\n");
  }

  A = Matrix_Read();
  B = Matrix_Read();
  TEST( test_Constraints_Remove_parm_eqs(A, B) )
  TEST( test_Polyhedron_Remove_parm_eqs(A, B) )
  TEST( test_Constraints_fullDimensionize(A, B, 4) )
  Matrix_Free(A);
  Matrix_Free(B);
  return (1-isOk);
}


/** extracts the equalities involving the parameters only, try to introduce
    them back and compare the two polyhedra.
    Reads a polyhedron and a context.
 */
int test_Constraints_Remove_parm_eqs(Matrix * A, Matrix * B) {
  int isOk = 1;
  Matrix * M, *C, *Cp, * Eqs, *M1, *C1;
  Polyhedron *Pm, *Pc, *Pcp, *Peqs, *Pint;
  unsigned int * elimParms;
  printf("----- test_Constraints_Remove_parm_eqs() -----\n");
  M1 = Matrix_Copy(A);
  C1 = Matrix_Copy(B);

  M = Matrix_Copy(M1);
  C = Matrix_Copy(C1);

   /* compute the combined polyhedron */
  Pm = Constraints2Polyhedron(M, maxRays);
  Pc = Constraints2Polyhedron(C, maxRays);
  Pcp = align_context(Pc, Pm->Dimension, maxRays);
  Polyhedron_Free(Pc);
  Pc = DomainIntersection(Pm, Pcp, maxRays);
  Polyhedron_Free(Pm);
  Polyhedron_Free(Pcp);
  Matrix_Free(M);
  Matrix_Free(C);

  /* extract the parm-equalities, expressed in the combined space */
  Eqs = Constraints_Remove_parm_eqs(&M1, &C1, 1, &elimParms);

  printf("Removed equalities: \n");
  show_matrix(Eqs); 
  printf("Polyhedron without equalities involving only parameters: \n");
  show_matrix(M1);  
  printf("Context without equalities: \n");
  show_matrix(C1);  
  
  /* compute the supposedly-same polyhedron, using the extracted equalities */
  Pm = Constraints2Polyhedron(M1, maxRays);
  Pcp = Constraints2Polyhedron(C1, maxRays);
  Peqs = align_context(Pcp, Pm->Dimension, maxRays);
  Polyhedron_Free(Pcp);
  Pcp = DomainIntersection(Pm, Peqs, maxRays);
  Polyhedron_Free(Peqs);
  Polyhedron_Free(Pm);
  Peqs = Constraints2Polyhedron(Eqs, maxRays);
  Matrix_Free(Eqs);
  Matrix_Free(M1);
  Matrix_Free(C1);
  Pint = DomainIntersection(Pcp, Peqs, maxRays);
  Polyhedron_Free(Pcp);
  Polyhedron_Free(Peqs);

  /* test their equality */
  if (!PolyhedronIncludes(Pint, Pc)) {
    isOk = 0;
  }
  else {
    if (!PolyhedronIncludes(Pc, Pint)) {
      isOk = 0;
    }
  }
  Polyhedron_Free(Pc);
  Polyhedron_Free(Pint);
  return isOk;
} /* test_Constraints_Remove_parm_eqs() */


/** extracts the equalities holding on the parameters only, try to introduce
    them back and compare the two polyhedra.
    Reads a polyhedron and a context.
 */
int test_Polyhedron_Remove_parm_eqs(Matrix * A, Matrix * B) {
  int isOk = 1;
  Matrix * M, *C;
  Polyhedron *Pm, *Pc, *Pcp, *Peqs, *Pint, *Pint1;
  unsigned int * elimParms;
  printf("----- test_Polyhedron_Remove_parm_eqs() -----\n");

  M = Matrix_Copy(A);
  C = Matrix_Copy(B);

   /* compute the combined polyhedron */
  Pm = Constraints2Polyhedron(M, maxRays);
  Pc = Constraints2Polyhedron(C, maxRays);
  Pcp = align_context(Pc, Pm->Dimension, maxRays);
  Polyhedron_Free(Pc);
  Pint1 = DomainIntersection(Pm, Pcp, maxRays);
  Polyhedron_Free(Pm);
  Polyhedron_Free(Pcp);
  Matrix_Free(M);
  Matrix_Free(C);

  M = Matrix_Copy(A);
  C = Matrix_Copy(B);
  /* extract the parm-equalities, expressed in the combined space */
  Pm = Constraints2Polyhedron(M, maxRays);
  Pc = Constraints2Polyhedron(C, maxRays);
  Matrix_Free(M);
  Matrix_Free(C);
  Peqs = Polyhedron_Remove_parm_eqs(&Pm, &Pc, 1, &elimParms, 200);
  
  /* compute the supposedly-same polyhedron, using the extracted equalities */
  Pcp = align_context(Pc, Pm->Dimension, maxRays);
  Polyhedron_Free(Pc);
  Pc = DomainIntersection(Pm, Pcp, maxRays);
  Polyhedron_Free(Pm);
  Polyhedron_Free(Pcp);
 
  Pint = DomainIntersection(Pc, Peqs, maxRays);
  Polyhedron_Free(Pc);
  Polyhedron_Free(Peqs);

  /* test their equality */
  if (!PolyhedronIncludes(Pint, Pint1)) {
    isOk = 0;
  }
  else {
    if (!PolyhedronIncludes(Pint1, Pint)) {
      isOk = 0;
    }
  }
  Polyhedron_Free(Pint1);
  Polyhedron_Free(Pint);
  return isOk;
} /* test_Polyhedron_remove_parm_eqs() */


/** 
 * Eliminates certain parameters from a vector of values for parameters
 * @param origParms the initial vector of values of parameters
 * @param elimParms the list of parameters to be eliminated in the vector
 * @param newParms the vector of values without the eliminated ones.
 */
void valuesWithoutElim(Matrix * origParms, unsigned int * elimParms, 
		       Matrix ** newParms) {
  unsigned int i, j=0;
  if (*newParms==NULL) {
    *newParms = Matrix_Alloc(1, origParms->NbColumns-elimParms[0]);
  } /* else assume enough space is allocated */
  if (elimParms[0] ==0) {
    for (i=0; i< origParms->NbColumns; i++) {
      value_assign((*newParms)->p[0][i], origParms->p[0][i]);
    }
  }
  for (i=0; i< origParms->NbColumns; i++) {
    if (i!=elimParms[j+1]) {
      value_assign((*newParms)->p[0][i-j], origParms->p[0][i]);
    }
    else {
      j++;
    }
  }
}/* valuesWithoutElim */


/**
 * takes a list of parameter names, a list ofparameters to eliminate, and
 * returns the list of parameters without the eliminated ones.
 * @param parms the original parameter names
 * @param nbParms the number of original parmaeters
 * @param elimParms the array-list of indices of parameters to eliminate (first
 * element set to the number of its elements)
 * @param newParms the returned list of parm names. Allocated if set to NULL,
 * reused if not.
 * @return the number of names in the returned list.
 */
unsigned int namesWithoutElim(char ** parms, unsigned nbParms,
			      unsigned int * elimParms,
			      char *** newParms) {
  unsigned int i, j=0;
  unsigned int newSize = nbParms -elimParms[0];
  if (dbg) {
    printf("Size of the new parm vector: %d\n", newSize);
  }
  if (*newParms==NULL) {
    *newParms = malloc(newSize*sizeof(char *));
  }
  if (elimParms[0]==0) {
    for (i=0; i< nbParms; i++) {
      (*newParms)[i] = strdup(parms[i]);
    }
    return newSize;
  }
  for (i=0; i< nbParms; i++) {
    if (i!=elimParms[j+1]) {
      (*newParms)[i-j] = strdup(parms[i]);
    }
    else {
      j++;
    }
  }
  return newSize;
}


/**
 * Tests Constraints_fullDimensionize by comparing the Ehrhart polynomials 
 * @param A the input set of constraints
 * @param B the corresponding context
 * @param the number of samples to generate for the test
 * @return 1 if the Ehrhart polynomial had the same value for the
 * full-dimensional and non-full-dimensional sets of constraints, for their
 * corresponding sample parameters values.
 */
int test_Constraints_fullDimensionize(Matrix * A, Matrix * B, 
				      unsigned int nbSamples) {
  Matrix * Eqs= NULL, *ParmEqs=NULL, *VL=NULL;
  unsigned int * elimVars=NULL, * elimParms=NULL;
  Matrix * sample, * smallerSample=NULL;
  Matrix * transfSample=NULL;
  Matrix * parmVL=NULL;
  unsigned int i, j, r, nbOrigParms, nbParms;
  Value div, mod, *origVal=NULL, *fullVal=NULL;
  Matrix * VLInv;
  Polyhedron * P, *PC;
  Matrix * M, *C;
  Enumeration * origEP, * fullEP=NULL;
  char ** fullNames = NULL;
  int isOk = 1; /* holds the result */

  /* compute the origial Ehrhart polynomial */
  M = Matrix_Copy(A);
  C = Matrix_Copy(B);
  P = Constraints2Polyhedron(M, maxRays);
  PC = Constraints2Polyhedron(C, maxRays);
  origEP = Polyhedron_Enumerate(P, PC, maxRays, origNames);
  Matrix_Free(M);
  Matrix_Free(C);
  Polyhedron_Free(P);
  Polyhedron_Free(PC);

  /* compute the full-dimensional polyhedron corresponding to A and its Ehrhart
     polynomial */
  M = Matrix_Copy(A);
  C = Matrix_Copy(B);
  nbOrigParms = B->NbColumns-2;
  Constraints_fullDimensionize(&M, &C, &VL, &Eqs, &ParmEqs, 
			       &elimVars, &elimParms, maxRays);
  if ((Eqs->NbRows==0) && (ParmEqs->NbRows==0)) {
    Matrix_Free(M);
    Matrix_Free(C);
    Matrix_Free(Eqs);
    Matrix_Free(ParmEqs);
    free(elimVars);
    free(elimParms);
    return 1;
  }
  nbParms = C->NbColumns-2;
  P = Constraints2Polyhedron(M, maxRays);
  PC = Constraints2Polyhedron(C, maxRays);
  namesWithoutElim(origNames, nbOrigParms, elimParms, &fullNames);
  fullEP = Polyhedron_Enumerate(P, PC, maxRays, fullNames);
  Matrix_Free(M);
  Matrix_Free(C);
  Polyhedron_Free(P);
  Polyhedron_Free(PC);
  
  /* make a set of sample parameter values and compare the corresponding
     Ehrhart polnomials */
  sample = Matrix_Alloc(1,nbOrigParms);
  transfSample = Matrix_Alloc(1, nbParms);
  Lattice_extractSubLattice(VL, nbParms, &parmVL);
  VLInv = Matrix_Alloc(parmVL->NbRows, parmVL->NbRows+1);
  MatInverse(parmVL, VLInv);
  if (dbg) {
    show_matrix(parmVL);
    show_matrix(VLInv);
  }
  srandom(nbSamples);
  value_init(mod);
  value_init(div);
  for (i = 0; i< nbSamples; i++) {
    /* create a random sample */
    for (j=0; j< nbOrigParms; j++) {
      value_set_si(sample->p[0][j], random()%100);
    }
    /* compute the corresponding value for the full-dimensional
       constraints */
    valuesWithoutElim(sample, elimParms, &smallerSample); 
    /* (N' i' 1)^T = VLinv.(N i 1)^T*/
    for (r = 0; r < nbParms; r++) {
      Inner_Product(&(VLInv->p[r][0]), smallerSample->p[0], nbParms,
		    &(transfSample->p[0][r]));
      /* add the constant part */
      value_addto(transfSample->p[0][r], transfSample->p[0][r], 
					 VLInv->p[r][VLInv->NbColumns-2]);
      value_pdivision(div, transfSample->p[0][r], 
			 VLInv->p[r][VLInv->NbColumns-1]);
      value_subtract(mod, transfSample->p[0][r], div);
      /* if the parameters value does not belong to the validity lattice, the
	 Ehrhart polynomial is zero. */
      if (!value_zero_p(mod)) {
	fullEP = Enumeration_zero(nbParms, maxRays);
	break;
      }
    }
    /* compare the two forms of the Ehrhart polynomial.*/
    if (origEP ==NULL) break; /* NULL has loose semantics for EPs */
    origVal = compute_poly(origEP, sample->p[0]);
    fullVal = compute_poly(fullEP, transfSample->p[0]);
    if (!value_eq(*origVal, *fullVal)) {
      isOk = 0;
      printf("EPs don't match. \n Original value = ");
      value_print(stdout, VALUE_FMT, *origVal);
      printf("\n Original sample = [");
      for (j=0; j<sample->NbColumns; j++) {
	value_print(stdout, VALUE_FMT, sample->p[0][j]);
	printf(" ");
      }
      printf("] \n EP = ");
      if(origEP!=NULL) {
	print_evalue(stdout, &(origEP->EP), origNames);
      }
      else {
	printf("NULL");
      }
      printf(" \n Full-dimensional value = ");
      value_print(stdout, P_VALUE_FMT, *fullVal);
      printf("\n full-dimensional sample = [");
      for (j=0; j<sample->NbColumns; j++) {
	value_print(stdout, VALUE_FMT, transfSample->p[0][j]);
	printf(" ");
      }
      printf("] \n EP = ");
      if(origEP!=NULL) {
	print_evalue(stdout, &(origEP->EP), fullNames);
      }
      else {
	printf("NULL");
      }
    }
    if (dbg) {
      printf("\nOriginal value = ");
      value_print(stdout, VALUE_FMT, *origVal);
      printf("\nFull-dimensional value = ");
      value_print(stdout, P_VALUE_FMT, *fullVal);
      printf("\n");
    }
    value_clear(*origVal);
    value_clear(*fullVal);
  }
  value_clear(mod);
  value_clear(div);
  Matrix_Free(sample);
  Matrix_Free(smallerSample);
  Matrix_Free(transfSample);
  Enumeration_Free(origEP);
  Enumeration_Free(fullEP);
  return isOk;
} /* test_Constraints_fullDimensionize */
