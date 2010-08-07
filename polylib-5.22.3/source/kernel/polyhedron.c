/* polyhedron.c
     COPYRIGHT
          Both this software and its documentation are

              Copyright 1993 by IRISA /Universite de Rennes I - France,
              Copyright 1995,1996 by BYU, Provo, Utah
                         all rights reserved.

          Permission is granted to copy, use, and distribute
          for any commercial or noncommercial purpose under the terms
          of the GNU General Public license, version 2, June 1991
          (see file : LICENSING).
*/

/*

1997/12/02 - Olivier Albiez
  Ce fichier contient les fonctions de la polylib de l'IRISA,
  passees en 64bits.
  La structure de la polylib a donc ete modifie pour permettre 
  le passage aux Value. La fonction Chernikova a ete reecrite.

*/

/*

1998/26/02 - Vincent Loechner
  Ajout de nombreuses fonctions, a la fin de ce fichier,
  pour les polyedres parametres 64 bits.
1998/16/03
  #define DEBUG  printf
  tests out of memory
  compatibilite avec la version de doran

*/

#undef POLY_DEBUG		/* debug printf: general functions */
#undef POLY_RR_DEBUG		/* debug printf: Remove Redundants */
#undef POLY_CH_DEBUG		/* debug printf: Chernikova */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <polylib/polylib.h>

#ifdef MAC_OS
  #define abs __abs
#endif

/* WSIZE is the number of bits in a word or int type */ 
#define WSIZE (8*sizeof(int)) 

#define bexchange(a, b, l)\
{\
  char *t = (char *)malloc(l*sizeof(char));\
  memcpy((t), (char *)(a), (int)(l));\
  memcpy((char *)(a), (char *)(b), (int)(l));\
  memcpy((char *)(b), (t), (int)(l));\
  free(t); \
}

#define exchange(a, b, t)\
{ (t)=(a); (a)=(b); (b)=(t); }

/*  errormsg1 is an external function which is usually supplied by the
    calling program (e.g. Domlib.c, ReadAlpha, etc...).
    See errormsg.c for an example of such a function.  */

void errormsg1(char *f , char *msgname, char *msg);

int Pol_status;                    /* error status after operations */

/*
 * The Saturation matrix is defined to be an integer (int type) matrix.
 * It is a boolean matrix which has a row for every constraint and a column
 * for every line or ray. The bits in the binary format of each integer in 
 * the stauration matrix stores the information whether the corresponding
 * constraint is saturated by ray(line) or not.   
 */

typedef struct {
  unsigned int NbRows;
  unsigned int NbColumns;
  int **p;
  int *p_init;
} SatMatrix;

/*
 * Allocate memory space for a saturation matrix. 
 */
static SatMatrix *SMAlloc(int rows,int cols) {
  
  int **q, *p, i;
  SatMatrix *result;
  
  result = (SatMatrix *) malloc (sizeof(SatMatrix));
  if(!result) { 
    errormsg1("SMAlloc", "outofmem", "out of memory space");
    return 0;
  }
  result->NbRows = rows;
  result->NbColumns = cols;
  if(rows == 0 || cols == 0) {
    result->p = NULL;
    return result;
  }
  result->p = q = (int **)malloc(rows * sizeof(int *));
  if(!result->p) {
    errormsg1("SMAlloc", "outofmem", "out of memory space");
    return 0;
  }
  result->p_init = p = (int *)malloc (rows * cols * sizeof (int));
  if(!result->p_init) {
    errormsg1("SMAlloc", "outofmem", "out of memory space");
    return 0;
  }
  for (i=0; i<rows; i++) {
    *q++ = p;
    p += cols;
  }
  return result;
} /* SMAlloc */

/* 
 * Free the memory space occupied by saturation matrix. 
 */ 
static void SMFree (SatMatrix **matrix) {
  SatMatrix *SM = *matrix;

  if (SM) { 
    if (SM->p) {
      free ((char *) SM->p_init);
      free ((char *) SM->p);
    }
    free ((char *) SM);
    *matrix = NULL;
  }
} /* SMFree */

/*
 * Print the contents of a saturation matrix.
 * This function is defined only for debugging purpose. 
 */
static void SMPrint (SatMatrix *matrix) {
  
  int *p;
  int i, j;
  unsigned NbRows, NbColumns;
  
  fprintf(stderr,"%d %d\n",NbRows=matrix->NbRows, NbColumns=matrix->NbColumns);
  for (i=0;i<NbRows;i++) {
    p = *(matrix->p+i);
    for (j=0;j<NbColumns;j++)
      fprintf(stderr, " %10X ", *p++);
    fprintf(stderr, "\n");
  }  
} /* SMPrint */

/* 
 * Compute the bitwise OR of two saturation matrices.
 */
static void SatVector_OR(int *p1,int *p2,int *p3,unsigned length) {
  
  int *cp1, *cp2, *cp3;
  int i;
  
  cp1=p1;
  cp2=p2;
  cp3=p3;
  for (i=0;i<length;i++) {
    *cp3 = *cp1 | *cp2;
    cp3++;
    cp1++;
    cp2++;
  }
} /* SatVector_OR */

/* 
 * Copy a saturation matrix to another (macro definition). 
 */
#define SMVector_Copy(p1, p2, length) \
  memcpy((char *)(p2), (char *)(p1), (int)((length)*sizeof(int)))

/*
 * Initialize a saturation matrix with zeros (macro definition) 
 */
#define SMVector_Init(p1, length) \
  memset((char *)(p1), 0, (int)((length)*sizeof(int)))

/*
 * Defining operations on polyhedron --
 */

/* 
 * Vector p3 is a linear combination of two vectors (p1 and p2) such that 
 * p3[pos] is zero. First element of each vector (p1,p2,p3) is a status 
 * element and is not changed in p3. The value of 'pos' may be 0 however.
 * The parameter 'length' does not include status element one. 
 */
static void Combine(Value *p1, Value *p2, Value *p3, int pos, unsigned length) { 

  Value a1, a2, gcd;
  Value abs_a1,abs_a2,neg_a1;

  /* Initialize all the 'Value' variables */
  value_init(a1); value_init(a2); value_init(gcd);
  value_init(abs_a1); value_init(abs_a2); value_init(neg_a1);
  
  /* a1 = p1[pos] */
  value_assign(a1,p1[pos]); 

  /* a2 = p2[pos] */
  value_assign(a2,p2[pos]); 

  /* a1_abs = |a1| */
  value_absolute(abs_a1,a1);
  
  /* a2_abs = |a2| */
  value_absolute(abs_a2,a2);

  /* gcd  = Gcd(abs(a1), abs(a2)) */
  Gcd(abs_a1,abs_a2,&gcd);

  /* a1 = a1/gcd */
  value_division (a1,a1,gcd);

  /* a2 = a2/gcd */
  value_division (a2,a2,gcd);

  /* neg_a1 = -(a1) */
  value_oppose(neg_a1,a1);

  Vector_Combine(p1+1,p2+1,p3+1,a2,neg_a1,length);
  Vector_Normalize(p3+1,length);

  /* Clear all the 'Value' variables */
  value_clear(a1); value_clear(a2); value_clear(gcd);
  value_clear(abs_a1); value_clear(abs_a2); value_clear(neg_a1);    
  
  return;
} /* Combine */

/* 
 * Return the transpose of the saturation matrix 'Sat'. 'Mat' is a matrix
 * of constraints and 'Ray' is a matrix of ray vectors and 'Sat' is the 
 * corresponding saturation matrix. 
 */
static SatMatrix *TransformSat(Matrix *Mat, Matrix *Ray, SatMatrix *Sat) { 
  
  int i, j, sat_nbcolumns;
  unsigned jx1, jx2, bx1, bx2;
  SatMatrix *result;

  if (Mat->NbRows != 0) 
    sat_nbcolumns = (Mat->NbRows-1) /(sizeof(int)*8) + 1;
  else                  
    sat_nbcolumns = 0;

  result = SMAlloc(Ray->NbRows, sat_nbcolumns);
  SMVector_Init(result->p_init, Ray->NbRows * sat_nbcolumns);

  for(i=0,jx1=0,bx1=MSB; i<Ray->NbRows; i++) { 
    for(j=0,jx2=0,bx2=MSB; j<Mat->NbRows; j++) { 
      if (Sat->p[j][jx1] & bx1) 
        result->p[i][jx2] |= bx2;
      NEXT(jx2,bx2);
    }
    NEXT(jx1, bx1);
  }
  return result;
} /* TransformSat */

/* 
 * Sort the rays (Ray, Sat) into three tiers as used in 'Chernikova' function:
 * NbBid         <= i <  equal_bound    : saturates the constraint
 * equal_bound   <= i <  sup_bound      : verifies the constraint
 * sup_bound     <= i <  NbRay          : does not verify 
 *
 * 'Ray' is the matrix of rays and 'Sat' is the corresponding saturation 
 * matrix. (jx,bx) pair specify the constraint in the saturation matrix. The
 * status element of the 'Ray' matrix holds the saturation value w.r.t the 
 * constraint specified by (jx,bx). Thus
 * Ray->p[i][0]  = 0 -> ray(i) saturates the constraint
 * Ray->p[i][0]  > 0 -> ray(i) verifies  the constraint
 * Ray->p[i][0]  < 0 -> ray(i) doesn't verify the constraint  
 */  
static void RaySort(Matrix *Ray,SatMatrix *Sat,int NbBid,int NbRay,int *equal_bound,int *sup_bound,unsigned RowSize1, unsigned RowSize2, unsigned bx, unsigned jx) {                     

  int inf_bound;
  Value **uni_eq, **uni_sup, **uni_inf;
  int **inc_eq, **inc_sup, **inc_inf;

  /* 'uni_eq' points to the first ray in the ray matrix which verifies a
   * constraint, 'inc_eq' is the corresponding pointer in saturation 
   * matrix. 'uni_inf' points to the first ray (from top) which doesn't 
   * verify a constraint, 'inc_inf' is the corresponding pointer in 
   * saturation matrix. 'uni_sup' scans the ray matrix and 'inc_sup' is 
   * the corresponding pointer in saturation matrix. 'inf_bound' holds the 
   * number of the first ray which does not verify the constraints. 
   */

  *sup_bound = *equal_bound = NbBid;
  uni_sup = uni_eq = Ray->p+NbBid;
  inc_sup = inc_eq = Sat->p+NbBid;
  inf_bound = NbRay;
  uni_inf = Ray->p+NbRay;
  inc_inf = Sat->p+NbRay;
  
  while (inf_bound>*sup_bound) {
    if (value_zero_p(**uni_sup)) {               /* status = satisfy */
      if (inc_eq != inc_sup) {
	Vector_Exchange(*uni_eq,*uni_sup,RowSize1);
	bexchange(*inc_eq,*inc_sup,RowSize2);
      }
      (*equal_bound)++; uni_eq++; inc_eq++;
      (*sup_bound)++; uni_sup++; inc_sup++;
    }
    else {
      *((*inc_sup)+jx)|=bx;
      
      /* if (**uni_sup<0) */
      if (value_neg_p(**uni_sup)) {             /* Status != verify  */
	inf_bound--; uni_inf--; inc_inf--;
	if (inc_inf != inc_sup) {
	  Vector_Exchange(*uni_inf,*uni_sup,RowSize1);
	  bexchange(*inc_inf,*inc_sup,RowSize2);
	}
      }
      else {                                     /* status == verify */
	 (*sup_bound)++; uni_sup++; inc_sup++;
      } 
    }
  }
} /* RaySort */ 

static void SatMatrix_Extend(SatMatrix *Sat, Matrix* Mat, unsigned rows)
{
  int i;
  unsigned cols;
  cols = (Mat->NbRows - 1)/(sizeof(int)*8) + 1;

  Sat->p = (int **)realloc(Sat->p, rows * sizeof(int *));
  if(!Sat->p) {
    errormsg1("SatMatrix_Extend", "outofmem", "out of memory space");
    return;
  }
  Sat->p_init = (int *)realloc(Sat->p_init, rows * cols * sizeof (int));
  if(!Sat->p_init) {
    errormsg1("SatMatrix_Extend", "outofmem", "out of memory space");
    return;
  }
  for (i = 0; i < rows; ++i)
    Sat->p[i] = Sat->p_init + (i * cols);
  Sat->NbRows = rows;
}

static void Matrix_Extend(Matrix *Mat, unsigned NbRows)
{
  Value *p, **q;
  int i,j;

  q = (Value **)realloc(Mat->p, NbRows * sizeof(*q));
  if(!q) {
    errormsg1("Matrix_Extend", "outofmem", "out of memory space");
    return;
  }
  Mat->p = q;
  if (Mat->p_Init_size < NbRows * Mat->NbColumns) {
    p = (Value *)realloc(Mat->p_Init, NbRows * Mat->NbColumns * sizeof(Value));
    if(!p) {
      errormsg1("Matrix_Extend", "outofmem", "out of memory space");
      return;
    }
    Mat->p_Init = p;
    Vector_Set(Mat->p_Init + Mat->NbRows*Mat->NbColumns, 0,
	       Mat->p_Init_size - Mat->NbRows*Mat->NbColumns);
    for (i = Mat->p_Init_size; i < Mat->NbColumns*NbRows; ++i)
	value_init(Mat->p_Init[i]);
    Mat->p_Init_size = Mat->NbColumns*NbRows;
  } else
    Vector_Set(Mat->p_Init + Mat->NbRows*Mat->NbColumns, 0,
	       (NbRows - Mat->NbRows) * Mat->NbColumns);
  for (i=0;i<NbRows;i++) {
    Mat->p[i] = Mat->p_Init + (i * Mat->NbColumns);
  }
  Mat->NbRows = NbRows;
}

/* 
 * Compute the dual of matrix 'Mat' and place it in matrix 'Ray'.'Mat' 
 * contains the constraints (equalities and inequalities) in rows and 'Ray' 
 * contains the ray space (lines and rays) in its rows. 'Sat' is a boolean 
 * saturation matrix defined as Sat(i,j)=0 if ray(i) saturates constraint(j),
 *  otherwise 1. The constraints in the 'Mat' matrix are processed starting at
 * 'FirstConstraint', 'Ray' and 'Sat' matrices are changed accordingly.'NbBid'
 * is the number of lines in the ray matrix and 'NbMaxRays' is the maximum 
 * number of rows (rays) permissible in the 'Ray' and 'Sat' matrix. Return 0 
 * if successful, otherwise return 1.  
 */     
static int Chernikova (Matrix *Mat,Matrix *Ray,SatMatrix *Sat, unsigned NbBid, unsigned NbMaxRays, unsigned FirstConstraint,unsigned dual) {

  unsigned NbRay, Dimension, NbConstraints, RowSize1, RowSize2, sat_nbcolumns;
  int sup_bound, equal_bound, index_non_zero, bound;
  int i, j, k, l, redundant, rayonly, nbcommonconstraints;
  int *Temp, aux;
  int *ip1, *ip2;
  unsigned bx, m, jx;
  Value *p1, *p2, *p3;

#ifdef POLY_CH_DEBUG
  fprintf(stderr, "[Chernikova: Input]\nRay = ");
  Matrix_Print(stderr,0,Ray);
  fprintf(stderr, "\nConstraints = ");
  Matrix_Print(stderr,0,Mat);
  fprintf(stderr, "\nSat = ");
  SMPrint(Sat);
#endif
  
  NbConstraints=Mat->NbRows;
  NbRay = Ray->NbRows;
  Dimension = Mat->NbColumns-1;         /* Homogeneous Dimension */
  sat_nbcolumns=Sat->NbColumns;
  
  RowSize1=(Dimension+1);
  RowSize2=sat_nbcolumns * sizeof(int);

  Temp=(int *)malloc(RowSize2);
  if(!Temp) {	
    errormsg1("Chernikova", "outofmem", "out of memory space");
    return 0;
  }
  CATCH(any_exception_error) {

  /* 
   * In case of overflow, free the allocated memory!
   * Rethrow upwards the stack to forward the exception.
   */
    free(Temp);
    RETHROW();
  }
  TRY {
    jx = FirstConstraint/WSIZE;    
    bx = MSB; bx >>= FirstConstraint%WSIZE;
    for (k=FirstConstraint; k<NbConstraints; k++) {
      
      /* Set the status word of each ray[i] to ray[i] dot constraint[k] */
      /* This is equivalent to evaluating each ray by the constraint[k] */
      /* 'index_non_zero' is assigned the smallest ray index which does */
      /* not saturate the constraint.                                   */
      
      index_non_zero = NbRay;
      for (i=0; i<NbRay; i++) { 
	p1 = Ray->p[i]+1;
	p2 = Mat->p[k]+1;
	p3 = Ray->p[i];
      	
	/* *p3 = *p1 * *p2 */     
	value_multiply(*p3,*p1,*p2);
	p1++; p2++;
	for (j=1; j<Dimension; j++) {	
	  
	  /* *p3 +=  *p1 * *p2 */
	  value_addmul(*p3, *p1, *p2);
	  p1++; p2++;
	}
	if (value_notzero_p(*p3) && (i<index_non_zero)) 
	  index_non_zero=i;
      }
      
#ifdef POLY_CH_DEBUG
      fprintf(stderr, "[Chernikova: A]\nRay = ");
      Matrix_Print(stderr,0,Ray);
      fprintf(stderr, "\nConstraints = ");
      Matrix_Print(stderr,0,Mat);
      fprintf(stderr, "\nSat = ");
      SMPrint (Sat);
#endif
      
      /* Find a bidirectional ray z such that cz <> 0 */  
      if (index_non_zero<NbBid) {
	
	/* Discard index_non_zero bidirectional ray */    
	NbBid--;
	if (NbBid!=index_non_zero) 
	  Vector_Exchange(Ray->p[index_non_zero],Ray->p[NbBid],RowSize1);	

#ifdef POLY_CH_DEBUG	
	fprintf(stderr,"************\n");
	for(i=0;i<RowSize1;i++) {
	  value_print(stderr,P_VALUE_FMT,Ray->p[index_non_zero][i]);
	}  
	fprintf(stderr,"\n******\n");
	for(i=0;i<RowSize1;i++) {
	  value_print(stderr,P_VALUE_FMT,Ray->p[NbBid][i]);
	}
	fprintf(stderr,"\n*******\n");
#endif

	/* Compute the new lineality space */    
	for (i=0; i<NbBid; i++)
	  if (value_notzero_p(Ray->p[i][0]))
	    Combine(Ray->p[i],Ray->p[NbBid],Ray->p[i],0,Dimension);

	/* Add the positive part of index_non_zero bidirectional ray to  */
	/* the set of unidirectional rays                                */
	
	if (value_neg_p(Ray->p[NbBid][0])) {
	  p1=Ray->p[NbBid]; 
	  for (j=0;j<Dimension+1; j++) { 
	    
	    /* *p1 = - *p1 */	
	    value_oppose(*p1,*p1);
	    p1++; 
	  }
	}
	
#ifdef POLY_CH_DEBUG
	fprintf(stderr, "[Chernikova: B]\nRay = ");
	Ray->NbRows=NbRay;
	Matrix_Print(stderr,0,Ray);
	fprintf(stderr, "\nConstraints = ");
	Matrix_Print(stderr,0,Mat);
	fprintf(stderr, "\nSat = ");
	SMPrint(Sat);
#endif
	
	/* Compute the new pointed cone */
	for (i=NbBid+1; i<NbRay; i++)
	  if (value_notzero_p(Ray->p[i][0]))
	    Combine(Ray->p[i],Ray->p[NbBid],Ray->p[i],0,Dimension);
	
	/* Add the new ray */
	if (value_notzero_p(Mat->p[k][0])) { /* Constraint is an inequality */ 
	  for (j=0;j<sat_nbcolumns;j++) {
	    Sat->p[NbBid][j] = 0;     /* Saturation vec for new ray */
	  }
	  /* The new ray saturates everything except last inequality */
	  Sat->p[NbBid][jx] |= bx;
	}
	else {                        /* Constraint is an equality */
	  if (--NbRay != NbBid) {
	    Vector_Copy(Ray->p[NbRay],Ray->p[NbBid],Dimension+1);
	    SMVector_Copy(Sat->p[NbRay],Sat->p[NbBid],sat_nbcolumns);
	  }
	}

#ifdef POLY_CH_DEBUG
	fprintf(stderr, "[Chernikova: C]\nRay = ");
	Ray->NbRows=NbRay;
	Matrix_Print(stderr,0,Ray);
	fprintf(stderr, "\nConstraints = ");
	Matrix_Print(stderr,0,Mat);
	fprintf(stderr, "\nSat = ");
	SMPrint (Sat);
#endif

      } 
      else {  /* If the new constraint satisfies all the rays */
	RaySort(Ray, Sat, NbBid, NbRay, &equal_bound, &sup_bound,
		RowSize1, RowSize2,bx,jx);
        
	/* Sort the unidirectional rays into R0, R+, R- */
        /* Ray 
	   NbRay-> bound-> ________
	                   |  R-  |    R- ==> ray.eq < 0  (outside domain)
		     sup-> |------|
		           |  R+  |    R+ ==> ray.eq > 0  (inside domain)
		   equal-> |------|
		           |  R0  |    R0 ==> ray.eq = 0  (on face of domain)
		   NbBid-> |______|
        */

#ifdef POLY_CH_DEBUG
	fprintf(stderr, "[Chernikova: D]\nRay = ");
	Ray->NbRows=NbRay;
	Matrix_Print(stderr,0,Ray);
	fprintf(stderr, "\nConstraints = ");
	Matrix_Print(stderr,0,Mat);
	fprintf(stderr, "\nSat = ");
	SMPrint (Sat);
#endif

	/* Compute only the new pointed cone */
	bound=NbRay;
	for (i=equal_bound; i<sup_bound; i++) /* for all pairs of R- and R+ */
	  for(j=sup_bound; j<bound; j++) {    
	    
	  /*--------------------------------------------------------------*/
	  /* Count the set of constraints saturated by R+ and R- */
	  /* Includes equalities, inequalities and the positivity constraint */
	  /*-----------------------------------------------------------------*/
	    
	    nbcommonconstraints = 0;
	    for (l=0; l<jx; l++) {
	      aux = Temp[l] = Sat->p[i][l] | Sat->p[j][l];
	      for (m=MSB; m!=0; m>>=1) 
		if (!(aux&m)) 
		  nbcommonconstraints++;
	    }
	    aux = Temp[jx] =  Sat->p[i][jx] | Sat->p[j][jx];
	    for (m=MSB; m!=bx; m>>=1) 
	      if (!(aux&m)) 
		nbcommonconstraints++;	    
	    rayonly = (value_zero_p(Ray->p[i][Dimension])  && 
		       value_zero_p(Ray->p[j][Dimension]) && 
		       (dual == 0));	      	    
	    if(rayonly)
	      nbcommonconstraints++;      /* account for pos constr */

          /*-----------------------------------------------------------------*/
          /* Adjacency Test : is combination [R-,R+] a non redundant ray?    */
          /*-----------------------------------------------------------------*/
          
	    if (nbcommonconstraints+NbBid>=Dimension-2) { /* Dimensionality check*/
	      /* Check whether a ray m saturates the same set of constraints */
	      redundant=0;
	      for (m=NbBid; m<bound; m++) 
		if ((m!=i)&&(m!=j)) {
		  
		  /* Two rays (r+ r-) are never made redundant by a vertex */
		  /* because the positivity constraint saturates both rays */
		  /* but not the vertex                                    */
		  
		  if (rayonly && value_notzero_p(Ray->p[m][Dimension]))
		    continue;

		  /* (r+ r-) is redundant if there doesn't exist an equation */
		  /* which saturates both r+ and r- but not rm.              */
		  
		  ip1 = Temp;
		  ip2 = Sat->p[m];
		  for (l=0; l<=jx; l++,ip2++,ip1++)
		    if (*ip2 & ~*ip1) 
		      break;
		  if (l>jx) { 
		    redundant=1;
		    break;
		  }
		}

#ifdef POLY_CH_DEBUG
	      fprintf(stderr, "[Chernikova: E]\nRay = ");
	      Ray->NbRows=NbRay;
	      Matrix_Print(stderr,0,Ray);
	      fprintf(stderr, "\nConstraints = ");
	      Matrix_Print(stderr,0,Mat);
	      fprintf(stderr, "\nSat = ");
	      SMPrint (Sat);
#endif
	      
	      /*------------------------------------------------------------*/
	      /* Add new ray generated by [R+,R-]                           */
	      /*------------------------------------------------------------*/
	    
	      if (!redundant) {
		if (NbRay==NbMaxRays) {
		  NbMaxRays *= 2;
		  Ray->NbRows = NbRay;
		  Matrix_Extend(Ray, NbMaxRays);
		  SatMatrix_Extend(Sat, Mat, NbMaxRays);
		}
		
		/* Compute the new ray */
		Combine(Ray->p[j],Ray->p[i],Ray->p[NbRay],0,Dimension);
		SatVector_OR(Sat->p[j],Sat->p[i],Sat->p[NbRay],sat_nbcolumns);
		Sat->p[NbRay][jx] &= ~bx;
		NbRay++;
	      }
	    }
	  }

#ifdef POLY_CH_DEBUG
	fprintf(stderr, 
		"[Chernikova: F]\n"
		"sup_bound=%d\n"
		"equal_bound=%d\n"
		"bound=%d\n"
		"NbRay=%d\n"
		"Dimension = %d\n"
		"Ray = ",sup_bound,equal_bound,bound,NbRay,Dimension);
#endif
#ifdef POLY_CH_DEBUG
	Ray->NbRows=NbRay;
	fprintf(stderr, "[Chernikova: F]:\nRay = ");
	Matrix_Print(stderr,0,Ray);
#endif
	
	/* Eliminates all non extremal rays */
	/* j = (Mat->p[k][0]) ? */
	
	j = (value_notzero_p(Mat->p[k][0])) ? 
	  sup_bound : equal_bound;
	
	i = NbRay;
#ifdef POLY_CH_DEBUG
	fprintf(stderr, "i = %d\nj = %d \n", i, j);
#endif
	while ((j<bound)&&(i>bound)) {
	  i--;
	  Vector_Copy(Ray->p[i],Ray->p[j],Dimension+1);
	  SMVector_Copy(Sat->p[i],Sat->p[j],sat_nbcolumns);
	  j++;
	}

#ifdef POLY_CH_DEBUG
	fprintf(stderr, "i = %d\nj = %d \n", i, j);
	fprintf(stderr, 
		"[Chernikova: F]\n"
		"sup_bound=%d\n"
		"equal_bound=%d\n"
		"bound=%d\n"
		"NbRay=%d\n"
		"Dimension = %d\n"
		"Ray = ",sup_bound,equal_bound,bound,NbRay, Dimension);
#endif
#ifdef POLY_CH_DEBUG
	Ray->NbRows=NbRay;
	fprintf(stderr, "[Chernikova: G]\nRay = "); 
	Matrix_Print(stderr,0,Ray);
#endif	
	if (j==bound) 
	  NbRay=i;
	else 
	  NbRay=j;
      }
      NEXT(jx,bx);
    }    
    Ray->NbRows=NbRay;
    Sat->NbRows=NbRay;
    
  } /* End of TRY */

  UNCATCH(any_exception_error);
  free(Temp);
  
#ifdef POLY_CH_DEBUG
  fprintf(stderr, "[Chernikova: Output]\nRay = ");
  Matrix_Print(stderr,0,Ray);
  fprintf(stderr, "\nConstraints = ");
  Matrix_Print(stderr,0,Mat);
  fprintf(stderr, "\nSat = ");
  SMPrint (Sat);
#endif
  
  return 0;
} /* Chernikova */

static int Gauss4(Value **p, int NbEq, int NbRows, int Dimension)
{
  int i, j, k, pivot, Rank;
  int *column_index = NULL;
  Value gcd, *cp;

  value_init(gcd);
  column_index=(int *)malloc(Dimension * sizeof(int));
  if(!column_index) {	
    errormsg1("Gauss","outofmem","out of memory space");
    value_clear(gcd);
    return 0;
  }
  Rank=0;
  
  CATCH(any_exception_error) {
    if (column_index)
      free(column_index);
    value_clear(gcd);
    RETHROW();
  }
  TRY {
    
    for (j=1; j<=Dimension; j++) {   /* for each column (except status) */  
      for (i=Rank; i<NbEq; i++)      /* starting at diagonal, look down */
	
	/* if (Mat->p[i][j] != 0) */
	if (value_notzero_p(p[i][j])) 
	  break;                    /* Find the first non zero element  */    
      if (i!=NbEq) {                /* If a non-zero element is found?  */
	if (i!=Rank)                /* If it is found below the diagonal*/
          Vector_Exchange(p[Rank]+1,p[i]+1,Dimension);
	
	/* Normalize the pivot row by dividing it by the gcd */      
	/* gcd = Vector_Gcd(p[Rank]+1,Dimension) */
	Vector_Gcd(p[Rank]+1,Dimension,&gcd);
	
	/* if (gcd >= 2) */
	if (value_cmp_si(gcd, 2) >= 0) { 
	  cp = &p[Rank][1];
	  for (k=0; k<Dimension; k++) {
	    value_division (*cp,*cp,gcd);       /* *cp /= gcd */    
	    cp++;
	  }
	}
	
	/* if (Mat->p[Rank][j] < 0) */
	if (value_neg_p(p[Rank][j])) { 
	  cp = p[Rank]+1;	
	  for (k=0; k<Dimension; k++) { 
	    value_oppose(*cp, *cp); /* *cp *= -1 */ 
	    cp++;
	  }
	}
	/* End of normalize */
	
	pivot=i;
	for (i=pivot+1; i<NbEq; i++) {  /* Zero out the rest of the column */
	  
	  /* if (Mat->p[i][j] != 0) */
	  if (value_notzero_p(p[i][j]))
	    Combine(p[i],p[Rank],p[i],j,Dimension);
	}
	
        /* For each row with non-zero entry Mat->p[Rank], store the column */
        /* number 'j' in 'column_index[Rank]'. This information will be    */
        /* useful in performing Gaussian elimination backward step.        */
	
	column_index[Rank]=j;
	Rank++;
      }
    } /* end of Gaussian elimination forward step */
    
    /* Back Substitution -- normalize the system of equations */
    for (k=Rank-1; k>=0; k--) { 
      j = column_index[k];

      /* Normalize the equations */
      for (i=0; i<k; i++) { 
	
	/* if (Mat->p[i][j] != 0) */
	if (value_notzero_p(p[i][j]))
	  Combine(p[i],p[k],p[i],j,Dimension);
      }
      
      /* Normalize the inequalities */
      for (i=NbEq;i<NbRows;i++) { 
	
	/* if (Mat->p[i][j] != 0) */
	if (value_notzero_p(p[i][j]))
	  Combine(p[i],p[k],p[i],j,Dimension);
      }
    }
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  free(column_index), column_index = NULL;
  
  value_clear(gcd);
  return Rank;
} /* Gauss */

/*  
 * Compute a minimal system of equations using Gausian elimination method.
 * 'Mat' is a matrix of constraints in which the first 'Nbeq' constraints
 * are equations. The dimension of the homogenous system is 'Dimension'. 
 * The function returns the rank of the matrix 'Mat'. 
 */
int Gauss(Matrix *Mat, int NbEq, int Dimension)
{
  int Rank;

#ifdef POLY_DEBUG
  fprintf(stderr, "[Gauss : Input]\nRay =");
  Matrix_Print(stderr,0,Mat);
#endif

  Rank = Gauss4(Mat->p, NbEq, Mat->NbRows, Dimension);

#ifdef POLY_DEBUG
  fprintf(stderr, "[Gauss : Output]\nRay =");
  Matrix_Print(stderr,0,Mat);
#endif

  return Rank;
}

/*
 * Given 'Mat' - a matrix of equations and inequalities, 'Ray' - a matrix of 
 * lines and rays, 'Sat' - the corresponding saturation matrix, and 'Filter'
 * - an array to mark (with 1) the non-redundant equalities and inequalities, 
 * compute a polyhedron composed of 'Mat' as constraint matrix and 'Ray' as 
 * ray matrix after reductions. This function is usually called as a follow
 * up to 'Chernikova' to remove redundant constraints or rays.
 * Note: (1) 'Chernikova' ensures that there are no redundant lines and rays. 
 *       (2) The same function can be used with constraint and ray matrix used
  interchangbly.
 */ 
static Polyhedron *Remove_Redundants(Matrix *Mat,Matrix *Ray,SatMatrix *Sat,unsigned *Filter) { 
  
  int i, j, k;
  unsigned Dimension, sat_nbcolumns, NbRay, NbConstraints, RowSize1, RowSize2, 
           *Trace = NULL, *bx = NULL, *jx = NULL, Dim_RaySpace, b;
  unsigned NbBid, NbUni, NbEq, NbIneq;
  unsigned NbBid2, NbUni2, NbEq2, NbIneq2;
  int Redundant;
  int aux, *temp2 = NULL;
  Polyhedron *Pol = NULL;
  Value *temp1 = NULL;
  Value *p, *q;
  Value Status,tmp1,tmp2,tmp3;
  
  Dimension = Mat->NbColumns-1;     /* Homogeneous Dimension */
  NbRay = Ray->NbRows;
  sat_nbcolumns = Sat->NbColumns;
  NbConstraints = Mat->NbRows;
  RowSize1=(Dimension+1);
  RowSize2=sat_nbcolumns * sizeof(int);
  
  temp1=(Value *)malloc(RowSize1*sizeof(Value));
  if(!temp1) {	
    errormsg1("Remove_Redundants", "outofmem", "out of memory space");
    return 0;
  }

  /* Initialize all the 'Value' variables */
  value_init(Status); value_init(tmp1);
  value_init(tmp2); value_init(tmp3);
  
  for(i=0;i<RowSize1;i++)
    value_init(temp1[i]);

  temp2=(int *)malloc(RowSize2);
  if(!temp2) {
    errormsg1("Remove_Redundants", "outofmem", "out of memory space");
    
    /* Clear all the 'Value' variables */
    value_clear(Status); value_clear(tmp1);
    value_clear(tmp2); value_clear(tmp3);
    for(i=0;i<RowSize1;i++)
      value_clear(temp1[i]);
    free(temp1);
    return 0;
  }
  
  /* Introduce indirections into saturation matrix 'Sat' to simplify */
  /* processing with 'Sat' and allow easy exchanges of columns.      */
  bx = (unsigned *)malloc(NbConstraints * sizeof(unsigned));
  if(!bx) {
    errormsg1("Remove_Redundants", "outofmem", "out of memory space");
    
    /* Clear all the 'Value' variables */
    value_clear(Status); value_clear(tmp1);
    value_clear(tmp2); value_clear(tmp3);
    for(i=0;i<RowSize1;i++)
      value_clear(temp1[i]);
    free(temp1); free(temp2);
    return 0;
  }
  jx = (unsigned *)malloc(NbConstraints * sizeof(unsigned));
  if(!jx) {
    errormsg1("Remove_Redundants", "outofmem", "out of memory space");
    
    /* Clear all the 'Value' variables */
    value_clear(Status); value_clear(tmp1);
    value_clear(tmp2); value_clear(tmp3);
    for(i=0;i<RowSize1;i++)
      value_clear(temp1[i]);
    free(temp1); free(temp2); free(bx);
    return 0;
  }
  CATCH(any_exception_error) {  
    
    if (temp1) {           
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      free(temp1);
    }  
    if (temp2) free(temp2);
    if (bx) free(bx);
    if (jx) free(jx);
    if (Trace) free(Trace);
    if (Pol) Polyhedron_Free(Pol);

    /* Clear all the 'Value' variables */
    value_clear(Status); value_clear(tmp1);
    value_clear(tmp2); value_clear(tmp3);
    
    RETHROW();
  }
  TRY {
    
    /* For each constraint 'j' following mapping is defined to facilitate  */
    /* data access from saturation matrix 'Sat' :-                         */
    /* (1) jx[j] -> floor[j/(8*sizeof(int))]                               */
    /* (2) bx[j] -> bin(00..10..0) where position of 1 = j%(8*sizeof(int)) */
    
    i = 0;
    b = MSB;
    for (j=0; j<NbConstraints; j++) { 
      jx[j] = i;
      bx[j] = b;
      NEXT(i,b);
    }
    
    /* 
     * STEP(0): Count the number of vertices among the rays while initializing
     * the ray status count to 0. If no vertices are found, quit the procedure
     * and return an empty polyhedron as the result. 
     */            
    
    /* Reset the status element of each ray to zero. Store the number of  */
    /* vertices in 'aux'.                                                 */
    aux = 0;
    for (i=0; i<NbRay; i++) {  
      
      /* Ray->p[i][0] = 0 */
      value_set_si(Ray->p[i][0],0);
      
      /* If ray(i) is a vertex of the Inhomogenous system */
      if (value_notzero_p(Ray->p[i][Dimension]))
	aux++;              
    }
    
    /* If no vertices, return an empty polyhedron. */
    if (!aux) { 
      
      /* Clear all the 'Value' variables */
      value_clear(Status); value_clear(tmp1);
      value_clear(tmp2); value_clear(tmp3);
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      
      /* Return an empty polyhedron */
      free(temp1); free(temp2); free(jx); free(bx);
      UNCATCH(any_exception_error);
      return Empty_Polyhedron(Dimension-1);
    }
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Init]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRays =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(1): Compute status counts for both rays and inequalities. For each
     * constraint, count the number of vertices/rays saturated by that 
     * constraint, and put the result in the status words. At the same time, 
     * for each vertex/ray, count the number of constraints saturated by it.
     * Delete any positivity constraints, but give rays credit in their status
     * counts for saturating the positivity constraint.
     */
    
    NbEq=0;
 #ifdef JUNK
    /* JUNK is a temporary flag, the code in the JUNK part, should probably
       be removed (see with Fabien and Doran) */
   memset((char *)temp2, 0, RowSize2);       
#endif
    
#ifdef POLY_RR_DEBUG
    fprintf (stderr, " j = ");
#endif
    
    for (j=0; j<NbConstraints; j++) {
      
#ifdef POLY_RR_DEBUG
      fprintf (stderr, " %i ", j);
      fflush (stderr);
#endif
      
#ifdef JUNK
      /* If constraint(j) is an equality, mark '1' in array 'temp2' */
      if (Filter && value_zero_p(Mat->p[j][0]))  
	temp2[jx[j]] |= bx[j]; 
#endif
      /* Reset the status element of each constraint to zero */
      value_set_si(Mat->p[j][0],0);
      
      /* Identify and remove the positivity constraint 1>=0 */
      for (i=1, p = &Mat->p[j][1]; i<Dimension; i++) { 
	
	/* if (*p) */
	if (value_notzero_p(*p)) {
	  p++; 
	  break;
	}
	else 
	  p++;
      }
      
#ifdef POLY_RR_DEBUG
      fprintf(stderr, "[Remove_redundants : IntoStep1]\nConstraints =");
      Matrix_Print(stderr,0,Mat);
      fprintf (stderr, " j = %i \n", j);
#endif
      
      /* Check if constraint(j) is a positivity constraint, 1 >= 0, or if it */
      /* is 1==0. If constraint(j) saturates all the rays of the matrix 'Ray'*/
      /* then it is an equality. in this case, return an empty polyhedron.   */
      
      if (i==Dimension) { 
	for (i=0; i<NbRay; i++)
	  if (!(Sat->p[i][jx[j]]&bx[j])) {
	    
	    /* Mat->p[j][0]++ */
	    value_increment(Mat->p[j][0],Mat->p[j][0]);
	  }
	
        /* if ((Mat->p[j][0] == NbRay) &&   : it is an equality
	   (Mat->p[j][Dimension] != 0)) : and its not 0=0 */
	value_set_si(tmp1,NbRay);
        if ((value_eq(Mat->p[j][0],tmp1)) &&  
            (value_notzero_p(Mat->p[j][Dimension]))) {
	  
	  /* Clear all the 'Value' variables */
	  value_clear(Status); value_clear(tmp1);
	  value_clear(tmp2); value_clear(tmp3);
	  for(i=0;i<RowSize1;i++)
	    value_clear(temp1[i]);
	  
	  /* Return an empty polyhedron */
	  free(temp1); free(temp2); free(jx); free(bx);
	  UNCATCH(any_exception_error);
	  return Empty_Polyhedron(Dimension-1);
        }
	
        /* Delete the positivity constraint */
        NbConstraints--;
        if (j==NbConstraints) continue;
        Vector_Exchange(Mat->p[j], Mat->p[NbConstraints],RowSize1);
        exchange(jx[j], jx[NbConstraints], aux);
        exchange(bx[j], bx[NbConstraints], aux);
        j--; continue;
      }
      
      /* Count the number of vertices/rays saturated by each constraint. At  */
      /* the same time, count the number of constraints saturated by each ray*/
      for (i=0; i<NbRay; i++) 
	if (!(Sat->p[i][jx[j]]&bx[j])) {  
	  
	  /* Mat->p[j][0]++ */
	  value_increment(Mat->p[j][0],Mat->p[j][0]);
	  
	  /* Ray->p[i][0]++ */
	  value_increment (Ray->p[i][0],Ray->p[i][0]);
	}
      
      /* if (Mat->p[j][0]==NbRay) then increment the number of eq. count */
      value_set_si(tmp1,NbRay);
      if (value_eq(Mat->p[j][0],tmp1)) 
	NbEq++;    /* all vertices/rays are saturated */
    }
    Mat->NbRows = NbConstraints;
    
    NbBid=0;
    for (i=0; i<NbRay; i++) {
      
      /* Give rays credit for saturating the positivity constraint */
      if (value_zero_p(Ray->p[i][Dimension]))	
	
	/* Ray->p[i][0]++ */
	value_increment(Ray->p[i][0],Ray->p[i][0]);
      
      /* If ray(i) saturates all the constraints including positivity  */
      /* constraint then it is a bi-directional ray or line. Increment */
      /* 'NbBid' by one.                                               */
      
      /* if (Ray->p[i][0]==NbConstraints+1) */
      value_set_si(tmp1,(NbConstraints+1));
      if (value_eq(Ray->p[i][0],tmp1))
	NbBid++;
    }
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step1]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(2): Sort equalities to the top of constraint matrix 'Mat'. Detect
     * implicit equations such as y>=3; y<=3. Keep Inequalities in same 
     * relative order. (Note: Equalities are constraints which saturate all of
     * the rays)              
     */
    
    for (i=0; i<NbEq; i++) {
      
      /* If constraint(i) doesn't saturate some ray, then it is an inequality*/
      value_set_si(tmp1,NbRay);
      if (value_ne(Mat->p[i][0],tmp1)) { 
	
	value_set_si(tmp1,NbRay);
	/* Skip over inequalities and find an equality */
	for (k=i+1;value_ne(Mat->p[k][0],tmp1) && k<NbConstraints;k++);
	if (k==NbConstraints) /* If none found then error */ break;
	
	/* Slide inequalities down the array 'Mat' and move equality up to */
	/* position 'i'.                                                   */
	Vector_Copy(Mat->p[k], temp1,RowSize1);
	aux = jx[k];
	j   = bx[k];
	for (;k>i;k--) {  
	  Vector_Copy(Mat->p[k-1],Mat->p[k],RowSize1);
	  jx[k] = jx[k-1];
	  bx[k] = bx[k-1];
	}
	Vector_Copy(temp1,Mat->p[i],RowSize1);
	jx[i] = aux;
	bx[i] = j;
      }
    }

#ifdef JUNK
    if (Filter)                     /* for SIMPLIFY */
      for (i=0; i<NbEq; i++) {
	
	/* Detect implicit constraints such as y>=3 and y<=3 */
	Redundant = 0;
	for (j=i+1; j<NbEq; j++) {
	  for (k=0, p=&Mat->p[i][1], q=&Mat->p[j][1]; k<Dimension; k++,p++,q++) {  
	    /* if (*p!=*q) */
	    if (value_ne(*p, *q)) 
	      break;
	  }
	  
	  /* Redundant if both are same `and' constraint(j) was equality. */
	  /* That is, 'temp2' has entry 1                                 */
	  if (k==Dimension && (temp2[jx[j]] & bx[j])) { 
	    Redundant=1; 
	    break;
	  }
	}
	
	/* Set 'Filter' entry to 1 corresponding to the irredundant equality*/
	if (!Redundant) Filter[jx[i]] |= bx[i];  /* set flag */
      }

#endif
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step2]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(3): Perform Gaussian elimiation on the list of equalities. Obtain
     * a minimal basis by solving for as many variables as possible. Use this 
     * solution to reduce the inequalities by eliminating as many variables as
     * possible. Set NbEq2 to the rank of the system of equalities.
     */
    
    NbEq2 = Gauss(Mat,NbEq,Dimension);
    
    /* If number of equalities is not less then the homogenous dimension, */
    /* return an empty polyhedron.                                        */
    
    if (NbEq2>=Dimension) {		
      
      /* Clear all the 'Value' variables */
      value_clear(Status); value_clear(tmp1);
      value_clear(tmp2); value_clear(tmp3);
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      
      free(temp1); free(temp2); free(jx); free(bx);
      UNCATCH(any_exception_error);
      return Empty_Polyhedron(Dimension-1);
    }
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step3]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /*
     * STEP(4): Sort lines to the top of ray matrix 'Ray', leaving rays
     * afterwards. Detect implicit lines such as ray(1,2) and ray(-1,-2). 
     * (Note: Lines are rays which saturate all of the constraints including
     * the positivity constraint 1>=0. 
     */
    
    
    for (i=0, k=NbRay; i<NbBid && k>i; i++) {
      value_set_si(tmp1,(NbConstraints+1));
      
      /* If ray(i) doesn't saturate some constraint then it is not a line */
      if (value_ne(Ray->p[i][0],tmp1)) { 
	
	value_set_si(tmp1,(NbConstraints+1));
	/* Skip over rays and vertices and find a line (bi-directional rays) */
	while (--k >i && value_ne(Ray->p[k][0],tmp1)) ;
	
	/* Exchange positions of ray(i) and line(k), thus sorting lines to */
	/* the top of matrix 'Ray'.                                        */
	Vector_Exchange(Ray->p[i], Ray->p[k], RowSize1);
	bexchange(Sat->p[i], Sat->p[k], RowSize2);
      }     
    }	

#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step4]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(5): Perform Gaussian elimination on the lineality space to obtain
     * a minimal basis of lines. Use this basis to reduce the representation
     * of the uniderectional rays. Set 'NbBid2' to the rank of the system of 
     * lines. 
     */
    
    NbBid2 = Gauss(Ray, NbBid, Dimension);
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : After Gauss]\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* If number of lines is not less then the homogenous dimension, return */
    /* an empty polyhedron.                                                 */
    if (NbBid2>=Dimension) {
      errormsg1("RemoveRedundants", "rmrdt", "dimension error");
      
      /* Clear all the 'Value' variables */
      value_clear(Status); value_clear(tmp1);
      value_clear(tmp2); value_clear(tmp3);
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      
      free(temp1); free(temp2); free(jx); free(bx);
      UNCATCH(any_exception_error);
      return Empty_Polyhedron(Dimension-1);
    }
    
    /* Compute dimension of non-homogenous ray space */
    Dim_RaySpace = Dimension-1-NbEq2-NbBid2;  
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step5]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(6): Do a first pass filter of inequalities and equality identifi-
     * cation. New positivity constraints may have been created by step(3). 
     * Check for and eliminate them. Count the irredundant inequalities and 
     * store count in 'NbIneq'.  
     */
 
    value_set_si(tmp2,Dim_RaySpace);
    value_set_si(tmp3,NbRay);
    NbIneq=0;
    for (j=0; j<NbConstraints; j++) {
      
      /* Identify and remove the positivity constraint 1>=0 */
      for (i=1, p = &Mat->p[j][1]; i<Dimension; i++)
	if (value_notzero_p (*p)) {
	  p++; 
	  break;
	}
	else
	  p++;
      
      /* Check if constraint(j) is a positivity constraint, 1>= 0, or if it */
      /* is 1==0.                                                           */
      if (i==Dimension) {  
	
	value_set_si(tmp1,NbRay);
	
        /* if ((Mat->p[j][0]==NbRay) &&   : it is an equality 
	   (Mat->p[j][Dimension]!=0))    : and its not 0=0 */
        if ((value_eq (Mat->p[j][0],tmp1)) &&
            (value_notzero_p(Mat->p[j][Dimension]))) {
	  
	  /* Clear all the 'Value' variables */
	  value_clear(Status); value_clear(tmp1);
	  value_clear(tmp2); value_clear(tmp3);
	  for(i=0;i<RowSize1;i++)
	    value_clear(temp1[i]);
	  
	  /* Return an empty polyhedron */
	  free(temp1); free(temp2); free(jx); free(bx);
	  UNCATCH(any_exception_error);
	  return Empty_Polyhedron(Dimension-1);
        }	
	
        /* Set the positivity constraint redundant by setting status element */
        /* equal to 2.                                                       */
	value_set_si(Mat->p[j][0],2);
        continue;
      }
      
      /* Status = Mat->p[j][0] */
      value_assign(Status, Mat->p[j][0]);
      
      /* if (Status == 0) then constraint is redundant */
      if (value_zero_p(Status)) 	
	
	/* Mat->p[j][0]=2 : redundant */
	value_set_si(Mat->p[j][0],2);	
      
      /* else if (Status<Dim_RaySpace) then constraint is redundant */   
      else if (value_lt(Status,tmp2)) 	
	
	/* Mat->p[j][0]=2 : redundant */
	value_set_si(Mat->p[j][0],2);	
      
      /* else if (Status==NbRay) then constraint is an equality */
      else if (value_eq(Status,tmp3)) 	
	
	/* Mat->p[j][0]=0 : equality */
	value_set_si(Mat->p[j][0],0);
      
      /* else constraint is an irredundant inequality */ 
      else { 
	NbIneq++; 	
	
	/* Mat->p[j][0]=1 : inequality */
	value_set_si(Mat->p[j][0],1);
      }
    }
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step6]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(7): Do a first pass filter of rays and identification of lines.
     * Count the irredundant Rays and store count in 'NbUni'. 
     */
    
    value_set_si(tmp2,Dim_RaySpace);
    value_set_si(tmp3,(NbConstraints+1));
    NbUni=0;
    for (j=0; j<NbRay; j++) { 
      
      /* Status = Ray->p[j][0] */
      value_assign(Status, Ray->p[j][0]);
      
      /* if (Status < Dim_RaySpace) the ray is redundant */
      if (value_lt(Status,tmp2)) 	
	
	/* Ray->p[j][0]=2 : redundant */
	value_set_si(Ray->p[j][0],2);	
      
      /* else if (Status == (NbConstraints+1)) then ray is a line */
      else if (value_eq(Status,tmp3)) 
	
	/* Ray->p[j][0]=0 : line */
	value_set_si(Ray->p[j][0],0);
      
      /* else ray is an irredundant unidirectional ray. */
      else {
	NbUni++; 
	
	/* Ray->p[j][0]=1 : ray */
	value_set_si(Ray->p[j][0],1);
      }
    }
    
    /*
     * STEP(8): Create the polyhedron (using approximate sizes).
     * Number of constraints = NbIneq + NbEq2 + 1
     * Number of rays = NbUni + NbBid2 
     * Partially fill the polyhedron structure with the lines computed in step
     * 3 and the equalities computed in step 5. 
     */
    
    Pol = Polyhedron_Alloc(Dimension-1, NbIneq+NbEq2+1, NbUni+NbBid2);
    if (!Pol) {
      errormsg1("Remove_redundants", "outofmem", "out of memory space");
      
      /* Clear all the 'Value' variables */
      value_clear(Status); value_clear(tmp1);
      value_clear(tmp2); value_clear(tmp3);
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      
      free(temp1); free(temp2); free(jx); free(bx);
      UNCATCH(any_exception_error);
      return 0;
    }
    Pol->NbBid = NbBid2;
    Pol->NbEq  = NbEq2;
    
    /* Partially fill the polyhedron structure */
    if (NbBid2) Vector_Copy(Ray->p[0], Pol->Ray[0], (Dimension+1)*NbBid2);
    if (NbEq2)  Vector_Copy(Mat->p[0], Pol->Constraint[0], (Dimension+1)*NbEq2);
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step7]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * STEP(9): Final Pass filter of inequalities and detection of redundant
     * inequalities. Redundant inequalities include: 
     * (1) Inequalities which are always true, such as 1>=0, 
     * (2) Redundant inequalities such as y>=4 given y>=3, or x>=1 given x=2. 
     * (3) Redundant inequalities such as x+y>=5 given x>=3 and y>=2.
     * Every 'good' inequality must saturate at least 'Dimension' rays and be 
     * unique.
     */
    
    /* 'Trace' is a (1 X sat_nbcolumns) row matrix to hold the union of all */
    /* rows (corresponding to irredundant rays) of saturation matrix 'Sat'  */
    /* which saturate some constraint 'j'. See figure below:-               */
    Trace=(unsigned *)malloc(sat_nbcolumns * sizeof(unsigned));
    if(!Trace) {
      errormsg1("Remove_Redundants", "outofmem", "out of memory space");
      
      /* Clear all the 'Value' variables */
      value_clear(Status); value_clear(tmp1);
      value_clear(tmp2); value_clear(tmp3);
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      
      free(temp1); free(temp2); free(jx); free(bx);
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /*                         NbEq      NbConstraints
			       |----------->
			        ___________j____
			       |           |   |
			       |      Mat  |   |
			       |___________|___|
		                           |                  
     NbRay  ^ ________         ____________|____
	    | |-------|--------|-----------0---|t1
	    |i|-------|--------|-----------0---|t2
	    | | Ray   |        |    Sat        |
     NbBid  - |-------|--------|-----------0---|tk
	      |_______|        |_______________|
			                   |
			                   |
			              -OR- (of rows t1,t2,...,tk)
			       ________|___|____
			       |_____Trace_0___|
			       
    */
    
    NbIneq2 = 0;
    for (j=NbEq; j<NbConstraints; j++) {
      
      /* if (Mat->p[j][0]==1) : non-redundant inequality */
      if (value_one_p (Mat->p[j][0])) { 
	for (k=0; k<sat_nbcolumns; k++) Trace[k]=0;  /* init Trace */
	
	/* Compute Trace: the union of all rows of Sat where constraint(j) */
	/* is saturated.                                                   */
	for (i=NbBid; i<NbRay; i++) 
	  
	  /* if (Ray->p[i][0]==1) */
	  if (value_one_p(Ray->p[i][0])) { 
	    if (!(Sat->p[i][jx[j]]&bx[j])) 
	      for (k=0; k<sat_nbcolumns; k++) Trace[k] |= Sat->p[i][k];
	  }
	
	/* Only constraint(j) should saturate this set of vertices/rays. */
	/* If another non-redundant constraint also saturates this set,  */
	/* then constraint(j) is redundant                               */
	Redundant=0;
	for (i=NbEq; i<NbConstraints; i++) {
	  
	  /* if ((Mat->p[i][0] ==1) && (i!=j) && !(Trace[jx[i]] & bx[i]) ) */
	  if (value_one_p(Mat->p[i][0]) && (i!=j) && !(Trace[jx[i]] & bx[i])) {
	    Redundant=1;
	    break;
	  }
	}
	if (Redundant) {
	  value_set_si(Mat->p[j][0],2);
	}	
	else {
	  Vector_Copy(Mat->p[j], Pol->Constraint[NbEq2+NbIneq2], Dimension+1);
	  if (Filter) Filter[jx[j]] |= bx[j];		/* for SIMPLIFY */
	  NbIneq2++;
	}
      }
    }
    free(Trace), Trace = NULL;
    
#ifdef POLY_RR_DEBUG
    fprintf(stderr, "[Remove_redundants : Step8]\nConstraints =");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nRay =");
    Matrix_Print(stderr,0,Ray);
#endif
    
    /* 
     * Step(10): Final pass filter of rays and detection of redundant rays.
     * The final list of rays is written to polyhedron.                     
     */
    
    /* Trace is a (NbRay x 1) column matrix to hold the union of all columns */
    /* (corresponding to irredundant inequalities) of saturation matrix 'Sat'*/
    /* which saturate some ray 'i'. See figure below:-                       */
    
    Trace=(unsigned *)malloc(NbRay * sizeof(unsigned));
    if(!Trace) {
      errormsg1("Remove_Redundants", "outofmem", "out of memory space");
      
      /* Clear all the 'Value' variables */
      value_clear(Status); value_clear(tmp1);
      value_clear(tmp2); value_clear(tmp3);
      for(i=0;i<RowSize1;i++)
	value_clear(temp1[i]);
      
      free(bx); free(jx); free(temp2); free(temp1);
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /* 			   NbEq     NbConstraints
                             |---------->
			___________j_____
			|      | |   |	|
			|      Mat   |	|
			|______|_|___|__|
                               | |   |
NbRay ^	_________	_______|_|___|__   ___
      |	|	|	|      | |   |	|  |T|
      | |  Ray	|	|   Sat| |   |	|  |r|
      | |	|	|      | |   |	|  |a|  Trace = Union[col(t1,t2,..,tk)]
      |i|-------|------>i      0 0   0  |  |c|
NbBid -	|	|       |      | |   |  |  |e|
	|_______|	|______|_|___|__|  |_|
                              t1 t2  tk
    */    
  
    NbUni2 = 0;
    
    /* Let 'aux' be the number of rays not vertices */ 
    aux = 0;     
    for (i=NbBid; i<NbRay; i++) {
      
      /* if (Ray->p[i][0]==1) */
      if (value_one_p (Ray->p[i][0])) { 
	
	/* if (Ray->p[i][Dimension]!=0) : vertex */
	if (value_notzero_p (Ray->p[i][Dimension]))
	  for (k=NbBid; k<NbRay; k++) Trace[k]=0;	  /* init Trace */
	else /* for ray */
	  
	  /* Include the positivity constraint incidences for rays. The */
	  /* positivity constraint saturates all rays and no vertices   */
	  
	  for (k=NbBid; k<NbRay; k++)
	    
	    /* Trace[k]=(Ray->p[k][Dimension]!=0); */
	    Trace[k] = (value_notzero_p (Ray->p[k][Dimension]));
	
	/* Compute Trace: the union of all columns of Sat where ray(i) is  */
	/* saturated.                                                      */
	for (j=NbEq; j<NbConstraints; j++)
	  
	  /* if (Mat->p[j][0]==1) : inequality */
	  if (value_one_p (Mat->p[j][0])) { 
	    if (!(Sat->p[i][jx[j]]&bx[j]))
	      for (k=NbBid; k<NbRay; k++) Trace[k] |= Sat->p[k][jx[j]]&bx[j];
	  }
	
	/* If ray i does not saturate any inequalities (other than the   */
	/* the positivity constraint, then it is the case that there is  */
	/* only one inequality and that ray is its orthogonal            */
	
	/* only ray(i) should saturate this set of inequalities. If      */
	/* another non-redundant ray also saturates this set, then ray(i)*/
	/* is redundant                                                  */
	
	Redundant = 0;
	for (j=NbBid; j<NbRay; j++) { 
	  
	  /* if ( (Ray->p[j][0]==1) && (i!=j) && !Trace[j] ) */
	  if (value_one_p (Ray->p[j][0]) && (i!=j) && !Trace[j]) { 
	    Redundant=1;
	    break;
	  }
	}
	if (Redundant) 
	  value_set_si(Ray->p[i][0],2);	
	else {
	  Vector_Copy(Ray->p[i], Pol->Ray[NbBid2+NbUni2], Dimension+1);
	  NbUni2++;  /* Increment number of uni-directional rays */
	  
	  /* if (Ray->p[i][Dimension]==0) */ 
	  if (value_zero_p (Ray->p[i][Dimension]))
	    aux++; /* Increment number of rays which are not vertices */
	}
      }
    }
    
    /* Include the positivity constraint */
    if (aux>=Dim_RaySpace) {
      Vector_Set(Pol->Constraint[NbEq2+NbIneq2],0,Dimension+1);
      value_set_si(Pol->Constraint[NbEq2+NbIneq2][0],1);
      value_set_si(Pol->Constraint[NbEq2+NbIneq2][Dimension],1);
      NbIneq2++;
    }   
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  
#ifdef POLY_RR_DEBUG
  fprintf(stderr, "[Remove_redundants : Step9]\nConstraints =");
  Matrix_Print(stderr,0,Mat);
  fprintf(stderr, "\nRay =");
  Matrix_Print(stderr,0,Ray);
#endif
  
  free(Trace);
  free(bx);
  free(jx);
  free(temp2);
  
  Pol->NbConstraints = NbEq2 + NbIneq2;
  Pol->NbRays = NbBid2 + NbUni2;
  
  /* Clear all the 'Value' variables */
  value_clear(Status); value_clear(tmp1);
  value_clear(tmp2); value_clear(tmp3);
  for(i=0;i<RowSize1;i++)
    value_clear(temp1[i]);
  free(temp1);
  F_SET(Pol, 
        POL_VALID | POL_INEQUALITIES | POL_FACETS | POL_POINTS | POL_VERTICES);
  return Pol;
} /* Remove_Redundants */

/*
 * Allocate memory space for polyhedron. 
 */
Polyhedron* Polyhedron_Alloc(unsigned Dimension,unsigned NbConstraints,unsigned NbRays) { 
  
  Polyhedron *Pol;
  unsigned NbRows,NbColumns;
  int i,j;
  Value *p, **q; 

  Pol=(Polyhedron *)malloc(sizeof(Polyhedron));
  if(!Pol) {
    errormsg1("Polyhedron_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  
  Pol->next          = (Polyhedron *)0;
  Pol->Dimension     = Dimension;
  Pol->NbConstraints = NbConstraints;
  Pol->NbRays        = NbRays;
  Pol->NbEq          = 0;
  Pol->NbBid         = 0;
  Pol->flags	     = 0;
  NbRows             = NbConstraints + NbRays;
  NbColumns          = Dimension + 2;
  
  q = (Value **)malloc(NbRows * sizeof(Value *));
  if(!q) {
    errormsg1("Polyhedron_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  p = value_alloc(NbRows*NbColumns, &Pol->p_Init_size);
  if(!p) {
    free(q);
    errormsg1("Polyhedron_Alloc", "outofmem", "out of memory space");
    return 0;
  }
  Pol->Constraint    = q;
  Pol->Ray           = q + NbConstraints;
  Pol->p_Init        = p;
  for (i=0;i<NbRows;i++) {
    *q++ = p;
    p += NbColumns;
  }
  return Pol;
} /* Polyhedron_Alloc */

/*    
 * Free the memory space occupied by the single polyhedron.
 */
void Polyhedron_Free(Polyhedron *Pol)
{
  if(!Pol)
    return;
  value_free(Pol->p_Init, Pol->p_Init_size);
  free(Pol->Constraint);
  free(Pol);
  return;
} /* Polyhedron_Free */

/*
 * Free the memory space occupied by the domain. 
 */
void Domain_Free(Polyhedron *Pol)
{
  Polyhedron *Next;
  
  for (; Pol; Pol = Next) {
    Next = Pol->next;
    Polyhedron_Free(Pol);
  }
  return;
} /* Domain_Free */

/*
 * Print the contents of a polyhedron. 
 */
void Polyhedron_Print(FILE *Dst,char *Format,Polyhedron *Pol) { 

  unsigned Dimension, NbConstraints, NbRays;
  int      i, j;
  Value    *p;
  
  if (!Pol) { 
    fprintf(Dst, "<null polyhedron>\n");
    return;
  }
  
  Dimension     = Pol->Dimension + 2;  /* Homogenous Dimension + status */
  NbConstraints = Pol->NbConstraints;
  NbRays        = Pol->NbRays;
  fprintf(Dst, "POLYHEDRON Dimension:%d\n", Pol->Dimension);
  fprintf(Dst,"           Constraints:%d  Equations:%d  Rays:%d  Lines:%d\n",
	  Pol->NbConstraints, Pol->NbEq, Pol->NbRays, Pol->NbBid);
  fprintf(Dst,"Constraints %d %d\n", NbConstraints, Dimension);
  
  for (i=0;i<NbConstraints;i++) {
    p=Pol->Constraint[i];
    
    /* if (*p) */
    if (value_notzero_p (*p))
      fprintf(Dst,"Inequality: [");
    else      
      fprintf(Dst,"Equality:   [");
    p++;
    for (j=1;j<Dimension;j++) {
      value_print(Dst,Format,*p++);
    }  
    (void)fprintf(Dst," ]\n");
  }

  (void)fprintf(Dst, "Rays %d %d\n", NbRays, Dimension);
  for (i=0;i<NbRays;i++) {
    p=Pol->Ray[i];
    
    /* if (*p) */
    if (value_notzero_p (*p)) {   
      p++;
      
      /* if ( p[Dimension-2] ) */
      if (value_notzero_p (p[Dimension-2]))
	fprintf(Dst, "Vertex: [");
      else                  
	fprintf(Dst, "Ray:    [");
    }
    else {
      p++;
      fprintf(Dst, "Line:   [");
    }
    for (j=1; j < Dimension-1; j++) {
      value_print(Dst,Format,*p++);  
    }  
    
    /* if (*p) */
    if (value_notzero_p (*p)) {
      fprintf( Dst, " ]/" );
      value_print(Dst,VALUE_FMT,*p);
      fprintf( Dst, "\n" );
    }
    else    
      fprintf(Dst, " ]\n");
  }
  if (Pol->next) {  
    fprintf(Dst, "UNION ");
    Polyhedron_Print(Dst,Format,Pol->next);
  }
} /* Polyhedron_Print */

/* 
 * Print the contents of a polyhedron 'Pol' (used for debugging purpose).
 */
void PolyPrint (Polyhedron *Pol) {
  Polyhedron_Print(stderr,"%4d",Pol);
} /* PolyPrint */

/* 
 * Create and return an empty polyhedron of non-homogenous dimension 
 * 'Dimension'. An empty polyhedron is characterized by :-
 *  (a) The dimension of the ray-space is -1.  
 *  (b) There is an over-constrained system of equations given by:
 *      x=0, y=0, ...... z=0, 1=0
 */
Polyhedron *Empty_Polyhedron(unsigned Dimension) {

  Polyhedron *Pol;
  int i;

  Pol = Polyhedron_Alloc(Dimension, Dimension+1, 0);
  if (!Pol) {
    errormsg1("Empty_Polyhedron", "outofmem", "out of memory space");
    return 0;
  }
  Vector_Set(Pol->Constraint[0],0,(Dimension+1)*(Dimension+2));
  for (i=0; i<=Dimension; i++) {
    
    /* Pol->Constraint[i][i+1]=1 */
    value_set_si(Pol->Constraint[i][i+1],1);
  }
  Pol->NbEq = Dimension+1;
  Pol->NbBid = 0;
  F_SET(Pol, 
        POL_VALID | POL_INEQUALITIES | POL_FACETS | POL_POINTS | POL_VERTICES);
  return Pol;
} /* Empty_Polyhedron */

/* 
 * Create and return a universe polyhedron of non-homogenous dimension
 * 'Dimension'. A universe polyhedron is characterized by :-
 * (a) The dimension of rayspace is zero. 
 * (b) The dimension of lineality space is the dimension of the polyhedron.
 * (c) There is only one constraint (positivity constraint) in the constraint
 *     set given by : 1 >= 0. 
 * (d) The bi-directional ray set is the canonical set of vectors. 
 * (e) The only vertex is the origin (0,0,0,....0).  
 */
Polyhedron *Universe_Polyhedron(unsigned Dimension) { 
  
  Polyhedron *Pol;
  int i;
  
  Pol = Polyhedron_Alloc(Dimension,1,Dimension+1);
  if (!Pol) {
    errormsg1("Universe_Polyhedron", "outofmem", "out of memory space");
    return 0;
  }
  Vector_Set(Pol->Constraint[0],0,(Dimension+2));
  
  /* Pol->Constraint[0][0] = 1 */
  value_set_si(Pol->Constraint[0][0],1);
  
  /* Pol->Constraint[0][Dimension+1] = 1 */
  value_set_si(Pol->Constraint[0][Dimension+1],1);
  Vector_Set(Pol->Ray[0],0,(Dimension+1)*(Dimension+2));
  for (i=0;i<=Dimension;i++) {
    
    /* Pol->Ray[i][i+1]=1 */
    value_set_si(Pol->Ray[i][i+1],1);
  }  
  
  /* Pol->Ray[Dimension][0] = 1 :  vertex status */
  value_set_si(Pol->Ray[Dimension][0],1);
  Pol->NbEq = 0;
  Pol->NbBid = Dimension;
  F_SET(Pol, 
        POL_VALID | POL_INEQUALITIES | POL_FACETS | POL_POINTS | POL_VERTICES);
  return Pol;
} /* Universe_Polyhedron */

/*

Sort constraints and remove trivially redundant constraints.

*/
static void SortConstraints(Matrix *Constraints, unsigned NbEq)
{
    int i, j, k;
    for (i = NbEq; i < Constraints->NbRows; ++i) {
	int max = i;
	for (k = i+1; k < Constraints->NbRows; ++k) {
	    for (j = 1; j < Constraints->NbColumns-1; ++j) {
		if (value_eq(Constraints->p[k][j],
			     Constraints->p[max][j]))
		    continue;
		if (value_abs_lt(Constraints->p[k][j],
				 Constraints->p[max][j]))
		    break;
		if (value_abs_eq(Constraints->p[k][j],
				 Constraints->p[max][j]) &&
		    value_pos_p(Constraints->p[max][j]))
			break;
		max = k;
		break;
	    }
	    /* equal, except for possibly the constant
	     * => remove constraint with biggest constant
	     */
	    if (j == Constraints->NbColumns-1) {
		if (value_lt(Constraints->p[k][j], Constraints->p[max][j]))
		    Vector_Exchange(Constraints->p[k], 
				    Constraints->p[max], 
				    Constraints->NbColumns);
		Constraints->NbRows--;
		if (k < Constraints->NbRows)
		    Vector_Exchange(Constraints->p[k], 
				    Constraints->p[Constraints->NbRows], 
				    Constraints->NbColumns);
		k--;
	    }
	}
	if (max != i)
	    Vector_Exchange(Constraints->p[max], Constraints->p[i], 
			    Constraints->NbColumns);
    }
}

/*

Search for trivial implicit equalities,
assuming the constraints have been sorted.

*/

static int ImplicitEqualities(Matrix *Constraints, unsigned NbEq)
{
    int row, nrow, k;
    int found = 0;
    Value tmp;
    for (row = NbEq; row < Constraints->NbRows; ++row) {
	int d = First_Non_Zero(Constraints->p[row]+1, Constraints->NbColumns-2);
	if (d == -1) {
	    /* -n >= 0 */
	    if (value_neg_p(Constraints->p[row][Constraints->NbColumns-1])) {
		value_set_si(Constraints->p[row][0], 0);
		found = 1;
	    }
	    break;
	}
	if (value_neg_p(Constraints->p[row][1+d]))
	    continue;
	for (nrow = row+1; nrow < Constraints->NbRows; ++nrow) {
	    if (value_zero_p(Constraints->p[nrow][1+d]))
		break;
	    if (value_pos_p(Constraints->p[nrow][1+d]))
		continue;
	    for (k = d; k < Constraints->NbColumns-1; ++k) {
		if (value_abs_ne(Constraints->p[row][1+k], 
				 Constraints->p[nrow][1+k]))
		    break;
		if (value_zero_p(Constraints->p[row][1+k]))
		    continue;
		if (value_eq(Constraints->p[row][1+k], Constraints->p[nrow][1+k]))
		    break;
	    }
	    if (k == Constraints->NbColumns-1) {
		value_set_si(Constraints->p[row][0], 0);
		value_set_si(Constraints->p[nrow][0], 0);
		found = 1;
		break;
	    }
	    if (k != Constraints->NbColumns-2)
		continue;
	    /* if the constants are such that 
	     * the sum c1+c2 is negative then the constraints conflict
	     */
	    value_init(tmp);
	    value_addto(tmp, Constraints->p[row][1+k], 
			     Constraints->p[nrow][1+k]);
	    if (value_sign(tmp) < 0) {
		Vector_Set(Constraints->p[row], 0, Constraints->NbColumns-1);
		Vector_Set(Constraints->p[nrow], 0, Constraints->NbColumns-1);
		value_set_si(Constraints->p[row][1+k], 1);
		value_set_si(Constraints->p[nrow][1+k], 1);
		found = 1;
	    }
	    value_clear(tmp);
	    if (found)
		break;
	}
    }
    return found;
}

/**

Given a matrix of constraints ('Constraints'), construct and return a 
polyhedron.

@param Constraints Constraints (may be modified!)
@param NbMaxRays Estimated number of rays in the ray matrix of the
polyhedron.
@return newly allocated Polyhedron

*/ 
Polyhedron *Constraints2Polyhedron(Matrix *Constraints,unsigned NbMaxRays) {
  
  Polyhedron *Pol = NULL;
  Matrix *Ray = NULL;
  SatMatrix *Sat = NULL;
  unsigned Dimension, nbcolumns;
  int i;

  Dimension = Constraints->NbColumns - 1;  /* Homogeneous Dimension */
  if (Dimension < 1) {
    errormsg1("Constraints2Polyhedron","invalidpoly","invalid polyhedron dimension");
    return 0;
  }
  
  /* If there is no constraint in the constraint matrix, return universe */
  /* polyhderon.                                                         */
  if (Constraints->NbRows==0) {  
    Pol = Universe_Polyhedron(Dimension-1);
    return Pol;
  }

  if (POL_ISSET(NbMaxRays, POL_NO_DUAL)) {
    unsigned NbEq;
    unsigned Rank;
    Value tmp;
    if (POL_ISSET(NbMaxRays, POL_INTEGER))
      value_init(tmp);
    do {
      NbEq = 0;
      /* Move equalities up */
      for (i = 0; i < Constraints->NbRows; ++i)
	if (value_zero_p(Constraints->p[i][0])) {
	  if (POL_ISSET(NbMaxRays, POL_INTEGER) &&
	    ConstraintSimplify(Constraints->p[i], 
			       Constraints->p[i], Dimension+1, &tmp)) {
	    value_clear(tmp);
	    return Empty_Polyhedron(Dimension-1);
	  }
	  /* detect 1 == 0, possibly created by ImplicitEqualities */
	  if (First_Non_Zero(Constraints->p[i]+1, Dimension-1) == -1 &&
	      value_notzero_p(Constraints->p[i][Dimension])) {
	    if (POL_ISSET(NbMaxRays, POL_INTEGER))
	      value_clear(tmp);
	    return Empty_Polyhedron(Dimension-1);
	  }
	  if (i != NbEq)
	    ExchangeRows(Constraints, i, NbEq);
	  ++NbEq;
	}
      Rank = Gauss(Constraints, NbEq, Dimension);
      if (POL_ISSET(NbMaxRays, POL_INTEGER))
	for (i = NbEq; i < Constraints->NbRows; ++i)
	  ConstraintSimplify(Constraints->p[i], 
			     Constraints->p[i], Dimension+1, &tmp);
      SortConstraints(Constraints, NbEq);
    } while (ImplicitEqualities(Constraints, NbEq));
    if (POL_ISSET(NbMaxRays, POL_INTEGER))
      value_clear(tmp);
    Pol = Polyhedron_Alloc(Dimension-1, Constraints->NbRows - (NbEq-Rank), 0);
    Vector_Copy(Constraints->p[0], Pol->Constraint[0], 
		Rank * Constraints->NbColumns);
    if (Constraints->NbRows > NbEq)
	Vector_Copy(Constraints->p[NbEq], Pol->Constraint[Rank], 
		    (Constraints->NbRows - NbEq) * Constraints->NbColumns);
    Pol->NbEq = Rank;
    /* Make sure nobody accesses the rays by accident */
    Pol->Ray = 0;
    F_SET(Pol, POL_VALID | POL_INEQUALITIES);
    return Pol;
  }

  if (Dimension > NbMaxRays)
    NbMaxRays = Dimension;
    
  /*
   * Rather than adding a 'positivity constraint', it is better to
   * initialize the lineality space with line in each of the index
   * dimensions, but no line in the lambda dimension. Then initialize
   * the ray space with an origin at 0.  This is what you get anyway,
   * after the positivity constraint has been processed by Chernikova
   * function.
   */

  /* Allocate and initialize the Ray Space */
  Ray = Matrix_Alloc(NbMaxRays, Dimension+1);
  if(!Ray) {
    errormsg1("Constraints2Polyhedron","outofmem","out of memory space");
    return 0;
  }
  Vector_Set(Ray->p_Init,0, NbMaxRays * (Dimension+1));
  for (i=0; i<Dimension; i++) {    
    
    /* Ray->p[i][i+1] = 1 */
    value_set_si(Ray->p[i][i+1],1);
  }

  /* Ray->p[Dimension-1][0] = 1 : mark for ray */
  value_set_si(Ray->p[Dimension-1][0],1);
  Ray->NbRows = Dimension; 

  /* Initialize the Sat Matrix */
  nbcolumns = (Constraints->NbRows - 1)/(sizeof(int)*8) + 1;
  Sat = SMAlloc(NbMaxRays, nbcolumns);
  SMVector_Init(Sat->p_init,Dimension*nbcolumns);
  Sat->NbRows = Dimension;

  CATCH(any_exception_error) {

    /* In case of overflow, free the allocated memory and forward. */    
    if (Sat) SMFree(&Sat);
    if (Ray) Matrix_Free(Ray);
    if (Pol) Polyhedron_Free(Pol);
    RETHROW();
  }
  TRY {

    /* Create ray matrix 'Ray' from constraint matrix 'Constraints' */
    Chernikova(Constraints,Ray,Sat,Dimension-1,NbMaxRays,0,0);

#ifdef POLY_DEBUG
    fprintf(stderr, "[constraints2polyhedron]\nConstraints = ");
    Matrix_Print(stderr,0,Constraints);
    fprintf(stderr, "\nRay = ");
    Matrix_Print(stderr,0,Ray);
    fprintf(stderr, "\nSat = ");
    SMPrint(Sat);
#endif
    
    /* Remove the redundant constraints and create the polyhedron */
    Pol = Remove_Redundants(Constraints,Ray,Sat,0);
  } /* end of TRY */

  UNCATCH(any_exception_error);
  
#ifdef POLY_DEBUG
  fprintf(stderr, "\nPol = ");
  Polyhedron_Print(stderr,"%4d",Pol);
#endif
  
  SMFree(&Sat), Sat = NULL;
  Matrix_Free(Ray), Ray = NULL; 
  return Pol;
} /* Constraints2Polyhedron */

#undef POLY_DEBUG

/* 
 * Given a polyhedron 'Pol', return a matrix of constraints. 
 */
Matrix *Polyhedron2Constraints(Polyhedron *Pol) {
  
  Matrix     *Mat;
  unsigned NbConstraints,Dimension;

  POL_ENSURE_INEQUALITIES(Pol);
  
  NbConstraints = Pol->NbConstraints;
  Dimension     = Pol->Dimension+2;
  Mat = Matrix_Alloc(NbConstraints,Dimension);
  if(!Mat) {
    errormsg1("Polyhedron2Constraints", "outofmem", "out of memory space");
    return 0;
  }
  Vector_Copy(Pol->Constraint[0],Mat->p_Init,NbConstraints * Dimension);
  return Mat;
} /* Polyhedron2Constraints */

/** 

Given a matrix of rays 'Ray', create and return a polyhedron. 

@param Ray Rays (may be modified!)
@param NbMaxConstrs Estimated number of constraints in the polyhedron.
@return newly allocated Polyhedron

*/ 
Polyhedron *Rays2Polyhedron(Matrix *Ray,unsigned NbMaxConstrs) {
  
  Polyhedron *Pol = NULL;
  Matrix *Mat = NULL;
  SatMatrix *Sat = NULL, *SatTranspose = NULL;
  unsigned Dimension, nbcolumns;
  int i;
  
  Dimension = Ray->NbColumns-1;        /* Homogeneous Dimension */
  Sat = NULL;
  SatTranspose = NULL;
  Mat = NULL;
  
  if (Ray->NbRows==0) {  
    
    /* If there is no ray in the matrix 'Ray', return an empty polyhedron */
    Pol = Empty_Polyhedron(Dimension-1);
    return(Pol);
  }

  /* Ignore for now */
  if (POL_ISSET(NbMaxConstrs, POL_NO_DUAL))
    NbMaxConstrs = 0;

  if (Dimension > NbMaxConstrs)
    NbMaxConstrs = Dimension;
  
  /* Allocate space for constraint matrix 'Mat' */
  Mat = Matrix_Alloc(NbMaxConstrs,Dimension+1);
  if(!Mat) {
    errormsg1("Rays2Polyhedron","outofmem","out of memory space");
    return 0;
  }
  
  /* Initialize the constraint matrix 'Mat' */
  Vector_Set(Mat->p_Init,0,NbMaxConstrs * (Dimension+1));
  for (i=0; i<Dimension; i++) {
    
    /* Mat->p[i][i+1]=1 */
    value_set_si(Mat->p[i][i+1],1);
  }
  
  /* Allocate and assign the saturation matrix. Remember we are using a */
  /* transposed saturation matrix referenced by (constraint,ray) pair.  */
  Mat->NbRows = Dimension;
  nbcolumns = (Ray->NbRows -1)/(sizeof(int)*8) + 1;
  SatTranspose = SMAlloc(NbMaxConstrs,nbcolumns);
  SMVector_Init(SatTranspose->p[0],Dimension * nbcolumns);
  SatTranspose->NbRows = Dimension;
  
#ifdef POLY_DEBUG
  fprintf(stderr, "[ray2polyhedron: Before]\nRay = ");
  Matrix_Print(stderr,0,Ray);
  fprintf(stderr, "\nConstraints = ");
  Matrix_Print(stderr,0,Mat);
  fprintf(stderr, "\nSatTranspose = ");
  SMPrint (SatTranspose);
#endif
  
  CATCH(any_exception_error) {
    
    /* In case of overflow, free the allocated memory before forwarding
     * the exception. 
     */
    if (SatTranspose) SMFree(&SatTranspose);
    if (Sat) SMFree(&Sat);
    if (Mat) Matrix_Free(Mat);
    if (Pol) Polyhedron_Free(Pol);
    RETHROW();
  }
  TRY {
    
    /* Create constraint matrix 'Mat' from ray matrix 'Ray' */ 
    Chernikova(Ray,Mat,SatTranspose,Dimension,NbMaxConstrs,0,1);
    
#ifdef POLY_DEBUG
    fprintf(stderr, "[ray2polyhedron: After]\nRay = ");
    Matrix_Print(stderr,0,Ray);
    fprintf(stderr, "\nConstraints = ");
    Matrix_Print(stderr,0,Mat);
    fprintf(stderr, "\nSatTranspose = ");
    SMPrint (SatTranspose);
#endif
    
    /* Transform the saturation matrix 'SatTranspose' in the standard  */
    /* format, that is, ray X constraint format.                       */
    Sat = TransformSat(Mat,Ray,SatTranspose);
    
#ifdef POLY_DEBUG
    fprintf(stderr, "\nSat =");
    SMPrint(Sat);
#endif
    
    SMFree(&SatTranspose), SatTranspose = NULL;
    
    /* Remove redundant rays from the ray matrix 'Ray' */
    Pol = Remove_Redundants(Mat,Ray,Sat,0);
  } /* of TRY */
  
  UNCATCH(any_exception_error);
  
#ifdef POLY_DEBUG
  fprintf(stderr, "\nPol = ");
  Polyhedron_Print(stderr,"%4d",Pol);
#endif
  
  SMFree(&Sat);
  Matrix_Free(Mat);
  return Pol;
} /* Rays2Polyhedron */

/*
 * Compute the dual representation if P only contains one representation.
 * Currently only handles the case where only the equalities are known.
 */
void Polyhedron_Compute_Dual(Polyhedron *P)
{
  if (!F_ISSET(P, POL_VALID))
    return;

  if (F_ISSET(P, POL_FACETS | POL_VERTICES))
    return;

  if (F_ISSET(P, POL_INEQUALITIES)) {
    Matrix M;
    Polyhedron tmp, *Q;
    /* Pretend P is a Matrix for a second */
    M.NbRows = P->NbConstraints;
    M.NbColumns = P->Dimension+2;
    M.p_Init = P->p_Init;
    M.p = P->Constraint;
    Q = Constraints2Polyhedron(&M, 0);

    /* Switch contents of P and Q ... */
    tmp = *Q;
    *Q = *P;
    *P = tmp;
    /* ... but keep the next pointer */
    P->next = Q->next;
    Polyhedron_Free(Q);
    return;
  }

  /* other cases not handled yet */
  assert(0);
}

/*   
 * Build a saturation matrix from constraint matrix 'Mat' and ray matrix 
 * 'Ray'. Only 'NbConstraints' constraint of matrix 'Mat' are considered 
 * in creating the saturation matrix. 'NbMaxRays' is the maximum number 
 * of rows (rays) allowed in the saturation matrix.
 * Vin100's stuff, for the polyparam vertices to work.
 */
static SatMatrix *BuildSat(Matrix *Mat,Matrix *Ray,unsigned NbConstraints,unsigned NbMaxRays) {
  
  SatMatrix *Sat = NULL;
  int i, j, k, jx;
  Value *p1, *p2, *p3;
  unsigned Dimension, NbRay, bx, nbcolumns;
  
  CATCH(any_exception_error) {
    if (Sat) 
      SMFree(&Sat);
    RETHROW();
  }
  TRY {
    NbRay = Ray->NbRows;
    Dimension = Mat->NbColumns-1;   /* Homogeneous Dimension */
    
    /* Build the Sat matrix */
    nbcolumns = (Mat->NbRows - 1)/(sizeof(int)*8) + 1;
    Sat = SMAlloc(NbMaxRays,nbcolumns);
    Sat->NbRows = NbRay;
    SMVector_Init(Sat->p_init, nbcolumns * NbRay);
    jx=0; bx=MSB;
    for (k=0; k<NbConstraints; k++) {
      for (i=0; i<NbRay; i++) {
	
	/* Compute the dot product of ray(i) and constraint(k) and */
	/* store in the status element of ray(i).                  */
	p1 = Ray->p[i]+1;
	p2 = Mat->p[k]+1;
	p3 = Ray->p[i];
	value_set_si(*p3,0);
	for (j=0; j<Dimension; j++) {
	  value_addmul(*p3, *p1, *p2);
	  p1++; p2++;
	}
      }
      for (j=0; j<NbRay; j++) {
	
	/* Set 1 in the saturation matrix if the ray doesn't saturate */
	/* the constraint, otherwise the entry is 0.                  */
	if (value_notzero_p(Ray->p[j][0]))
	  Sat->p[j][jx]|=bx;
      }
      NEXT(jx, bx);
    }
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  return Sat;
} /* BuildSat */

/* 
 * Add 'Nbconstraints' new constraints to polyhedron 'Pol'. Constraints are 
 * pointed by 'Con' and the maximum allowed rays in the new polyhedron is
 * 'NbMaxRays'.   
 */
Polyhedron *AddConstraints(Value *Con,unsigned NbConstraints,Polyhedron *Pol,unsigned NbMaxRays) {

  Polyhedron *NewPol = NULL;
  Matrix   *Mat = NULL, *Ray = NULL;
  SatMatrix *Sat = NULL;
  unsigned NbRay, NbCon, Dimension;

  if (NbConstraints == 0)
    return Polyhedron_Copy(Pol);
  
  POL_ENSURE_FACETS(Pol);
  POL_ENSURE_VERTICES(Pol);

  CATCH(any_exception_error) {
    if (NewPol) Polyhedron_Free(NewPol);
    if (Mat) Matrix_Free(Mat);
    if (Ray) Matrix_Free(Ray);
    if (Sat) SMFree(&Sat);
    RETHROW();
  }
  TRY {
    NbRay	= Pol->NbRays;
    NbCon      	= Pol->NbConstraints + NbConstraints;
    Dimension	= Pol->Dimension + 2;	/* Homogeneous Dimension + Status */

    /* Ignore for now */
    if (POL_ISSET(NbMaxRays, POL_NO_DUAL))
      NbMaxRays = 0;

    if (NbRay > NbMaxRays)
      NbMaxRays = NbRay;
    
    Mat = Matrix_Alloc(NbCon, Dimension);
    if(!Mat) {
      errormsg1("AddConstraints", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /* Copy constraints of polyhedron 'Pol' to matrix 'Mat' */
    Vector_Copy(Pol->Constraint[0], Mat->p[0], Pol->NbConstraints * Dimension);
    
    /* Add the new constraints pointed by 'Con' to matrix 'Mat' */
    Vector_Copy(Con, Mat->p[Pol->NbConstraints], NbConstraints * Dimension);  
    
    /* Allocate space for ray matrix 'Ray' */
    Ray = Matrix_Alloc(NbMaxRays, Dimension);
    if(!Ray) {
      errormsg1("AddConstraints", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    Ray->NbRows = NbRay;

    /* Copy rays of polyhedron 'Pol' to matrix 'Ray' */
    if (NbRay)
	Vector_Copy(Pol->Ray[0], Ray->p[0], NbRay * Dimension);  
    
    /* Create the saturation matrix 'Sat' from constraint matrix 'Mat' and */
    /* ray matrix 'Ray' .                                                  */ 
    Sat = BuildSat(Mat, Ray, Pol->NbConstraints, NbMaxRays);
    
    /* Create the ray matrix 'Ray' from the constraint matrix 'Mat' */
    Pol_status = Chernikova(Mat, Ray, Sat, Pol->NbBid, NbMaxRays, Pol->NbConstraints,0);
    
    /* Remove redundant constraints from matrix 'Mat' */
    NewPol = Remove_Redundants(Mat, Ray, Sat, 0);
    
  } /* end of TRY */
  
  UNCATCH(any_exception_error);  
  SMFree(&Sat);
  Matrix_Free(Ray);
  Matrix_Free(Mat);  
  return NewPol;
} /* AddConstraints */

/* 
 * Return 1 if 'Pol1' includes (covers) 'Pol2', 0 otherwise. 
 * Polyhedron 'A' includes polyhedron 'B' if the rays of 'B' saturate
 * the equalities and verify the inequalities of 'A'. Both 'Pol1' and 
 * 'Pol2' have same dimensions. 
 */
int PolyhedronIncludes(Polyhedron *Pol1,Polyhedron *Pol2) {
	
  int Dimension = Pol1->Dimension + 1;   /* Homogenous Dimension */
  int i, j, k;
  Value *p1, *p2, p3;
  
  POL_ENSURE_FACETS(Pol1);
  POL_ENSURE_VERTICES(Pol1);
  POL_ENSURE_FACETS(Pol2);
  POL_ENSURE_VERTICES(Pol2);

  value_init(p3); 
  for (k=0; k<Pol1->NbConstraints; k++) {
    for (i=0;i<Pol2->NbRays;i++) {
      
      /* Compute the dot product of ray(i) and constraint(k) and store in p3 */
      p1 = Pol2->Ray[i]+1;
      p2 = Pol1->Constraint[k]+1;
      value_set_si(p3,0);
      for(j=0;j<Dimension;j++) {
	value_addmul(p3, *p1,*p2);
	p1++; p2++;
      }
     
      /* If (p3 < 0) or (p3 > 0 and (constraint(k) is equality
                                     or ray(i) is a line)), return 0 */
      if(value_neg_p(p3) ||
          (value_notzero_p(p3)
             && (value_zero_p(Pol1->Constraint[k][0]) || (value_zero_p(Pol2->Ray[i][0]))   ) )) {
	value_clear(p3); 
	return 0;
      }
    }
  } 
  value_clear(p3); 
  return 1;
} /* PolyhedronIncludes */

/*
 * Add Polyhedron 'Pol' to polhedral domain 'PolDomain'. If 'Pol' covers
 * some polyhedron in the domain 'PolDomain', it is removed from the list.
 * On the other hand if some polyhedron in the domain covers polyhedron 
 * 'Pol' then 'Pol' is not included in the domain.   
 */
Polyhedron *AddPolyToDomain(Polyhedron *Pol,Polyhedron *PolDomain) {
  
  Polyhedron *p, *pnext, *p_domain_end = (Polyhedron *) 0;
  int Redundant;
  
  if (!Pol) 
    return PolDomain;
  if (!PolDomain)	
    return Pol;

  POL_ENSURE_FACETS(Pol);
  POL_ENSURE_VERTICES(Pol);
  
  /* Check for emptiness of polyhedron 'Pol' */
  if (emptyQ(Pol)) {
    Polyhedron_Free(Pol);
    return PolDomain;
  }
  
  POL_ENSURE_FACETS(PolDomain);
  POL_ENSURE_VERTICES(PolDomain);

  /* Check for emptiness of polyhedral domain 'PolDomain' */
  if (emptyQ(PolDomain)) {
    Polyhedron_Free(PolDomain);
    return Pol;
  }
  
  /* Test 'Pol' against the domain 'PolDomain' */
  Redundant = 0;
  for (p=PolDomain,PolDomain=(Polyhedron *)0; p; p=pnext) {
    
    /* If 'Pol' covers 'p' */    
    if (PolyhedronIncludes(Pol, p))
    {
       /* free p */
		 pnext = p->next;
       Polyhedron_Free( p );
       continue;
    }

    /* Add polyhedron p to the new domain list */
    if (!PolDomain) PolDomain = p; else p_domain_end->next = p;
    p_domain_end = p;
    
    /* If p covers Pol */
    if (PolyhedronIncludes(p,Pol)) {
      Redundant = 1;
      break;
    }
    pnext = p->next;
  }
  if (!Redundant) {  
    
    /* The whole list has been checked. Add new polyhedron 'Pol' to the */
    /* new domain list.                                                 */ 
    if (!PolDomain) PolDomain = Pol; else p_domain_end->next = Pol;
  }
  else {
    
    /* The rest of the list is just inherited from p */
    Polyhedron_Free(Pol);
  }
  return PolDomain;
} /* AddPolyToDomain */

/* 
 * Given a polyhedra 'Pol' and a single constraint 'Con' and an integer 'Pass' 
 * whose value ranges from 0 to 3, add the inverse of constraint 'Con' to the 
 * constraint set of 'Pol' and return the new polyhedron. 'NbMaxRays' is the 
 * maximum allowed rays in the new generated polyhedron. 
 * If Pass == 0, add ( -constraint -1) >= 0
 * If Pass == 1, add ( +constraint -1) >= 0
 * If Pass == 2, add ( -constraint   ) >= 0
 * If Pass == 3, add ( +constraint   ) >= 0
 */
Polyhedron *SubConstraint(Value *Con,Polyhedron *Pol,unsigned NbMaxRays,int Pass) {
  
  Polyhedron *NewPol = NULL;
  Matrix   *Mat = NULL, *Ray = NULL;
  SatMatrix *Sat = NULL;
  unsigned NbRay, NbCon, NbEle1, Dimension;
  int i;

  POL_ENSURE_FACETS(Pol);
  POL_ENSURE_VERTICES(Pol);
  
  CATCH(any_exception_error) {
    if (NewPol) Polyhedron_Free(NewPol);
    if (Mat) Matrix_Free(Mat);
    if (Ray) Matrix_Free(Ray);
    if (Sat) SMFree(&Sat);
    RETHROW();
  }
  TRY {
    
    /* If 'Con' is the positivity constraint, return Null */
    Dimension  = Pol->Dimension+1;      /* Homogeneous Dimension */
    for (i=1; i<Dimension; i++)
      if (value_notzero_p(Con[i])) break;
    if (i==Dimension) {
      UNCATCH(any_exception_error);
      return (Polyhedron *) 0;
    }
    
    NbRay     = Pol->NbRays;
    NbCon     = Pol->NbConstraints;
    Dimension = Pol->Dimension+2;	/* Homogeneous Dimension + Status */
    NbEle1    = NbCon * Dimension;
    
    /* Ignore for now */
    if (POL_ISSET(NbMaxRays, POL_NO_DUAL))
      NbMaxRays = 0;

    if (NbRay > NbMaxRays)
      NbMaxRays = NbRay;
    
    Mat = Matrix_Alloc(NbCon + 1, Dimension);
    if(!Mat) {
      errormsg1("SubConstraint", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /* Set the constraints of Pol */
    Vector_Copy(Pol->Constraint[0], Mat->p[0], NbEle1);
    
    /* Add the new constraint */
    value_set_si(Mat->p[NbCon][0],1);
    if (!(Pass&1))
      for(i=1; i<Dimension; i++) 
	value_oppose(Mat->p[NbCon][i],Con[i]);
    else
      for(i=1; i<Dimension; i++)
	value_assign(Mat->p[NbCon][i],Con[i]);
    if (!(Pass&2))
      value_decrement(Mat->p[NbCon][Dimension-1],Mat->p[NbCon][Dimension-1]);
   
    /* Allocate the ray matrix. */
    Ray = Matrix_Alloc(NbMaxRays, Dimension);
    if(!Ray) {
      errormsg1("SubConstraint", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /* Initialize the ray matrix with the rays of polyhedron 'Pol' */
    Ray->NbRows = NbRay;
    if (NbRay)
	Vector_Copy(Pol->Ray[0], Ray->p[0], NbRay * Dimension);   
    
    /* Create the saturation matrix from the constraint matrix 'mat' and */
    /* ray matrix 'Ray'.                                                 */
    Sat = BuildSat(Mat, Ray, NbCon, NbMaxRays);
    
    /* Create the ray matrix 'Ray' from consraint matrix 'Mat'           */
    Pol_status = Chernikova(Mat, Ray, Sat, Pol->NbBid, NbMaxRays, NbCon,0);
    
    /* Remove redundant constraints from matrix 'Mat' */ 
    NewPol = Remove_Redundants(Mat, Ray, Sat, 0);
    
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  
  SMFree(&Sat);
  Matrix_Free(Ray);
  Matrix_Free(Mat);
  return NewPol;
} /* SubConstraint */

/*
 * Return the intersection of two polyhedral domains 'Pol1' and 'Pol2'. 
 * The maximum allowed rays in the new polyhedron generated is 'NbMaxRays'. 
 */
Polyhedron *DomainIntersection(Polyhedron *Pol1,Polyhedron *Pol2,unsigned NbMaxRays) {
  
  Polyhedron *p1, *p2, *p3, *d;
  
  if (!Pol1 || !Pol2) return (Polyhedron*) 0;
  if (Pol1->Dimension != Pol2->Dimension) {
    errormsg1( "DomainIntersection", "diffdim",
	       "operation on different dimensions");
    return (Polyhedron*) 0;
  }
  
  /* For every polyhedron pair (p1,p2) where p1 belongs to domain Pol1 and */
  /* p2 belongs to domain Pol2, compute the intersection and add it to the */
  /* new domain 'd'.                                                       */
  d = (Polyhedron *)0;
  for (p1=Pol1; p1; p1=p1->next) {
    for (p2=Pol2; p2; p2=p2->next) {
      p3 = AddConstraints(p2->Constraint[0],
			  p2->NbConstraints, p1, NbMaxRays);	  
      d = AddPolyToDomain(p3,d);
    }
  }
  if (!d)
    return Empty_Polyhedron(Pol1->Dimension);
  else
    return d;
  
} /* DomainIntersection */

/*
 * Given a polyhedron 'Pol', return a matrix of rays. 
 */
Matrix *Polyhedron2Rays(Polyhedron *Pol) {
  
  Matrix     *Ray;
  unsigned NbRays, Dimension;

  POL_ENSURE_POINTS(Pol);
  
  NbRays    = Pol->NbRays;
  Dimension = Pol->Dimension+2;		/* Homogeneous Dimension + Status */
  Ray = Matrix_Alloc(NbRays, Dimension);
  if(!Ray) {
    errormsg1("Polyhedron2Rays", "outofmem", "out of memory space");
    return 0;
  }
  Vector_Copy(Pol->Ray[0], Ray->p_Init, NbRays*Dimension);
  return Ray;
} /* Polyhedron2Rays */

/*
 * Add 'NbAddedRays' rays to polyhedron 'Pol'. Rays are pointed by 'AddedRays'
 * and the maximum allowed constraints in the new polyhedron is 'NbMaxConstrs'.
 */ 
Polyhedron *AddRays(Value *AddedRays,unsigned NbAddedRays,Polyhedron *Pol,unsigned NbMaxConstrs) {

  Polyhedron *NewPol = NULL;
  Matrix   *Mat = NULL, *Ray = NULL;
  SatMatrix *Sat = NULL, *SatTranspose = NULL;
  unsigned NbCon, NbRay,NbEle1, Dimension;
  
  POL_ENSURE_FACETS(Pol);
  POL_ENSURE_VERTICES(Pol);

  CATCH(any_exception_error) {
    if (NewPol) Polyhedron_Free(NewPol);
    if (Mat) Matrix_Free(Mat);
    if (Ray) Matrix_Free(Ray);
    if (Sat) SMFree(&Sat);
    if (SatTranspose) SMFree(&SatTranspose);
    RETHROW();
  }
  TRY {
    
    NbCon      = Pol->NbConstraints;
    NbRay      = Pol->NbRays;
    Dimension  = Pol->Dimension + 2;	/* Homogeneous Dimension + Status */
    NbEle1     = NbRay * Dimension;
    
    Ray = Matrix_Alloc(NbAddedRays + NbRay, Dimension);
    if(!Ray) {
      errormsg1("AddRays", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /* Copy rays of polyhedron 'Pol' to matrix 'Ray' */
    if (NbRay)
      Vector_Copy(Pol->Ray[0], Ray->p_Init, NbEle1);
    
    /* Add the new rays pointed by 'AddedRays' to matrix 'Ray' */
    Vector_Copy(AddedRays, Ray->p_Init+NbEle1, NbAddedRays * Dimension);
    
    /* Ignore for now */
    if (POL_ISSET(NbMaxConstrs, POL_NO_DUAL))
      NbMaxConstrs = 0;

    /* We need at least NbCon rows */
    if (NbMaxConstrs < NbCon)
	NbMaxConstrs = NbCon;

    /* Allocate space for constraint matrix 'Mat' */
    Mat = Matrix_Alloc(NbMaxConstrs, Dimension);
    if(!Mat) {
      errormsg1("AddRays", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    Mat->NbRows = NbCon;
    
    /* Copy constraints of polyhedron 'Pol' to matrix 'Mat' */
    Vector_Copy(Pol->Constraint[0], Mat->p_Init, NbCon*Dimension);

    /* Create the saturation matrix 'SatTranspose' from constraint matrix */
    /* 'Mat' and ray matrix 'Ray'. Remember the saturation matrix is      */
    /* referenced by (constraint,ray) pair                                */ 
    SatTranspose = BuildSat(Ray, Mat, NbRay, NbMaxConstrs);
    
    /* Create the constraint matrix 'Mat' from the ray matrix 'Ray' */
    Pol_status = Chernikova(Ray, Mat, SatTranspose, Pol->NbEq, NbMaxConstrs, NbRay,1);
    
    /* Transform the saturation matrix 'SatTranspose' in the standard format */
    /* , that is, (ray X constraint) format.                                 */
    Sat = TransformSat(Mat, Ray, SatTranspose);
    SMFree(&SatTranspose), SatTranspose = NULL;
    
    /* Remove redundant rays from the ray matrix 'Ray' */
    NewPol = Remove_Redundants(Mat, Ray, Sat, 0);
    
    SMFree(&Sat), Sat = NULL;
    Matrix_Free(Mat), Mat = NULL;
    Matrix_Free(Ray), Ray = NULL;  
  } /* end of TRY */
  
  UNCATCH(any_exception_error);  
  return NewPol;
} /* AddRays */

/* 
 * Add rays pointed by 'Ray' to each and every polyhedron in the polyhedral 
 * domain 'Pol'. 'NbMaxConstrs' is maximum allowed constraints in the 
 * constraint set of a polyhedron.                         
 */ 
Polyhedron *DomainAddRays(Polyhedron *Pol,Matrix *Ray,unsigned NbMaxConstrs) {
  
  Polyhedron *PolA, *PolEndA, *p1, *p2, *p3;
  int Redundant;
  
  if (!Pol) return (Polyhedron*) 0;
  if (!Ray || Ray->NbRows == 0)
    return Domain_Copy(Pol);
  if (Pol->Dimension != Ray->NbColumns-2) {
    errormsg1(
	      "DomainAddRays", "diffdim", "operation on different dimensions");
    return (Polyhedron*) 0;
  }
  
  /* Copy Pol to PolA */
  PolA = PolEndA = (Polyhedron *)0;
  for (p1=Pol; p1; p1=p1->next) {
    p3 = AddRays(Ray->p[0], Ray->NbRows, p1, NbMaxConstrs);
    
    /* Does any component of PolA cover 'p3' ? */
    Redundant = 0;
    for (p2=PolA; p2; p2=p2->next) {
      if (PolyhedronIncludes(p2, p3)) { /* If p2 covers p3 */
	Redundant = 1;
	break;
      }
    }
    
    /* If the new polyheron is not redundant, add it ('p3') to the list */
    if (Redundant)
      Polyhedron_Free(p3);
    else { 
      if (!PolEndA)
	PolEndA = PolA = p3;
      else {
	PolEndA->next = p3;
	PolEndA = PolEndA->next;
      }
    }
  }
  return PolA;
} /* DomainAddRays */

/*
 * Create a copy of the polyhedron 'Pol' 
 */
Polyhedron *Polyhedron_Copy(Polyhedron *Pol) {
  
  Polyhedron *Pol1;
  
  if (!Pol) return (Polyhedron *)0;
  
  /* Allocate space for the new polyhedron */
  Pol1 = Polyhedron_Alloc(Pol->Dimension, Pol->NbConstraints, Pol->NbRays);
  if (!Pol1) {
    errormsg1("Polyhedron_Copy", "outofmem", "out of memory space");
    return 0;
  }
  if( Pol->NbConstraints )
    Vector_Copy(Pol->Constraint[0], Pol1->Constraint[0],
	      Pol->NbConstraints*(Pol->Dimension+2));
  if( Pol->NbRays )
    Vector_Copy(Pol->Ray[0], Pol1->Ray[0],
	      Pol->NbRays*(Pol->Dimension+2));
  Pol1->NbBid = Pol->NbBid;
  Pol1->NbEq = Pol->NbEq;
  Pol1->flags = Pol->flags;
  return Pol1;
} /* Polyhedron_Copy */

/* 
 * Create a copy of a polyhedral domain. 
 */
Polyhedron *Domain_Copy(Polyhedron *Pol) {
  
  Polyhedron *Pol1;
  
  if (!Pol) return (Polyhedron *) 0;
  Pol1 = Polyhedron_Copy(Pol);
  if (Pol->next) Pol1->next = Domain_Copy(Pol->next);
  return Pol1;
} /* Domain_Copy */

/*
 * Given constraint number 'k' of a polyhedron, and an array 'Filter' to store
 * the non-redundant constraints of the polyhedron in bit-wise notation, and
 * a Matrix 'Sat', add the constraint 'k' in 'Filter' array. tmpR[i] stores the
 * number of constraints, other than those in 'Filter', which ray(i) saturates 
 * or verifies. In case, ray(i) does not saturate or verify a constraint in
 * array 'Filter', it is assigned to -1. Similarly, tmpC[j] stores the number
 * of rays which constraint(j), if it doesn't belong to Filter, saturates or 
 * verifies. If constraint(j) belongs to 'Filter', then tmpC[j] is assigned to
 * -1. 'NbConstraints' is the number of constraints in the constraint matrix of
 * the polyhedron. 
 * NOTE: (1) 'Sat' is not the saturation matrix of the polyhedron. In fact, 
 *           entry in 'Sat' is set to 1 if ray(i) of polyhedron1 verifies or 
 *           saturates the constraint(j) of polyhedron2 and otherwise it is set
 *           to 0. So here the difference with saturation matrix is in terms 
 *           definition and entries(1<->0) of the matrix 'Sat'.    
 *       
 * ALGORITHM:->
 * (1) Include constraint(k) in array 'Filter'. 
 * (2) Set tmpC[k] to -1.
 * (3) For all ray(i) {
 *        If ray(i) saturates or verifies constraint(k) then --(tmpR[i])
 *        Else {
 *           Discard ray(i) by assigning tmpR[i] = -1
 *           Decrement tmpC[j] for all constraint(j) not in array 'Filter'.
 *        }
 *     }
 */
static void addToFilter(int k, unsigned *Filter, SatMatrix *Sat,Value *tmpR,Value *tmpC,int NbRays,int NbConstraints) {
  
  int kj, i,j, jx;
  unsigned kb, bx;
  
  /* Remove constraint k */
  kj =   k/WSIZE; kb = MSB; kb >>= k%WSIZE;
  Filter[kj]|=kb;
  value_set_si(tmpC[k],-1);
  
  /* Remove rays excluded by constraint k */
  for(i=0; i<NbRays; i++)
    if (value_posz_p(tmpR[i])) {
      if (Sat->p[i][kj]&kb)
	value_decrement(tmpR[i],tmpR[i]);  /* adjust included ray */
      else {
	
	/* Constraint k excludes ray i -- delete ray i */
	value_set_si(tmpR[i],-1);
	
	/* Adjust non-deleted constraints */
	jx=0; bx=MSB;
	for(j=0; j<NbConstraints; j++) {
	  if (value_posz_p(tmpC[j]) && (Sat->p[i][jx]&bx) )
	    value_decrement(tmpC[j],tmpC[j]);
	  NEXT(jx,bx);
	}
      }
    } 
} /* addToFilter */

/*
 * Given polyhedra 'P1' and 'P2' such that their intersection is an empty
 * polyhedron, find the minimal set of constraints of 'P1' which contradict
 * all of the constraints of 'P2'. This is believed to be an NP-hard problem
 * and so a heuristic is employed to solve it in worst case. The heuristic is 
 * to select in every turn that constraint of 'P1' which excludes most rays of
 * 'P2'. A bit in the binary format of an element of array 'Filter' is set to
 * 1 if the corresponding constraint is to be included in the minimal set of 
 * constraints otherwise it is set to 0.
 */
static void FindSimple(Polyhedron *P1,Polyhedron *P2,unsigned *Filter,unsigned NbMaxRays) {
  
  Matrix *Mat = NULL;
  SatMatrix *Sat = NULL;
  int i, j, k, jx, found;
  Value *p1, *p2, p3;
  unsigned Dimension, NbRays, NbConstraints, bx, nc;
  Value NbConstraintsLeft, tmp;
  Value *tmpC = NULL, *tmpR = NULL;
  Polyhedron *Pol = NULL, *Pol2 = NULL;
  
  /* Initialize all the 'Value' variables */
  value_init(p3); value_init(NbConstraintsLeft);
  value_init(tmp);
 
  CATCH(any_exception_error) {
    if (tmpC) free(tmpC);
    if (tmpR) free(tmpR);
    if (Mat) Matrix_Free(Mat);
    if (Sat) SMFree(&Sat);
    if (Pol2 && Pol2!=P2) Polyhedron_Free(Pol2);
    if (Pol && Pol!=Pol2 && Pol!=P2) Polyhedron_Free(Pol);
    
    /* Clear all the 'Value' variables */
    value_clear(p3); value_clear(NbConstraintsLeft);
    value_clear(tmp);
    RETHROW();
  }
  TRY {
    
    Dimension = P1->Dimension+2;       /* status + homogeneous Dimension */
    Mat = Matrix_Alloc(P1->NbConstraints, Dimension);
    if(!Mat) {
      errormsg1("FindSimple", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      
      /* Clear all the 'Value' variables */
      value_clear(p3); value_clear(NbConstraintsLeft); value_clear(tmp);  
      return;
    }
    
    /* Post constraints in P1 already included by Filter */
    jx = 0; bx = MSB; Mat->NbRows=0;
    value_set_si(NbConstraintsLeft,0);
    for (k=0; k<P1->NbConstraints; k++) {
      if (Filter[jx]&bx) {
	Vector_Copy(P1->Constraint[k], Mat->p[Mat->NbRows], Dimension);
	Mat->NbRows++;
      }
      else
	value_increment(NbConstraintsLeft,NbConstraintsLeft);
      NEXT(jx,bx);
    }
    Pol2 = P2;
    
    for (;;) {
      if (Mat->NbRows==0)
	Pol = Polyhedron_Copy(Pol2);
      else {
	Pol = AddConstraints(Mat->p_Init, Mat->NbRows, Pol2, NbMaxRays);
	if (Pol2 != P2) Polyhedron_Free(Pol2), Pol2 = NULL;
      }
      if (emptyQ(Pol)) {
	Matrix_Free(Mat), Mat = NULL;
	Polyhedron_Free(Pol), Pol = NULL;
	UNCATCH(any_exception_error);
	
	/* Clear all the 'Value' variables */
	value_clear(p3); value_clear(NbConstraintsLeft); value_clear(tmp);
	return;
      }
      Mat->NbRows = 0;        /* Reset Mat */
      Pol2 = Pol;
      
      /* Its not enough-- find some more constraints */
      Dimension         = Pol->Dimension+1;       /* homogeneous Dimension */
      NbRays            = Pol->NbRays;
      NbConstraints     = P1->NbConstraints;
      tmpR = (Value *)malloc(NbRays*sizeof(Value));
      if(!tmpR) {
	errormsg1("FindSimple", "outofmem", "out of memory space");
	UNCATCH(any_exception_error);
	
	/* Clear all the 'Value' variables */
	value_clear(p3); value_clear(NbConstraintsLeft); value_clear(tmp);  
	return;
      }
      for(i=0;i<NbRays;i++)
	value_init(tmpR[i]);
      tmpC = (Value *)malloc(NbConstraints*sizeof(Value));
      if(!tmpC) {
	errormsg1("FindSimple", "outofmem", "out of memory space");
	UNCATCH(any_exception_error);
	
	/* Clear all the 'Value' variables */
	value_clear(p3); value_clear(NbConstraintsLeft);
	for(i=0;i<NbRays;i++)
	  value_clear(tmpR[i]);
	free(tmpR);
	return;
      }
      for(i=0;i<NbConstraints;i++)
	value_init(tmpC[i]);
      Vector_Set(tmpR,0,NbRays);
      Vector_Set(tmpC,0,NbConstraints);
      
      /* Build the Sat matrix */
      nc      = (NbConstraints - 1) / (sizeof(int)*8) + 1;
      Sat     = SMAlloc(NbRays, nc);
      Sat->NbRows = NbRays;
      SMVector_Init(Sat->p_init, nc*NbRays);
      
      jx=0; bx=MSB;
      for (k=0; k<NbConstraints; k++) {
	if (Filter[jx]&bx)
	  value_set_si(tmpC[k],-1);
	else
	  for (i=0; i<NbRays; i++) {
	    p1 = Pol->Ray[i]+1;
	    p2 = P1->Constraint[k]+1;
	    value_set_si(p3,0);
	    for (j=0; j<Dimension; j++) {
	      value_addmul(p3, *p1, *p2);
	      p1++; p2++;
	    }
	    if(value_zero_p(p3) ||
	       (value_pos_p(p3) && value_notzero_p(P1->Constraint[k][0]))) {
	      Sat->p[i][jx]|=bx;  /* constraint includes ray, set flag */
	      value_increment(tmpR[i],tmpR[i]);
	      value_increment(tmpC[k],tmpC[k]);
	    }
	  }
	NEXT(jx, bx);
      }
      
      do { /* find all of the essential constraints */
	found = 0;
	for(i=0; i<NbRays; i++)
	  if(value_posz_p(tmpR[i])) {
	    value_add_int(tmp,tmpR[i],1);
	    if(value_eq(tmp,NbConstraintsLeft)) {
	      
	      /* Ray i is excluded by only one constraint... find it */
	      jx = 0; bx = MSB;
	      for(k=0; k<NbConstraints; k++) {
		if(value_posz_p(tmpC[k]) && ((Sat->p[i][jx]&bx)==0)) {
		  addToFilter(k, Filter, Sat, tmpR, tmpC,
			      NbRays, NbConstraints);
		  Vector_Copy(P1->Constraint[k],
			      Mat->p[Mat->NbRows],Dimension+1);
		  Mat->NbRows++;
		  value_decrement(NbConstraintsLeft,NbConstraintsLeft);
		  found=1;
		  break;
		}
		NEXT(jx,bx);
	      }
	      break;
	    }
	  }
      }
      while (found);
     
      if (!Mat->NbRows) { /* Well then, just use a stupid heuristic */
	/* find the constraint which excludes the most */
	Value cmax;
	value_init(cmax);
	
#ifndef LINEAR_VALUE_IS_CHARS
        value_set_si(cmax,(NbRays * NbConstraints+1));
#else
	value_set_si(cmax,1);
#endif
	
	j = -1;
	for(k=0; k<NbConstraints; k++)
	  if (value_posz_p(tmpC[k])) {
	    if (value_gt(cmax,tmpC[k])) {
	      value_assign(cmax,tmpC[k]);
	      j = k;
	    }
	  }
	value_clear(cmax);
	if (j<0) {
	  errormsg1("DomSimplify","logerror","logic error");
	}
	else {
	  addToFilter(j, Filter, Sat, tmpR, tmpC, NbRays, NbConstraints);
	  Vector_Copy(P1->Constraint[j],Mat->p[Mat->NbRows],Dimension+1);
	  Mat->NbRows++;
	  value_decrement(NbConstraintsLeft,NbConstraintsLeft);
	}
      }
      SMFree(&Sat), Sat = NULL;
      free(tmpC), tmpC = NULL;
      free(tmpR), tmpR = NULL;
    }   
  } /* end of TRY */
  
  /* Clear all the 'Value' variables */
  value_clear(p3); value_clear(NbConstraintsLeft);
  value_clear(tmp);
  for(i=0;i<NbRays;i++)
    value_clear(tmpR[i]);
  for(i=0;i<NbRays;i++)
    value_clear(tmpC[i]);
  
  UNCATCH(any_exception_error);
} /* FindSimple */

/* 
 * Return 0 if the intersection of Pol1 and Pol2 is empty, otherwise return 1.
 * If the intersection is non-empty, store the non-redundant constraints in 
 * 'Filter' array. If the intersection is empty then store the smallest set of
 * constraints of 'Pol1' which on intersection with 'Pol2' gives empty set, in
 * 'Filter' array. 'NbMaxRays' is the maximum allowed rays in the intersection
 *  of 'Pol1' and 'Pol2'.   
 */
static int SimplifyConstraints(Polyhedron *Pol1,Polyhedron *Pol2,unsigned *Filter,unsigned NbMaxRays) {
  
  Polyhedron *Pol = NULL;
  Matrix   *Mat = NULL, *Ray = NULL;
  SatMatrix *Sat = NULL;
  unsigned NbRay, NbCon, NbCon1, NbCon2, NbEle1, Dimension, notempty;

  CATCH(any_exception_error) {
    if (Pol) Polyhedron_Free(Pol);
    if (Mat) Matrix_Free(Mat);
    if (Ray) Matrix_Free(Ray);
    if (Sat) SMFree(&Sat);
    RETHROW();
  }
  TRY {

    NbRay         = Pol1->NbRays;
    NbCon1        = Pol1->NbConstraints;
    NbCon2        = Pol2->NbConstraints;
    NbCon         = NbCon1 + NbCon2;
    Dimension     = Pol1->Dimension+2;    /* Homogeneous Dimension + Status */
    NbEle1        = NbCon1*Dimension;
    
    /* Ignore for now */
    if (POL_ISSET(NbMaxRays, POL_NO_DUAL))
      NbMaxRays = 0;

    if (NbRay > NbMaxRays)
      NbMaxRays = NbRay;

    /* Allocate space for constraint matrix 'Mat' */
    Mat = Matrix_Alloc(NbCon, Dimension);
    if(!Mat) {
      errormsg1("SimplifyConstraints", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }

    /* Copy constraints of 'Pol1' to matrix 'Mat' */
    Vector_Copy(Pol1->Constraint[0], Mat->p_Init, NbEle1);
    
    /* Add constraints of 'Pol2' to matrix 'Mat'*/
    Vector_Copy(Pol2->Constraint[0], Mat->p_Init+NbEle1, NbCon2*Dimension);

    /* Allocate space for ray matrix 'Ray' */
    Ray = Matrix_Alloc(NbMaxRays, Dimension);
    if(!Ray) {
      errormsg1("SimplifyConstraints", "outofmem", "out of memory space");
      UNCATCH(any_exception_error);
      return 0;
    }
    Ray->NbRows = NbRay;
	
    /* Copy rays of polyhedron 'Pol1' to matrix 'Ray' */
    Vector_Copy(Pol1->Ray[0], Ray->p_Init, NbRay*Dimension);

    /* Create saturation matrix from constraint matrix 'Mat' and ray matrix */
    /* 'Ray'.                                                               */
    Sat = BuildSat(Mat, Ray, NbCon1, NbMaxRays);

    /* Create the ray matrix 'Ray' from the constraint matrix 'Mat' */
    Pol_status = Chernikova(Mat, Ray, Sat, Pol1->NbBid, NbMaxRays, NbCon1,0);

    /* Remove redundant constraints from the constraint matrix 'Mat' */
    Pol = Remove_Redundants(Mat, Ray, Sat, Filter);
    notempty = 1;
    if (Filter && emptyQ(Pol)) {
      notempty = 0;
      FindSimple(Pol1, Pol2, Filter, NbMaxRays);
    }
    /* Polyhedron_Print(stderr,"%4d",Pol1); */

    Polyhedron_Free(Pol), Pol = NULL;
    SMFree(&Sat), Sat = NULL;
    Matrix_Free(Ray), Ray = NULL;
    Matrix_Free(Mat), Mat = NULL;
    
  } /* end of TRY */

  UNCATCH(any_exception_error);  
  return notempty;
} /* SimplifyConstraints */

/* 
 * Eliminate equations of Pol1 using equations of Pol2. Mark as needed, 
 * equations of Pol1 that are not eliminated. Or info into Filter vector. 
 */
static void SimplifyEqualities(Polyhedron *Pol1, Polyhedron *Pol2, unsigned *Filter) {

  int i,j;
  unsigned ix, bx, NbEqn, NbEqn1, NbEqn2, NbEle2, Dimension;
  Matrix   *Mat;

  NbEqn1        = Pol1->NbEq;
  NbEqn2	= Pol2->NbEq;
  NbEqn         = NbEqn1 + NbEqn2;
  Dimension     = Pol1->Dimension+2;    /* Homogeneous Dimension + Status */
  NbEle2        = NbEqn2*Dimension;

  Mat = Matrix_Alloc(NbEqn, Dimension);
  if (!Mat) {
    errormsg1("SimplifyEqualities", "outofmem", "out of memory space");
    Pol_status = 1;
    return;
  }

  /* Set the equalities of Pol2 */
  Vector_Copy(Pol2->Constraint[0], Mat->p_Init, NbEle2);

  /* Add the equalities of Pol1 */
  Vector_Copy(Pol1->Constraint[0], Mat->p_Init+NbEle2, NbEqn1*Dimension);

  Gauss(Mat, NbEqn2, Dimension-1);

  ix = 0;
  bx = MSB;
  for (i=NbEqn2; i<NbEqn; i++) {
    for (j=1; j<Dimension; j++) {
      if (value_notzero_p(Mat->p[i][j])) { 
	/* If any coefficient of the equation is non-zero */	
	/* Set the filter bit for the equation */
	
	Filter[ix] |= bx;
	break;
      }
    }
    NEXT(ix,bx);
  }
  Matrix_Free(Mat);
  return;
} /* SimplifyEqualities */

 
/* 
 * Given two polyhedral domains 'Pol1' and 'Pol2', find the largest domain
 * set (or the smallest list of non-redundant constraints), that when 
 * intersected with polyhedral domain 'Pol2' equals (Pol1)intersect(Pol2).
 * The output is a polyhedral domain with the "redundant" constraints removed.
 * 'NbMaxRays' is the maximium allowed rays in the new polyhedra. 
 */
Polyhedron *DomainSimplify(Polyhedron *Pol1, Polyhedron *Pol2, unsigned NbMaxRays) {
  
  Polyhedron *p1, *p2, *p3, *d;
  unsigned k, jx, bx, nbentries, NbConstraints, Dimension, NbCon, empty;
  unsigned *Filter;
  Matrix *Constraints;
  

  if (!Pol1 || !Pol2) return Pol1;
  if (Pol1->Dimension != Pol2->Dimension) {
    errormsg1("DomSimplify","diffdim","operation on different dimensions");
    Pol_status = 1;
    return 0;
  }
  POL_ENSURE_VERTICES(Pol1);
  POL_ENSURE_VERTICES(Pol2);
  if (emptyQ(Pol1)||emptyQ(Pol2)) 
    return Empty_Polyhedron(Pol1->Dimension);

  /* Find the maximum number of constraints over all polyhedron in the  */
  /* polyhedral domain 'Pol2' and store in 'NbCon'.                     */
  NbCon = 0;
  for (p2=Pol2; p2; p2=p2->next)
    if (p2->NbConstraints > NbCon) 
      NbCon = p2->NbConstraints;
  
  Dimension = Pol1->Dimension+2;     /* Homogenous Dimension + Status  */
  d = (Polyhedron *)0;
  for (p1=Pol1; p1; p1=p1->next) { 
    
    /* Filter is an array of integers, each bit in an element of Filter */
    /* array corresponds to a constraint. The bit is marked 1 if the    */
    /* corresponding constraint is non-redundant and is 0 if it is      */
    /* redundant.                                                       */
    
    NbConstraints = p1->NbConstraints;
    nbentries = (NbConstraints + NbCon - 1) / (sizeof(int)*8) + 1;

    /* Allocate space for array 'Filter' */
    Filter  = (unsigned *)malloc(nbentries * sizeof(int));
    if (!Filter) {
      errormsg1("DomSimplify", "outofmem", "out of memory space\n");
      Pol_status = 1;
      return 0;
    } 
    
    /* Initialize 'Filter' with zeros */
    SMVector_Init(Filter, nbentries);
    
    /* Filter the constraints of p1 in context of polyhedra p2(s) */
    empty = 1;
    for (p2=Pol2; p2; p2=p2->next) {
      
      /* Store the non-redundant constraints in array 'Filter'. With    */
      /* successive loops, the array 'Filter' holds the union of all    */
      /* non-redundant constraints. 'empty' is set to zero if the       */
      /* intersection of two polyhedra is non-empty and Filter is !Null */
      
      SimplifyEqualities(p1, p2, Filter);
      if (SimplifyConstraints(p1, p2, Filter, NbMaxRays)) 
	empty=0;      
            
      /* takes the union of all non redundant constraints */
    }

    if (!empty) {
      
      /* Copy all non-redundant constraints to matrix 'Constraints' */
      Constraints = Matrix_Alloc(NbConstraints, Dimension);
      if (!Constraints) {
	errormsg1("DomSimplify", "outofmem", "out of memory space\n");
	Pol_status = 1;
	return 0;
      }
      Constraints->NbRows = 0;
      for (k=0, jx=0, bx=MSB; k<NbConstraints; k++) {

	/* If a bit entry in Filter[jx] is marked 1, copy the correspond- */
	/* ing constraint in matrix 'Constraints'.                        */
	if (Filter[jx]&bx) { 
	  Vector_Copy(p1->Constraint[k],
		      Constraints->p[Constraints->NbRows],
		      Dimension);
	  Constraints->NbRows++;
	}
	NEXT(jx,bx);
      }
      
      /* Create the polyhedron 'p3' corresponding to the constraints in   */
      /* matrix 'Constraints'.                                            */
      p3 = Constraints2Polyhedron(Constraints,NbMaxRays);
      Matrix_Free(Constraints);
      
      /* Add polyhedron 'p3' in the domain 'd'. */
      d = AddPolyToDomain (p3, d);
      p3 = NULL;
    }
    free(Filter);
  }
  if (!d) 
    return Empty_Polyhedron(Pol1->Dimension); 
  else return d;

} /* DomainSimplify */

/*
 * Domain Simplify as defined in Strasborg Polylib version. 
 */
Polyhedron *Stras_DomainSimplify(Polyhedron *Pol1,Polyhedron *Pol2,unsigned NbMaxRays) {

  Polyhedron *p1, *p2, *p3 = NULL, *d = NULL;
  unsigned k, jx, bx, nbentries, NbConstraints, Dimension, NbCon, empty;
  unsigned  *Filter = NULL;
  Matrix *Constraints = NULL;
  
  CATCH(any_exception_error) {
    if (Constraints) Matrix_Free(Constraints);
    if (Filter) free(Filter);
    if (d) Polyhedron_Free(d);
    if (p2) Polyhedron_Free(p3);
    RETHROW();
  }
  TRY {
    if (!Pol1 || !Pol2) {
      UNCATCH(any_exception_error);
      return Pol1;
    }
    if (Pol1->Dimension != Pol2->Dimension) {
      errormsg1("DomainSimplify","diffdim","operation on different dimensions");
      UNCATCH(any_exception_error);
      return 0;
    }
    POL_ENSURE_VERTICES(Pol1);
    POL_ENSURE_VERTICES(Pol2);
    if (emptyQ(Pol1)||emptyQ(Pol2)) {
      UNCATCH(any_exception_error);
      return Empty_Polyhedron(Pol1->Dimension);
    }
    
    /* Find the maximum number of constraints over all polyhedron in the  */
    /* polyhedral domain 'Pol2' and store in 'NbCon'.                     */
    NbCon = 0;
    for (p2=Pol2; p2; p2=p2->next)
      if (p2->NbConstraints > NbCon)
	NbCon = p2->NbConstraints;
    
    Dimension = Pol1->Dimension+2;      /* Homogenous Dimension + Status  */
    d = (Polyhedron *)0;
    for (p1=Pol1; p1; p1=p1->next) { 

      /* Filter is an array of integers, each bit in an element of Filter */
      /* array corresponds to a constraint. The bit is marked 1 if the    */
      /* corresponding constraint is non-redundant and is 0 if it is      */
      /* redundant.                                                       */
      
      NbConstraints = p1->NbConstraints;
      nbentries = (NbConstraints + NbCon - 1)/(sizeof(int)*8) + 1;
      
      /* Allocate space for array 'Filter' */
      Filter  = (unsigned *)malloc(nbentries * sizeof(int));
      if(!Filter) {
	errormsg1("DomainSimplify", "outofmem", "out of memory space");
	UNCATCH(any_exception_error);
	return 0;
      }
      
      /* Initialize 'Filter' with zeros */
      SMVector_Init(Filter, nbentries);
      
      /* Filter the constraints of p1 in context to the polyhedra p2(s)   */
      empty = 1;
      for (p2=Pol2; p2; p2=p2->next) {
	
	/* Store the non-redundant constraints in array 'Filter'. With    */
	/* successive loops, the array 'Filter' holds the union of all    */
        /* non-redundant constraints. 'empty' is set to zero if the       */
        /* intersection of two polyhedra is non-empty and Filter is !Null */
   
	if (SimplifyConstraints(p1, p2, Filter, NbMaxRays))
	  empty=0;
      }
      
      if (!empty) {
	
	/* Copy all non-redundant constraints to matrix 'Constraints' */
	Constraints = Matrix_Alloc(NbConstraints,Dimension);
	if(!Constraints) {
	  errormsg1("DomainSimplify", "outofmem", "out of memory space");
	  UNCATCH(any_exception_error);
	  return 0;
	}
	Constraints->NbRows = 0;
	for (k=0, jx=0, bx=MSB; k<NbConstraints; k++) {
	  
	  /* If a bit entry in Filter[jx] is marked 1, copy the correspond- */
	  /* ing constraint in matrix 'Constraints'.                        */
	  if (Filter[jx]&bx) { 
	    Vector_Copy(p1->Constraint[k],
			Constraints->p[Constraints->NbRows],
			Dimension);
	    Constraints->NbRows++;
	  }
	  NEXT(jx,bx);
	}
	
	/* Create the polyhedron 'p3' corresponding to the constraints in   */
	/* matrix 'Constraints'.                                            */
	p3 = Constraints2Polyhedron(Constraints,NbMaxRays);
	Matrix_Free(Constraints), Constraints = NULL;
	
	/* Add polyhedron 'p3' in the domain 'd'. */
	d = AddPolyToDomain (p3, d);
	p3 = NULL;
      }
      free(Filter), Filter = NULL;
    }
  } /* end of TRY */
	
  UNCATCH(any_exception_error);  
  if (!d)
    return Empty_Polyhedron(Pol1->Dimension);
  else
    return d;
} /* DomainSimplify */

/*
 * Return the Union of two polyhedral domains 'Pol1' and Pol2'. The result is
 * a new polyhedral domain.
 */
Polyhedron *DomainUnion(Polyhedron *Pol1,Polyhedron *Pol2,unsigned NbMaxRays) {

  Polyhedron *PolA, *PolEndA, *PolB, *PolEndB, *p1, *p2;
  int Redundant;
  
  if (!Pol1 || !Pol2) return (Polyhedron *) 0;
  if (Pol1->Dimension != Pol2->Dimension) {
    errormsg1("DomainUnion","diffdim","operation on different dimensions");
    return (Polyhedron*) 0;
  }






  /* Copy 'Pol1' to 'PolA' */
  PolA = PolEndA = (Polyhedron *)0;
  for (p1=Pol1; p1; p1=p1->next) {
    
    /* Does any component of 'Pol2' cover 'p1' ? */
    Redundant = 0;
    for (p2=Pol2; p2; p2=p2->next) {
      if (PolyhedronIncludes(p2, p1) ) { /* p2 covers p1 */ 
	Redundant = 1;
	

	break;

      }
    }
    if (!Redundant) {
      
      /* Add 'p1' to 'PolA' */
      if (!PolEndA)
	PolEndA = PolA = Polyhedron_Copy(p1);
      else {
	PolEndA->next = Polyhedron_Copy(p1);
	PolEndA = PolEndA->next;
      }

    }
  }

  /* Copy 'Pol2' to PolB */
  PolB = PolEndB = (Polyhedron *)0;
  for (p2=Pol2; p2; p2=p2->next) {

    /* Does any component of PolA cover 'p2' ? */
    Redundant = 0;
    for (p1=PolA; p1; p1=p1->next) {
      if (PolyhedronIncludes(p1, p2)) { /* p1 covers p2 */
	Redundant = 1;
	break;
      }
    }
    if (!Redundant) {
      
      /* Add 'p2' to 'PolB' */
      if (!PolEndB)
	PolEndB = PolB = Polyhedron_Copy(p2);
      else {
	PolEndB->next = Polyhedron_Copy(p2);
	PolEndB = PolEndB->next;
      }


    }
  }

  if (!PolA) return PolB;
  PolEndA->next = PolB;
  return PolA;
} /* DomainUnion */

/* 
 * Given a polyhedral domain 'Pol', concatenate the lists of rays and lines 
 * of the two (or more) polyhedra in the domain into one combined list, and 
 * find the set of constraints which tightly bound all of those objects. 
 * 'NbMaxConstrs' is the maximum allowed constraints in the new polyhedron. 
 */ 
Polyhedron *DomainConvex(Polyhedron *Pol,unsigned NbMaxConstrs) {
  
  Polyhedron *p, *q, *NewPol = NULL;
  
  CATCH(any_exception_error) {
    if (NewPol) Polyhedron_Free(NewPol);
    RETHROW();
  }
  TRY {
    
    if (!Pol) {
      UNCATCH(any_exception_error);
      return (Polyhedron*) 0;
    }
    
    NewPol = Polyhedron_Copy(Pol);
    for (p=Pol->next; p; p=p->next) {
      q = AddRays(p->Ray[0], p->NbRays, NewPol, NbMaxConstrs);
      Polyhedron_Free(NewPol);
      NewPol = q;
    }
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  
  return NewPol;
} /* DomainConvex */

/*
 * Given polyhedral domains 'Pol1' and 'Pol2', create a new polyhedral 
 * domain which is mathematically the differnce of the two domains. 
 */
Polyhedron *DomainDifference(Polyhedron *Pol1,Polyhedron *Pol2,unsigned NbMaxRays) {

  Polyhedron *p1, *p2, *p3, *d;
  int i;
  
  if (!Pol1 || !Pol2) return (Polyhedron*) 0;
  if (Pol1->Dimension != Pol2->Dimension) {
    errormsg1("DomainDifference", 
	      "diffdim", "operation on different dimensions");
    return (Polyhedron*) 0;
  }
  POL_ENSURE_FACETS(Pol1);
  POL_ENSURE_VERTICES(Pol1);
  POL_ENSURE_FACETS(Pol2);
  POL_ENSURE_VERTICES(Pol2);
  if (emptyQ(Pol1) || emptyQ(Pol2))
    return (Domain_Copy(Pol1));
  d = (Polyhedron *)0;
  for (p2=Pol2; p2; p2=p2->next) {
    for (p1=Pol1; p1; p1=p1->next) {
      for (i=0; i<p2->NbConstraints; i++) {
	
	/* Add the constraint ( -p2->constraint[i] -1) >= 0 in 'p1' */
	/* and create the new polyhedron 'p3'.                      */
	p3 = SubConstraint(p2->Constraint[i], p1, NbMaxRays,0);
	
	/* Add 'p3' in the new domain 'd' */
	d = AddPolyToDomain (p3, d);
	
	/* If the constraint p2->constraint[i][0] is an equality, then  */
	/* add the constraint ( +p2->constraint[i] -1) >= 0  in 'p1' and*/
	/* create the new polyhedron 'p3'.                              */
	
	if( value_notzero_p(p2->Constraint[i][0]) ) /* Inequality */
	  continue;  
	p3 = SubConstraint(p2->Constraint[i], p1, NbMaxRays,1);
	
	/* Add 'p3' in the new domain 'd' */
	d = AddPolyToDomain (p3, d);
      }
    }
    if (p2 != Pol2)
	Domain_Free(Pol1);
    Pol1 = d;
    d = (Polyhedron *)0;
  }
  if (!Pol1)
    return Empty_Polyhedron(Pol2->Dimension);
  else
    return Pol1;
} /* DomainDifference */

/*
 * Given a polyhedral domain 'Pol', convert it to a new polyhedral domain 
 * with dimension expanded to 'align_dimension'. 'NbMaxRays' is the maximum
 * allowed rays in the new polyhedra.
 */
Polyhedron *align_context(Polyhedron *Pol,int align_dimension,int NbMaxRays) {
  
  int i, j, k;
  Polyhedron *p = NULL, **next, *result = NULL;

  CATCH(any_exception_error) {
    if (result) Polyhedron_Free(result);
    RETHROW();
  }
  TRY {
    
    if (!Pol) return Pol;
    if (align_dimension < Pol->Dimension) {
      errormsg1("align_context", "diffdim", "context dimension exceeds data");
      UNCATCH(any_exception_error);
      return Pol;
    }
    if (align_dimension == Pol->Dimension) {
      UNCATCH(any_exception_error);
      return Polyhedron_Copy(Pol);
    }

    /* 'k' is the dimension increment */
    k = align_dimension - Pol->Dimension;
    next = &result;

    /* Expand the dimension of all polyhedron in the polyhedral domain 'Pol' */
    for (; Pol; Pol=Pol->next) {
      int have_cons = !F_ISSET(Pol, POL_VALID) || F_ISSET(Pol, POL_INEQUALITIES);
      int have_rays = !F_ISSET(Pol, POL_VALID) || F_ISSET(Pol, POL_POINTS);
      unsigned NbCons = have_cons ? Pol->NbConstraints : 0;
      unsigned NbRays = have_rays ? Pol->NbRays + k : 0;

      p = Polyhedron_Alloc(align_dimension, NbCons, NbRays);
      if (have_cons) {
	for (i = 0; i < NbCons; ++i) {
	  value_assign(p->Constraint[i][0], Pol->Constraint[i][0]);  /* Status bit */
	  Vector_Copy(Pol->Constraint[i]+1, p->Constraint[i]+k+1, Pol->Dimension+1);
	}
	p->NbEq = Pol->NbEq;
      }

      if (have_rays) {
	for (i = 0; i < k; ++i)
	  value_set_si(p->Ray[i][1+i], 1);			    /* A line */
	for (i = 0; i < Pol->NbRays; ++i) {
	  value_assign(p->Ray[k+i][0], Pol->Ray[i][0]);  	    /* Status bit */
	  Vector_Copy(Pol->Ray[i]+1, p->Ray[i+k]+k+1, Pol->Dimension+1);
	}
	p->NbBid = Pol->NbBid + k;
      }
      p->flags = Pol->flags;
      
      *next = p;
      next = &p->next;
    }
  } /* end of TRY */
  
  UNCATCH(any_exception_error); 
  return result;
} /* align_context */

/*----------------------------------------------------------------------*/
/* Polyhedron *Polyhedron_Scan(D, C, NbMaxRays)                         */
/*       D : Domain to be scanned (single polyhedron only)              */
/*       C : Context domain                                             */
/*       NbMaxRays : Workspace size                                     */
/* Returns a linked list of scan domains, outer loop first              */
/*----------------------------------------------------------------------*/
Polyhedron *Polyhedron_Scan(Polyhedron *D, Polyhedron *C,unsigned NbMaxRays) {
  
  int i, j, dim ;
  Matrix *Mat;
  Polyhedron *C1, *C2, *D1, *D2;
  Polyhedron *res, *last, *tmp;
  
  dim = D->Dimension - C->Dimension;
  res = last = (Polyhedron *) 0;
  if (dim==0) return (Polyhedron *)0;

  assert(!D->next);

  POL_ENSURE_FACETS(D);
  POL_ENSURE_VERTICES(D);
  POL_ENSURE_FACETS(C);
  POL_ENSURE_VERTICES(C);

  /* Allocate space for constraint matrix. */
  Mat   = Matrix_Alloc(D->Dimension, D->Dimension+2);
  if(!Mat) {
    errormsg1("Polyhedron_Scan", "outofmem", "out of memory space");
    return 0;
  }
  C1  = align_context(C,D->Dimension,NbMaxRays);
  if(!C1) {
    return 0;
  }
  /* Vin100, aug 16, 2001:  The context is intersected with D */
  D2 = DomainIntersection( C1, D, NbMaxRays);

  for (i=0; i<dim; i++)
  {
    Vector_Set(Mat->p_Init,0,D2->Dimension*(D2->Dimension + 2));
    for (j=i+1; j<dim; j++) {
      value_set_si(Mat->p[j-i-1][j+1],1);
    }
    Mat->NbRows = dim-i-1;
    D1 = Mat->NbRows ? DomainAddRays(D2, Mat, NbMaxRays) : D2;
    tmp = DomainSimplify(D1, C1, NbMaxRays);
    if (!last) res = last = tmp;
    else { last->next = tmp; last = tmp; }
    C2 = DomainIntersection(C1, D1, NbMaxRays);
    Domain_Free(C1);
    C1 = C2;
    if (Mat->NbRows) Domain_Free(D1);
  }
  Domain_Free(D2);
  Domain_Free(C1);
  Matrix_Free(Mat);
  return res;
} /* Polyhedron_Scan */

/*---------------------------------------------------------------------*/
/* int lower_upper_bounds(pos,P,context,LBp,UBp)                       */
/*    pos : index position of current loop index (1..hdim-1)           */
/*    P: loop domain                                                   */
/*    context : context values for fixed indices                       */
/*              notice that context[hdim] must be 1                    */
/*    LBp, UBp : pointers to resulting bounds                          */
/* returns the flag = (UB_INFINITY, LB_INFINITY)                       */
/*---------------------------------------------------------------------*/
int lower_upper_bounds(int pos,Polyhedron *P,Value *context,Value *LBp,Value *UBp) {
  
  Value LB, UB;
  int flag, i;
  Value n, n1, d, tmp;
  
  POL_ENSURE_FACETS(P);
  POL_ENSURE_VERTICES(P);

  /* Initialize all the 'Value' variables */
  value_init(LB); value_init(UB); value_init(tmp);
  value_init(n); value_init(n1); value_init(d);
  
  value_set_si(LB,0);
  value_set_si(UB,0);
  
  /* Compute Upper Bound and Lower Bound for current loop */
  flag = LB_INFINITY | UB_INFINITY;
  for (i=0; i<P->NbConstraints; i++) {
    value_assign(d,P->Constraint[i][pos]);
    if (value_zero_p(d)) continue;    
    Inner_Product(&context[1],&(P->Constraint[i][1]),P->Dimension+1,&n);
    value_oppose(n,n);
    
    /*---------------------------------------------------*/
    /* Compute n/d        n/d<0              n/d>0       */
    /*---------------------------------------------------*/
    /*  n%d == 0    floor   = n/d      floor   = n/d     */
    /*              ceiling = n/d      ceiling = n/d     */
    /*---------------------------------------------------*/
    /*  n%d != 0    floor   = n/d - 1  floor   = n/d     */
    /*              ceiling = n/d      ceiling = n/d + 1 */
    /*---------------------------------------------------*/

    /* Check to see if constraint is inequality */
    /* if constraint is equality, both upper and lower bounds are fixed */
    if(value_zero_p(P->Constraint[i][0])) {	/* Equality */
      value_modulus(tmp,n,d);
      
      /* if not integer, return 0; */
      if(value_notzero_p(tmp)) {
	value_set_si(*LBp,1);
	value_set_si(*UBp,0);	/* empty loop */
	
	/* Clear all the 'Value' variables */
	value_clear(LB); value_clear(UB); value_clear(tmp);
	value_clear(n); value_clear(n1); value_clear(d);
	return 0;
      }
      value_division(n1,n,d);
      
      /* Upper and Lower bounds found */
      if((flag&LB_INFINITY) || value_gt(n1,LB))
	value_assign(LB,n1);
      if((flag&UB_INFINITY) || value_lt(n1,UB))
	value_assign(UB,n1);
      flag = 0;
    }
    
    if (value_pos_p(d)) {  /* Lower Bound */
      value_modulus(tmp,n,d);
      
      /* n1 = ceiling(n/d) */
      if (value_pos_p(n) && value_notzero_p(tmp)) {
	value_division(n1,n,d);
	value_add_int(n1,n1,1);
      }
      else
	value_division(n1,n,d);
      if (flag&LB_INFINITY) {
	value_assign(LB,n1); 
	flag^=LB_INFINITY; 
      }
      else if(value_gt(n1,LB))
	value_assign(LB,n1);
    }
    
    if (value_neg_p(d)) {   /* Upper Bound */
      value_modulus(tmp,n,d);
      
      /* n1 = floor(n/d) */
      if (value_pos_p(n) && value_notzero_p(tmp)) {
	value_division(n1,n,d);
	value_sub_int(n1,n1,1);
      }
      else
	value_division(n1,n,d);
      
      if (flag&UB_INFINITY) {
	value_assign(UB,n1); 
	flag^=UB_INFINITY; 
      }
      else if (value_lt(n1,UB))
	value_assign(UB, n1);
    }
  }
  if ((flag & LB_INFINITY)==0) value_assign(*LBp,LB);
  if ((flag & UB_INFINITY)==0) value_assign(*UBp,UB);
  
  /* Clear all the 'Value' variables */
  value_clear(LB); value_clear(UB); value_clear(tmp);
  value_clear(n); value_clear(n1); value_clear(d);
  return flag;
} /* lower_upper_bounds */

/*
 *  C = A x B
 */
static void Rays_Mult(Value **A, Matrix *B, Value **C, unsigned NbRays)
{
  int i, j, k;
  unsigned Dimension1, Dimension2;
  Value Sum, tmp;

  value_init(Sum); value_init(tmp);

  CATCH(any_exception_error) {
    value_clear(Sum); value_clear(tmp);
    RETHROW();
  }
  TRY {
    Dimension1 = B->NbRows;
    Dimension2 = B->NbColumns;

    for (i=0; i<NbRays; i++) {
      value_assign(C[i][0],A[i][0]);
      for (j=0; j<Dimension2; j++) {
	value_set_si(Sum,0);
	for (k=0; k<Dimension1; k++) {
	  
	  /* Sum+=A[i][k+1] * B->p[k][j]; */
	  value_addmul(Sum, A[i][k+1], B->p[k][j]);
	}
	value_assign(C[i][j+1],Sum);
      }
      Vector_Gcd(C[i]+1, Dimension2, &tmp);
      if (value_notone_p(tmp))
	  Vector_AntiScale(C[i]+1, C[i]+1, tmp, Dimension2);
    }
  }
  UNCATCH(any_exception_error);
  value_clear(Sum); value_clear(tmp);
}

/*
 *  C = A x B^T
 */
static void Rays_Mult_Transpose(Value **A, Matrix *B, Value **C, 
                                unsigned NbRays)
{
  int i, j, k;
  unsigned Dimension1, Dimension2;
  Value Sum, tmp;

  value_init(Sum); value_init(tmp);

  CATCH(any_exception_error) {
    value_clear(Sum); value_clear(tmp);
    RETHROW();
  }
  TRY {
    Dimension1 = B->NbColumns;
    Dimension2 = B->NbRows;

    for (i=0; i<NbRays; i++) {
      value_assign(C[i][0],A[i][0]);
      for (j=0; j<Dimension2; j++) {
	value_set_si(Sum,0);
	for (k=0; k<Dimension1; k++) {
	  
	  /* Sum+=A[i][k+1] * B->p[j][k]; */
	  value_addmul(Sum, A[i][k+1], B->p[j][k]);
	}
	value_assign(C[i][j+1],Sum);
      }
      Vector_Gcd(C[i]+1, Dimension2, &tmp);
      if (value_notone_p(tmp))
	  Vector_AntiScale(C[i]+1, C[i]+1, tmp, Dimension2);
    }
  }
  UNCATCH(any_exception_error);
  value_clear(Sum); value_clear(tmp);
}

/*
 * Given a polyhedron 'Pol' and a transformation matrix 'Func', return the 
 * polyhedron which when transformed by mapping function 'Func' gives 'Pol'. 
 * 'NbMaxRays' is the maximum number of rays that can be in the ray matrix 
 * of the resulting polyhedron.
 */
Polyhedron *Polyhedron_Preimage(Polyhedron *Pol,Matrix *Func,unsigned NbMaxRays) {

  Matrix *Constraints = NULL;
  Polyhedron *NewPol = NULL;
  unsigned NbConstraints, Dimension1, Dimension2;

  POL_ENSURE_INEQUALITIES(Pol);

  CATCH(any_exception_error) {
    if (Constraints) Matrix_Free(Constraints);
    if (NewPol) Polyhedron_Free(NewPol);
    RETHROW();
  }
  TRY {
    
    NbConstraints = Pol->NbConstraints;
    Dimension1    = Pol->Dimension+1;	/* Homogeneous Dimension */
    Dimension2    = Func->NbColumns;	/* Homogeneous Dimension */
    if (Dimension1!=(Func->NbRows)) {
      errormsg1("Polyhedron_Preimage", "dimincomp", "incompatable dimensions");
      UNCATCH(any_exception_error);
      return Empty_Polyhedron(Dimension2-1);
    }
    
    /*            Dim1           Dim2            Dim2
	          __k__          __j__           __j__	
	    NbCon |   |  X   Dim1|   |  =  NbCon |   |
	      i   |___|       k  |___|       i   |___|
	    Pol->Constraints Function        Constraints
    */
  
    /* Allocate space for the resulting constraint matrix */
    Constraints = Matrix_Alloc(NbConstraints, Dimension2+1);
    if (!Constraints) { 
      errormsg1("Polyhedron_Preimage", "outofmem", "out of memory space\n");
      Pol_status = 1;
      UNCATCH(any_exception_error);
      return 0;
    }
    
    /* The new constraint matrix is the product of constraint matrix of the */
    /* polyhedron and the function matrix.                                  */
    Rays_Mult(Pol->Constraint, Func, Constraints->p, NbConstraints);
    NewPol = Constraints2Polyhedron(Constraints, NbMaxRays);
    Matrix_Free(Constraints), Constraints = NULL;
    
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  
  return NewPol;
} /* Polyhedron_Preimage */

/*
 * Given a polyhedral domain 'Pol' and a transformation matrix 'Func', return 
 * the polyhedral domain which when transformed by mapping function 'Func' 
 * gives 'Pol'. 'NbMaxRays' is the maximum number of rays that can be in the 
 * ray matrix of the resulting domain.
 */
Polyhedron *DomainPreimage(Polyhedron *Pol,Matrix *Func,unsigned NbMaxRays) {
  
  Polyhedron *p, *q, *d = NULL;
  
  CATCH(any_exception_error) {
    if (d) Polyhedron_Free(d);
    RETHROW();
  }
  TRY {
    if (!Pol || !Func) {
      UNCATCH(any_exception_error);
      return (Polyhedron *) 0;
    }
    d = (Polyhedron *) 0;
    for (p=Pol; p; p=p->next) {
      q = Polyhedron_Preimage(p, Func, NbMaxRays);
      d = AddPolyToDomain (q, d);
    } 
  } /* end of TRY */
  UNCATCH(any_exception_error);
  return d;
} /* DomainPreimage */

/*
 * Transform a polyhedron 'Pol' into another polyhedron according to a given
 * affine mapping function 'Func'. 'NbMaxConstrs' is the maximum number of 
 * constraints that can be in the constraint matrix of the new polyhedron. 
 */
Polyhedron *Polyhedron_Image(Polyhedron *Pol, Matrix *Func,unsigned NbMaxConstrs) {
  
  Matrix *Rays = NULL;
  Polyhedron *NewPol = NULL;
  unsigned NbRays, Dimension1, Dimension2;
  
  POL_ENSURE_FACETS(Pol);
  POL_ENSURE_VERTICES(Pol);

  CATCH(any_exception_error) {
    if (Rays) Matrix_Free(Rays);
    if (NewPol) Polyhedron_Free(NewPol);
    RETHROW();
  }
  TRY {
  
    NbRays     = Pol->NbRays;
    Dimension1 = Pol->Dimension+1;	/* Homogeneous Dimension */
    Dimension2 = Func->NbRows;		/* Homogeneous Dimension */
    if (Dimension1!=Func->NbColumns) {
      errormsg1("Polyhedron_Image", "dimincomp", "incompatible dimensions");
      UNCATCH(any_exception_error);
      return Empty_Polyhedron(Dimension2-1);
    }
    
    /*   
        Dim1     /      Dim1  \Transpose      Dim2
        __k__    |      __k__ |              __j__
  NbRays|   |  X | Dim2 |   | |     =  NbRays|   |
    i   |___|    |   j  |___| |          i   |___|
     Pol->Rays  \       Func /               Rays

    */

    if (Dimension1 == Dimension2) {
	Matrix *M, *M2;
	int ok;
	M = Matrix_Copy(Func);
	M2 = Matrix_Alloc(Dimension2, Dimension1);
	if (!M2) {
	  errormsg1("Polyhedron_Image", "outofmem", "out of memory space\n");
	  UNCATCH(any_exception_error);
	  return 0;
	}

	ok = Matrix_Inverse(M, M2);
	Matrix_Free(M);
	if (ok) {
	    NewPol = Polyhedron_Alloc(Pol->Dimension, Pol->NbConstraints,
				      Pol->NbRays);
	    if (!NewPol) {
	      errormsg1("Polyhedron_Image", "outofmem", 
			"out of memory space\n");
	      UNCATCH(any_exception_error);
	      return 0;
	    }
	    Rays_Mult_Transpose(Pol->Ray, Func, NewPol->Ray, NbRays);
	    Rays_Mult(Pol->Constraint, M2, NewPol->Constraint, 
		      Pol->NbConstraints);
	    NewPol->NbEq = Pol->NbEq;
	    NewPol->NbBid = Pol->NbBid;
	    if (NewPol->NbEq)
	      Gauss4(NewPol->Constraint, NewPol->NbEq, NewPol->NbConstraints,
		     NewPol->Dimension+1);
	    if (NewPol->NbBid)
	      Gauss4(NewPol->Ray, NewPol->NbBid, NewPol->NbRays,
		     NewPol->Dimension+1);
	}
	Matrix_Free(M2);
    }
    
    if (!NewPol) {
	/* Allocate space for the resulting ray matrix */
	Rays = Matrix_Alloc(NbRays, Dimension2+1);
	if (!Rays) {
	  errormsg1("Polyhedron_Image", "outofmem", "out of memory space\n");
	  UNCATCH(any_exception_error);
	  return 0;
	}
	
	/* The new ray space is the product of ray matrix of the polyhedron and */
	/* the transpose matrix of the mapping function.                        */
	Rays_Mult_Transpose(Pol->Ray, Func, Rays->p, NbRays);
	NewPol = Rays2Polyhedron(Rays, NbMaxConstrs);
	Matrix_Free(Rays), Rays = NULL;
    }
    
  } /* end of TRY */

  UNCATCH(any_exception_error);
  return NewPol;
} /* Polyhedron_Image */

/* 
 *Transform a polyhedral domain 'Pol' into another domain according to a given
 * affine mapping function 'Func'. 'NbMaxConstrs' is the maximum number of 
 * constraints that can be in the constraint matrix of the resulting domain. 
 */
Polyhedron *DomainImage(Polyhedron *Pol,Matrix *Func,unsigned NbMaxConstrs) {

  Polyhedron *p, *q, *d = NULL;

  CATCH(any_exception_error) {
    if (d) Polyhedron_Free(d);
    RETHROW();
  }
  TRY {
    
    if (!Pol || !Func) {
      UNCATCH(any_exception_error);
      return (Polyhedron *) 0;
    }
    d = (Polyhedron *) 0;
    for (p=Pol; p; p=p->next) { 
      q = Polyhedron_Image(p, Func, NbMaxConstrs);
      d = AddPolyToDomain (q, d);
    }
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  
  return d;
} /* DomainImage */

/* 
 * Given a polyhedron 'Pol' and an affine cost function 'Cost', compute the 
 * maximum and minimum value of the function over set of points representing
 * polyhedron. 
 * Note: If Polyhedron 'Pol' is empty, then there is no feasible solution. 
 * Otherwise, if there is a bidirectional ray with Sum[cost(i)*ray(i)] != 0 or
 * a unidirectional ray with Sum[cost(i)*ray(i)] >0, then the maximum is un-
 * bounded else the finite optimal solution occurs at one of the vertices of
 * the polyhderon. 
 */
Interval *DomainCost(Polyhedron *Pol,Value *Cost) {
  
  int i, j, NbRay, Dim;
  Value *p1, *p2, p3, d, status;
  Value tmp1, tmp2, tmp3;
  Value **Ray;
  Interval *I = NULL;

  value_init(p3); value_init(d); value_init(status);
  value_init(tmp1); value_init(tmp2); value_init(tmp3);

  POL_ENSURE_FACETS(Pol);
  POL_ENSURE_VERTICES(Pol);

  CATCH(any_exception_error) {
    if (I) free(I);
    RETHROW();
    value_clear(p3); value_clear(d); value_clear(status);
    value_clear(tmp1); value_clear(tmp2); value_clear(tmp3);
  }
  TRY {
    
    Ray = Pol->Ray;
    NbRay = Pol->NbRays;
    Dim = Pol->Dimension+1;		/* Homogenous Dimension */
    I = (Interval *) malloc (sizeof(Interval));
    if (!I) {
      errormsg1("DomainCost", "outofmem", "out of memory space\n");
      UNCATCH(any_exception_error);
      value_clear(p3); value_clear(d); value_clear(status);
      value_clear(tmp1); value_clear(tmp2); value_clear(tmp3);
      return 0;
    }
    
    /* The maximum and minimum values of the cost function over polyhedral  */
    /* domain is stored in 'I'. I->MaxN and I->MaxD store the numerator and */
    /* denominator of the maximum value. Likewise,I->MinN and I->MinD store */
    /* the numerator and denominator of the minimum value. I->MaxI and      */
    /* I->MinI store the ray indices corresponding to the max and min values*/
    /* of the function.                                                     */
    
    value_set_si(I->MaxN,-1);
    value_set_si(I->MaxD,0);	      /* Actual cost is MaxN/MaxD */
    I->MaxI = -1;
    value_set_si(I->MinN,1);
    value_set_si(I->MinD,0);
    I->MinI = -1;
    
    /* Compute the cost of each ray[i] */
    for (i=0; i<NbRay; i++) {
      p1 = Ray[i];
      value_assign(status, *p1);
      p1++;
      p2 = Cost;
      
      /* p3 = *p1++ * *p2++; */
      value_multiply(p3,*p1,*p2);
      p1++; p2++;
      for (j=1; j<Dim; j++) {
	value_multiply(tmp1,*p1,*p2);
	
	/* p3 += *p1++ * *p2++; */
	value_addto(p3,p3,tmp1);
	p1++; p2++;
      }
      
      /* d = *--p1; */
      p1--;
      value_assign(d,*p1); /* d == 0 for lines and ray, non-zero for vertex */
      value_multiply(tmp1,p3,I->MaxD); 
      value_multiply(tmp2,I->MaxN,d);
      value_set_si(tmp3,1);
      
      /* Compare p3/d with MaxN/MaxD to assign new maximum cost value */
      if (I->MaxI==-1 ||
	  value_gt(tmp1,tmp2) ||
	  (value_eq(tmp1,tmp2) &&
	   value_eq(d,tmp3) && value_ne(I->MaxD,tmp3))) {
	value_assign(I->MaxN,p3);
	value_assign(I->MaxD,d);
	I->MaxI = i;
      }
      value_multiply(tmp1,p3,I->MinD);
      value_multiply(tmp2,I->MinN,d);
      value_set_si(tmp3,1);
      
      /* Compare p3/d with MinN/MinD to assign new minimum cost value */
      if (I->MinI==-1 ||
	  value_lt(tmp1,tmp2) ||
	  (value_eq(tmp1,tmp2) &&
	   value_eq(d,tmp3) && value_ne(I->MinD,tmp3))) {
	value_assign(I->MinN, p3);
	value_assign(I->MinD, d);
	I->MinI = i;
      }
      value_multiply(tmp1,p3,I->MaxD);
      value_set_si(tmp2,0);
      
      /* If there is a line, assign max to +infinity and min to -infinity */
      if (value_zero_p(status)) { /* line , d is 0 */
	if (value_lt(tmp1,tmp2)) {
          value_oppose(I->MaxN,p3);
          value_set_si(I->MaxD,0);
          I->MaxI = i;
        }
	value_multiply(tmp1,p3,I->MinD);
	value_set_si(tmp2,0);

        if (value_gt(tmp1,tmp2)) {
          value_oppose(I->MinN,p3);
          value_set_si(I->MinD,0);
          I->MinI = i;
        }
      }
    }
  } /* end of TRY */
  
  UNCATCH(any_exception_error);
  value_clear(p3); value_clear(d); value_clear(status);
  value_clear(tmp1); value_clear(tmp2); value_clear(tmp3);
  return I;
} /* DomainCost */

/* 
 * Add constraints pointed by 'Mat' to each and every polyhedron in the 
 * polyhedral domain 'Pol'. 'NbMaxRays' is maximum allowed rays in the ray 
 * matrix of a polyhedron.
 */
Polyhedron *DomainAddConstraints(Polyhedron *Pol,Matrix *Mat,unsigned NbMaxRays) {

  Polyhedron *PolA, *PolEndA, *p1, *p2, *p3;
  int Redundant;
  
  if (!Pol) return (Polyhedron*) 0;
  if (!Mat) return Pol;
  if (Pol->Dimension != Mat->NbColumns-2) {
    errormsg1("DomainAddConstraints",
	      "diffdim", "operation on different dimensions");
    return (Polyhedron*) 0;
  }
  
  /* Copy 'Pol' to 'PolA' */
  PolA = PolEndA = (Polyhedron *)0;
  for (p1=Pol; p1; p1=p1->next) {
    p3 = AddConstraints(Mat->p_Init, Mat->NbRows, p1, NbMaxRays);
    
    /* Does any component of 'PolA' cover 'p3' */
    Redundant = 0;
    for (p2=PolA; p2; p2=p2->next) {
      if (PolyhedronIncludes(p2, p3)) { /* 'p2' covers 'p3' */
	Redundant = 1;
	break;
      }
    }
    
    /* If the new polyhedron 'p3' is not redundant, add it to the domain */
    if (Redundant)
      Polyhedron_Free(p3);
    else { 
      if (!PolEndA)
	PolEndA = PolA = p3;
      else {
	PolEndA->next = p3;
	PolEndA = PolEndA->next;
      }
    }
  }
  return PolA;
} /* DomainAddConstraints */


/* 
 * Computes the disjoint union of a union of polyhedra.
 * If flag = 0 the result is such that there are no intersections
 *                   between the resulting polyhedra,
 * if flag = 1 it computes a joint union, the resulting polyhedra are
 *                   adjacent (they have their facets in common).
 *
 * WARNING: if all polyhedra are not of same geometrical dimension
 *          duplicates may appear.
 */
Polyhedron *Disjoint_Domain( Polyhedron *P, int flag, unsigned NbMaxRays )
{
	Polyhedron *lP, *tmp, *Result, *lR, *prec, *reste;
	Polyhedron *p1, *p2, *p3, *Pol1, *dx, *d1, *d2, *pi, *newpi;
	int i;

	if( flag!=0 && flag!=1 )
	{
		errormsg1("Disjoint_Domain",
			"invalidarg", "flag should be equal to 0 or 1");
		return (Polyhedron*) 0;
	}
	if(!P) return (Polyhedron*) 0;
	if(!P->next) return Polyhedron_Copy(P);

	Result = (Polyhedron *)0;

	for(lP=P;lP;lP=lP->next)
	{
		reste = Polyhedron_Copy(lP);
		prec = (Polyhedron *)0; /* preceeding lR */
		/* Intersection with each polyhedron of the current Result */
		lR=Result;
		while( lR && reste )
		{
			/* dx = DomainIntersection(reste,lR->P,WS); */
			dx = (Polyhedron *)0;
			for( p1=reste; p1; p1=p1->next )
			{
				p3 = AddConstraints(lR->Constraint[0], lR->NbConstraints, p1,
						NbMaxRays);
				dx = AddPolyToDomain(p3,dx);
			}

			/* if empty intersection, continue */
			if(!dx)
			{	prec = lR;
				lR=lR->next;
				continue;
			}
			if (emptyQ(dx)) {	
				Domain_Free(dx);
				prec = lR;
				lR=lR->next;
				continue;
   		}

			/* intersection is not empty, we need to compute the differences */
			/* between the intersection and the two polyhedra, such that the */
			/* results are disjoint unions (according to flag)               */
			/* d1 = reste \ P = DomainDifference(reste,lR->P,WS);	*/
			/* d2 = P \ reste = DomainDifference(lR->P,reste,WS); */

			/* compute d1 */
			d1 = (Polyhedron *)0;
			for (p1=reste; p1; p1=p1->next)
			{
				pi = p1;
				for (i=0; i<P->NbConstraints && pi ; i++)
				{

					/* Add the constraint ( -P->constraint[i] [-1 if flag=0]) >= 0 in 'p1' */
					/* and create the new polyhedron 'p3'.                      */
					p3 = SubConstraint(P->Constraint[i], pi, NbMaxRays,2*flag);
					/* Add 'p3' in the new domain 'd1' */
					d1 = AddPolyToDomain (p3, d1);

					/* If the constraint P->constraint[i][0] is an equality, then add   */
					/* the constraint ( +P->constraint[i] [-1 if flag=0]) >= 0  in 'pi' */
					/* and create the new polyhedron 'p3'.                              */
					if( value_zero_p(P->Constraint[i][0]) ) /* Inequality */
					{
						p3 = SubConstraint(P->Constraint[i], pi, NbMaxRays,1+2*flag);
						/* Add 'p3' in the new domain 'd1' */
						d1 = AddPolyToDomain (p3, d1);

						/* newpi : add constraint P->constraint[i]==0 to pi */
						newpi = AddConstraints( P->Constraint[i], 1, pi, NbMaxRays);
					}
					else
					{
						/* newpi : add constraint +P->constraint[i] >= 0 in pi */
						newpi = SubConstraint(P->Constraint[i], pi, NbMaxRays,3);
					}
					if( newpi && emptyQ( newpi ) )
					{
						Domain_Free( newpi );
						newpi = (Polyhedron *)0;
					}
					if( pi != p1 )
						Domain_Free( pi );
					pi = newpi;
				}
				if( pi != p1 )
					Domain_Free( pi );
			}

			/* and now d2 */
			Pol1 = Polyhedron_Copy( lR );
			for (p2=reste; p2; p2=p2->next)
			{
				d2 = (Polyhedron *)0;
				for (p1=Pol1; p1; p1=p1->next)
				{
					pi = p1;
					for (i=0; i<p2->NbConstraints && pi ; i++)
					{

						/* Add the constraint ( -p2->constraint[i] [-1 if flag=0]) >= 0 in 'pi' */
						/* and create the new polyhedron 'p3'.                      */
						p3 = SubConstraint(p2->Constraint[i], pi, NbMaxRays,2*flag);
						/* Add 'p3' in the new domain 'd2' */
						d2 = AddPolyToDomain (p3, d2);

						/* If the constraint p2->constraint[i][0] is an equality, then add   */
						/* the constraint ( +p2->constraint[i] [-1 if flag=0]) >= 0  in 'pi' */
						/* and create the new polyhedron 'p3'.                              */
						if( value_zero_p(p2->Constraint[i][0]) ) /* Inequality */
						{
							p3 = SubConstraint(p2->Constraint[i], pi, NbMaxRays,1+2*flag);
							/* Add 'p3' in the new domain 'd2' */
							d2 = AddPolyToDomain (p3, d2);

							/* newpi : add constraint p2->constraint[i]==0 to pi */
							newpi = AddConstraints( p2->Constraint[i], 1, pi, NbMaxRays);
						}
						else
						{
							/* newpi : add constraint +p2->constraint[i] >= 0 in pi */
							newpi = SubConstraint(p2->Constraint[i], pi, NbMaxRays,3);
						}
						if( newpi && emptyQ( newpi ) )
						{
							Domain_Free( newpi );
							newpi = (Polyhedron *)0;
						}
						if( pi != p1 )
							Domain_Free( pi );
						pi = newpi;
					}
					if( pi && pi!=p1 )
						Domain_Free( pi );
				}
				if( Pol1 )
					Domain_Free( Pol1 );
				Pol1 = d2;
			}
			/* ok, d1 and d2 are computed */

			/* now, replace lR by d2+dx (at least dx is nonempty) and set reste to d1 */
			if( d1 && emptyQ(d1) )
			{
				Domain_Free( d1 );
				d1 = NULL;
			}
			if( d2 && emptyQ(d2) )
			{
				Domain_Free( d2 );
				d2 = NULL;
			}

			/* set reste */
			Domain_Free( reste );
			reste = d1;

			/* add d2 at beginning of Result */
			if( d2 )
			{
				for( tmp=d2 ; tmp->next ; tmp=tmp->next )
						;
				tmp->next = Result;
				Result = d2;
				if( !prec )
					prec = tmp;
			}

			/* add dx at beginning of Result */
			for( tmp=dx ; tmp->next ; tmp=tmp->next )
				;
			tmp->next = Result;
			Result = dx;
			if( !prec )
				prec = tmp;

			/* suppress current lR */
			if( !prec )
				errormsg1( "Disjoint_Domain","internalerror","internal error");
			prec->next = lR->next;
			Polyhedron_Free( lR );
			lR = prec->next;
		} /* end for result */

		  /* if there is something left, add it to Result : */
		if(reste)
		{
			if(emptyQ(reste))
			{
				Domain_Free( reste );
				reste = NULL;
			}
			else
			{
				Polyhedron *tnext;
				for( tmp=reste ; tmp ; tmp=tnext )
				{
					tnext = tmp->next;
					tmp->next = NULL;
					Result = AddPolyToDomain(tmp, Result);
				}
   		}
		}
	}

	return( Result );
}



/* Procedure to print constraint matrix of a polyhedron */
void Polyhedron_PrintConstraints(FILE *Dst,char *Format,Polyhedron *Pol)
{
	int i,j;

	fprintf( Dst, "%d %d\n", Pol->NbConstraints, Pol->Dimension+2 );
	for( i=0 ; i<Pol->NbConstraints ; i++ )
	{
		for( j=0 ; j<Pol->Dimension+2 ; j++ )
			value_print( Dst, Format, Pol->Constraint[i][j] );
		fprintf( Dst, "\n" );
	}

}

/* Procedure to print constraint matrix of a domain */
void Domain_PrintConstraints(FILE *Dst,char *Format,Polyhedron *Pol)
{
    Polyhedron *Q;
    for (Q = Pol; Q; Q = Q->next)
	Polyhedron_PrintConstraints(Dst, Format, Q);
}

static Polyhedron *p_simplify_constraints(Polyhedron *P, Vector *row,
					  Value *g, unsigned MaxRays)
{
    Polyhedron *T, *R = P;
    int len = P->Dimension+2;
    int r;

    /* Also look at equalities.
     * If an equality can be "simplified" then there
     * are no integer solutions anyway and the following loop
     * will add a conflicting constraint
     */
    for (r = 0; r < R->NbConstraints; ++r) {
	if (ConstraintSimplify(R->Constraint[r], row->p, len, g)) {
	    T = R;
	    R = AddConstraints(row->p, 1, R, MaxRays);
	    if (T != P)
		Polyhedron_Free(T);
	    r = -1;
	}
    }
    if (R != P)
	Polyhedron_Free(P);
    return R;
}

/*
 * Replaces constraint a x >= c by x >= ceil(c/a)
 * where "a" is a common factor in the coefficients
 * Destroys P and returns a newly allocated Polyhedron
 * or just returns P in case no changes were made
 */
Polyhedron *DomainConstraintSimplify(Polyhedron *P, unsigned MaxRays)
{
    Polyhedron **prev;
    int len = P->Dimension+2;
    Vector *row = Vector_Alloc(len);
    Value g;
    Polyhedron *R = P, *N;
    value_set_si(row->p[0], 1);
    value_init(g);

    for (prev = &R; P; P = N) {
	Polyhedron *T;
	N = P->next;
	T = p_simplify_constraints(P, row, &g, MaxRays);

	if (emptyQ(T) && prev != &R) {
	    Polyhedron_Free(T);
	    *prev = NULL;
	    continue;
	}

	if (T != P)
	    T->next = N;
	*prev = T;
	prev = &T->next;
    }

    if (R->next && emptyQ(R)) {
	N = R->next;
	Polyhedron_Free(R);
	R = N;
    }

    value_clear(g);
    Vector_Free(row);
    return R;
}
