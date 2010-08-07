/***********************************************************************/
/*                Parametrized polyhedra V4.20                         */
/*                copyright 1995-2000 Vincent Loechner                 */
/*                copyright 1996-1997, Doran Wilde                     */
/*       Permission is granted to copy, use, and distribute            */
/*       for any commercial or noncommercial purpose under the terms   */
/*       of the GNU General Public license, version 2, June 1991       */
/*       (see file : LICENSING).                                       */
/***********************************************************************/

/********************* -----------USER #DEFS-------- ***********************/
/* These are mainly for debug purposes. You shouldn't need to change       */
/* anything for daily usage...                                             */
/***************************************************************************/

/* you may define each macro independently */
/* #define DEBUGPP 	*/
/* #define DEBUGPP3	*/		/* initialization of domain, context, ... */
/* #define DEBUGPP31	*/		/* even more init-domains */
/* #define DEBUGPP32	*/		/* even even more... (Elim_Columns) */
/* #define DEBUGPP4	*/		/* m-faces scan */
/* #define DEBUGPP41	*/		/* inverse Di in scan */
/* #define DEBUGPP5	*/		/* Compute_PDomains */
/********************* ---------END USER #DEFS------ ***********************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifdef DEBUGPP
#include <time.h>
#endif

#include <polylib/polylib.h>

#ifdef __STDC__
static void traite_m_face(Polyhedron *, unsigned int *);
static void scan_m_face(int,int,Polyhedron *,unsigned int *);
#else
static void traite_m_face();
static void scan_m_face();
#endif /* __STDC__ */

/*
 * Return the intersection of two polyhedral domains 'Pol1' and 'Pol2' such 
 * that if the intersection is a polyhedron of lower dimension (a degenerate
 * polyhedron) than the operands, it is discarded from the resulting polyhedra
 * list. 
 */
Polyhedron *PDomainIntersection(Polyhedron *Pol1,Polyhedron *Pol2,unsigned NbMaxRays) {
	
  Polyhedron *p1, *p2, *p3, *d;
  
  if (!Pol1 || !Pol2) return (Polyhedron*) 0;
  if((Pol1->Dimension != Pol2->Dimension) || (Pol1->NbEq != Pol2->NbEq)) {
    fprintf(stderr,
	    "? PDomainIntersection: operation on different dimensions\n");
    return (Polyhedron*) 0;
  }
  
  POL_ENSURE_FACETS(Pol1);
  POL_ENSURE_VERTICES(Pol1);
  POL_ENSURE_FACETS(Pol2);
  POL_ENSURE_VERTICES(Pol2);
 
  d = (Polyhedron *)0;
  for (p1=Pol1; p1; p1=p1->next) {
    for (p2=Pol2; p2; p2=p2->next) {
      p3 = AddConstraints(p2->Constraint[0],
			  p2->NbConstraints,p1,NbMaxRays);
      if (!p3) continue;
      
      /* If the new polyhedron 'p3' has lower dimension, discard it */
      if (p3->NbEq!=Pol1->NbEq)
	Polyhedron_Free(p3) ;
      
      /* Otherwise add it to the new polyhderal domain 'd'. */
      else
	d = AddPolyToDomain(p3,d);
    }
  }
  return d;
} /* PDomainIntersection */

/* 
 * Given polyhderal domains 'Pol1' and 'Pol2', return the difference of the 
 * two domains with a modification that the resulting polyhedra in the new 
 * domain don't have a 1 unit space around cut and the degenrate results are 
 * discarded. 
 */
Polyhedron *PDomainDifference(Polyhedron *Pol1,Polyhedron *Pol2,unsigned NbMaxRays) {
  
  Polyhedron *p1, *p2, *p3, *d;
  int i;
  
  if (!Pol1 || !Pol2)
    return (Polyhedron*) 0;
  if((Pol1->Dimension != Pol2->Dimension) || (Pol1->NbEq != Pol2->NbEq)) {
    fprintf(stderr,
	    "? PDomainDifference: operation on different dimensions\n");
    return (Polyhedron*) 0;
  }

  POL_ENSURE_FACETS(Pol1);
  POL_ENSURE_VERTICES(Pol1);
  POL_ENSURE_FACETS(Pol2);
  POL_ENSURE_VERTICES(Pol2);
 
  d = (Polyhedron *)0;
  for (p2=Pol2; p2; p2=p2->next) {
    for (p1=Pol1; p1; p1=p1->next) {
      for (i=0; i<p2->NbConstraints; i++) {
	
	/* Add the constraint (-p2->Constraint[i]) >= 0 in 'p1' */
	p3 = SubConstraint(p2->Constraint[i],p1,NbMaxRays,2);
	if (!p3) continue;
	
	/* If the new polyhedron 'p3' is empty or is a polyhedron of lower */
	/* dimension, discard it.                                          */
	if (emptyQ(p3) || p3->NbEq!=Pol1->NbEq)
	  Polyhedron_Free(p3);
	
	/* Otherwise add 'p3' to the new polyhderal domain 'd' */
	else
	  d = AddPolyToDomain(p3,d);
      }
    }
    if (p2 != Pol2)
	Domain_Free(Pol1);
    Pol1 = d;
    d = (Polyhedron *)0;
  }
  return Pol1;
} /* PDomainDifference */

/* 
 * Return 1 if matrix 'Mat' is full column ranked, otherwise return 0. 
 */  
static int TestRank(Matrix *Mat) {
  
  int i,j,k;
  Value m1,m2,m3,gcd,tmp;

  /* Initialize all the 'Value' variables */
  value_init(m1); value_init(m2); 
  value_init(m3); value_init(gcd); value_init(tmp);
  
  for(k=0;k<Mat->NbColumns;++k) {
    
    /* If the digonal entry (k,k) is zero, search down the column(k) */
    /* starting from row(k) to find a non-zero entry                 */
    if(value_zero_p(Mat->p[k][k])) {
      for(j=k+1;j<Mat->NbRows;++j) {
	
	/* If a non-zero entry (j,k) is found */
	if(value_notzero_p(Mat->p[j][k])) {
	  
	  /* Exchange row(k) and row(j) */
	  for(i=k;i<Mat->NbColumns;++i) {
	    value_assign(tmp,Mat->p[j][i]);
	    value_assign(Mat->p[j][i],Mat->p[k][i]);
	    value_assign(Mat->p[k][i],tmp);
	  }
	  break;
	}
      }
      
      /* If no non-zero entry is found then the matrix 'Mat' is not full */
      /* ranked. Return zero.                                            */
      if(j>=Mat->NbRows) {
	
	/* Clear all the 'Value' variables */
	value_clear(m1); value_clear(m2); 
	value_clear(m3); value_clear(gcd); value_clear(tmp);
	return 0;
      }	
    }
    
    /* Now Mat[k][k] is the pivot element */
    for(j=k+1;j<Mat->NbRows;++j) {
      
      /* Make every other entry (below row(k)) in column(k) zero */
      Gcd(Mat->p[j][k],Mat->p[k][k],&gcd);
      for(i=k+1;i<Mat->NbColumns;++i) {
	
	/* pour tous les indices i > k */
	value_multiply(m1,Mat->p[j][i],Mat->p[k][k]);
	value_multiply(m2,Mat->p[j][k],Mat->p[k][i]);
	value_subtract(m3,m1,m2);
	value_division(Mat->p[j][i],m3,gcd);
      }
    }   
  }

  /* Clear all the 'Value' variables */
  value_clear(m1); value_clear(m2); 
  value_clear(m3); value_clear(gcd); value_clear(tmp);
  
  /* The matrix 'Mat' is full ranked, return 1 */
  return 1;
} /* TestRank */

/*
 * The Saturation matrix is defined to be an integer (int type) matrix. It is
 * a boolean matrix which has a row for every constraint and a column for 
 * every line or ray. The bits in the binary format of each integer in the 
 * saturation matrix stores the information whether the corresponding constr-
 * aint is saturated by ray(line) or not. 
 */
typedef struct {
  unsigned int NbRows;
  unsigned int NbColumns;
  unsigned int **p;
  unsigned int *p_init;
} SatMatrix; 

static SatMatrix *SMAlloc(int rows,int cols) {

  unsigned int **q, *p;
  int i;

  SatMatrix *result = (SatMatrix *)malloc(sizeof(SatMatrix));
  assert (result != NULL);

  result->NbRows = rows;
  result->NbColumns = cols;
  result->p = q = (unsigned int **)malloc(rows * sizeof(unsigned int *));
  assert (result->p != NULL);  
  result->p_init = p = (unsigned int *)malloc(rows * cols * sizeof(unsigned int));
  assert (result->p_init != NULL);  
  
  for (i=0;i<rows;i++) {
    *q++ = p;
    p += cols;
  }
  
  return result;
} /* SMAlloc */

static void SMPrint (SatMatrix *matrix) {
  
  unsigned int *p;
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


static void SMFree (SatMatrix *matrix) {
  
  free ((char *) matrix->p_init);
  free ((char *) matrix->p);
  free ((char *) matrix);
  return;
} /* SMFree */

/* -------------------------------------------------------------------------
 * Shared Global Variables:
 * Used by procedures: Find_m_face, scan_m_face, Poly2Sat, traite_m_face,
 *                     count_sat 
 * -------------------------------------------------------------------------
 */
static int m;			/* number of parameters */
static int m_dim;		/* dimension of m-face */
static int n;			/* dimension (not including parameters) */
static int ws;			/* Working Space size */
static int nr;			/* (NbRays-1)/32 + 1 */

static Polyhedron *CEqualities;/* Equalities in the context */
static SatMatrix   *Sat;       /* Saturation Matrix (row=constraint, col=ray)*/
static int *egalite;	       /* Bool vector marking constraints in m-face  */
static Matrix *Xi, *Pi;	       /* Xi and Pi */
static Matrix *PiTest;	       /* Matrix used to test if Pi is full ranked? */
static Matrix *CTest;
static Matrix *PiInv;	       /* Matrix inverse Pi, with the last col of   */
			       /* each line = denominator of the line       */
static Matrix *RaysDi;	       /* Constraint matrix for computing Di */

static int KD;			 /* Flag : keep the full domains in memory ? */
				 /* 1 = yes; 0 = no, keep constraints only   */

static int nbPV;		  /* The number of parameterized vertices */
static Param_Vertices *PV_Result; /* List of parameterized vertices */
static Param_Domain *PDomains;    /* List of domains. */

#ifdef DEBUGPP
static int nbfaces;
#endif

/*
 * Add the constraints from the context polyhedron 'CEqualities' to the 
 * constraints of polyhedra in the polyhedral domain 'D' and return the new
 * polyhedral domain. Polyhedral domain 'D' is subsequently deleted from memory
 */
static Polyhedron *Add_CEqualities(Polyhedron *D) {
  
  Polyhedron *d,*r,*tmp;

  if(!CEqualities)
    return D;
  else {
    if(!D || emptyQ(D)) {
      if(D)
	Domain_Free(D);
      return(Polyhedron_Copy(CEqualities));
    }
    r = AddConstraints(D->Constraint[0],D->NbConstraints,
		       CEqualities,ws);
    tmp = r;
    for(d=D->next;d;d=d->next) {
      tmp->next = AddConstraints(d->Constraint[0],d->NbConstraints,
				 CEqualities,ws);
      tmp = tmp->next;
    }
    Domain_Free(D);
    return(r);
  }
} /* Add_CEqualities */

/*----------------------------------------------------------------------*/
/* traite_m_face                                                        */
/*       Given an m-face, compute the parameterized vertex              */
/*----------------------------------------------------------------------*/
static void traite_m_face(Polyhedron *D,unsigned int *mf) {
     /* D  - The entire domain */
     /* mf - Bit vector marking the lines/rays in the m-face */

  Matrix *Si;				/* Solution */
  Polyhedron *PDi;		        /* polyhedron Di */
  Param_Vertices *PV;
  int j,k,c,r;
  unsigned kx, bx;

#ifdef DEBUGPP
  ++nbfaces;
#endif
  
  /* Extract  Xi, Pi, and RaysDi from D */
  RaysDi->NbRows = 0;
  for(k=0,c=0,kx=0,bx=MSB;k<D->NbRays;++k) {
    if(mf[kx]&bx) {      /* this ray is in the current m-face */      
      if(c<m+1) {	
	int i;
	
	/* tester si cette nouvelle colonne est lin. indep. des autres */
	/* i.e. si gauss ne donne pas de '0' sur la colonne Pi */
	/* jusqu'a l'indice 'c'                                */
	
	/* construit PiTest */
	for(j=0;j<m+1;++j) {
	  for(i=0;i<c;++i)
	    
	    /* les c premieres colonnes */
	    value_assign(PiTest->p[j][i],Pi->p[j][i]);
	  
	  /* la nouvelle */
	  value_assign(PiTest->p[j][c],D->Ray[k][j+1+n]);
	}
	PiTest->NbColumns = c+1;
	r = TestRank(PiTest);
	if(r /* TestRank(PiTest) */) {
				
	  /* Ok, c'est lin. indep. */
	  for (j=0;j<n;j++)
	    value_assign(Xi->p[j][c],D->Ray[k][j+1]);	/* Xi */
	  for (j=0;j<m;j++)
	    value_assign(Pi->p[j][c],D->Ray[k][j+1+n]);	/* Pi */
	  value_assign(Xi->p[n][c],D->Ray[k][n+m+1]);	/* const */
	  value_assign(Pi->p[m][c],D->Ray[k][n+m+1]);	/* const */
	  c++;
	}
      }
      
      /* Status bit */
      value_assign(RaysDi->p[RaysDi->NbRows][0],D->Ray[k][0]);     
      Vector_Copy(&D->Ray[k][n+1],&RaysDi->p[RaysDi->NbRows][1],(m+1));
      ++RaysDi->NbRows;
    }
    NEXT(kx,bx);
  }
  
#ifdef DEBUGPP41
  fprintf(stderr, "\nRaysDi=\n");
  Matrix_Print(stderr,P_VALUE_FMT,RaysDi);
  if(c < m+1)
    fprintf(stderr, "Invalid ");
  fprintf(stderr, "Pi=\n");
  Matrix_Print(stderr,P_VALUE_FMT,Pi);
#endif
  
#ifdef DEBUGPP4
  if(c < m+1)
    fprintf(stderr,"Eliminated because of no vertex\n");
#endif

  if(c < m+1)	
    return;

  /* RaysDi->numColumns = m+2; */  /* stays the same */

  /*	Xi->NbColumns = m+1;*/	/* VIN100: stays the same. was 'c'! */
  /*	Xi->NbRows = n+1; */ 	/* stays the same */
  /*	Pi->NbColumns = m+1;*/	/* VIN100: stays the same. was 'c'! */
  /*	Pi->NbRows = m+1; */		/* stays the same */
  
#ifdef DEBUGPP4
  fprintf(stderr,"Xi = ");
  Matrix_Print(stderr,P_VALUE_FMT,Xi);
  fprintf(stderr,"Pi = ");
  Matrix_Print(stderr,P_VALUE_FMT,Pi);
#endif
  
  /* (Right) invert Pi if POSSIBLE, if not then next m-face */
  /* Pi is destroyed                                        */
  if(!MatInverse(Pi,PiInv)) {
    
#ifdef DEBUGPP4
    fprintf(stderr, "Eliminated because of no inverse Pi\n");
#endif
    
    return;
  }
  
#ifdef DEBUGPP4
  fprintf(stderr,"FACE GENERATED!\n");
  fprintf(stderr,"PiInv = ");
  Matrix_Print(stderr,P_VALUE_FMT,PiInv);
#endif
  
  /* Compute  Si (now called Ti in the paper) */
  Si = Matrix_Alloc(Xi->NbRows,PiInv->NbColumns);
  rat_prodmat(Si,Xi,PiInv);
  
#ifdef DEBUGPP4
  fprintf(stderr,"Si = ");
  Matrix_Print(stderr,P_VALUE_FMT,Si);
#endif
  
  Si->NbRows--;      /* throw out the last row = 0 ... 0 1 */
  
  /* Copy all of that into the PV structure */
  PV = (Param_Vertices *) malloc(sizeof(Param_Vertices));
  PV->next = PV_Result;
  PV->Vertex = Si;
  PV->Domain = NULL;
  PV_Result = PV;
  nbPV++;         /* increment vertex count */
  
  /* Ok... now compute the parameter domain */
  PDi = Rays2Polyhedron(RaysDi,ws);
  
#ifdef DEBUGPP3
  fprintf(stderr,"RaysDi = ");
  Matrix_Print(stderr,P_VALUE_FMT,RaysDi);
  fprintf(stderr,"PDi = ");
  Polyhedron_Print(stderr,P_VALUE_FMT,PDi);
#endif
  
  if(KD==0) {
    
    /* Add the equalities again to the domain */
    PDi = Add_CEqualities(PDi);
    PV->Domain = Polyhedron2Constraints(PDi);
    Polyhedron_Free(PDi);
  }
  else {
    Param_Domain *PD;
    PD = (Param_Domain *) malloc(sizeof(Param_Domain));
    PD->Domain = PDi;
    PD->F = NULL;
    PD->next = PDomains;
    PDomains = PD;
  }
  return;
} /* traite_m_face */

/*----------------------------------------------------------------------*/
/* count_sat                                                            */
/*      count the number of saturated rays in the bit vector mf         */
/*      Uses nr from global area                                        */
/*----------------------------------------------------------------------*/
int cntbit[256] = {				/* counts for 8 bits */
0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,

1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,

1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,

2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8 };

static int count_sat (unsigned int *mf) {
  
  register unsigned int i, tmp, cnt=0;
  
  for (i=0; i<nr; i++) {
    tmp = mf[i];
    cnt = cnt
      + cntbit[ tmp & 0xff ]
      + cntbit[ (tmp>>8) & 0xff ]
      + cntbit[ (tmp>>16) & 0xff ]
      + cntbit[ (tmp>>24) & 0xff ]
      ;
  }
  return cnt;
} /* count_sat */

/*----------------------------------------------------------------------*/
/* let D + E + L be the dimension of the polyhedron                     */
/* D = Dimension of polytope (ray space)                                */
/* L = Dimension of Lineality space (number of lines, bid)              */
/* E = Dimension of Affine hull (number of equations)                   */
/* n = number of data indices                                           */
/* m = number of parameters                                             */
/* full domain:                                                         */
/*     n + m = D + E + L                                                */
/* projected domains:                                                   */
/*     m = D_m + E_m + L_m                                              */
/*     n = D_n + E_n + L_n                                              */
/* What dimension M-face, when projected onto parameter space,          */
/* will give an L_m-face?                                               */
/* What are the conditions?                                             */
/*   - at least one vertex                                              */
/*   - number of rays >= D_m+1 after removal of redundants              */
/*                                                                      */
/* dim of face    nb saturated constraints   nb saturated lines,rays    */
/* -----------    ------------------------   -----------------------    */
/* (0+L)-face     all E eqns + >=D ineq      all L lines + 1 ray        */
/* (M+L)-face     all E eqns + >=(D-M) ineq  all L lines + >=(M+1) rays */
/* (D+L)-face     all E eqns + 0 ineq        all L lines + >=(D+1) rays */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/* scan_m_face                                                          */
/*      pos : the next candidate constraint position                    */
/*    nb-un : number of saturated constraints needed to finish a face   */
/*        D : the source polyhedron (context included )                 */
/*       mf : bit-array marking rays which are saturated so far         */
/* From Global area:                                                    */
/* ----------------                                                     */
/*        n : number of data indices                                    */
/*        m : number of parameters                                      */
/*  egalite : boolean vector marking saturated constraints in m-face    */
/*      Sat : Saturation Matrix (row=constraints, col=rays)             */
/*       ws : working space size                                        */
/*       nr : (NbRays-1)/32 + 1                                         */
/*                                                                      */
/* Recursive function to find the rays and vertices of each m-face      */
/*----------------------------------------------------------------------*/
static void scan_m_face(int pos,int nb_un,Polyhedron *D,unsigned int *mf) {
  /* pos   - the next candidate constraint position */
  /* nb_un - the number of constraints needed to finish a face */
  /* D     - the source polyhedron */
  /* mf    - (bit vector) marks rays that are saturated so far */

  unsigned int *new_mf;

#ifdef DEBUGPP4
  fprintf(stderr,"Start scan_m_face(pos=%d, nb_un=%d, n=%d, m=%d\n",
	  pos,nb_un,n,m);
  fprintf(stderr,"mf = ");
  {
    int i;
    for(i=0;i<nr;i++)
      fprintf(stderr,"%08X", mf[i]);
    fprintf(stderr,"\nequality = [");
    for(i=0;i<D->NbConstraints;i++)
      fprintf(stderr," %1d",egalite[i]);
    fprintf(stderr,"]\n");
  }
#endif

  if(nb_un == 0) {	/* Base case */   
    int i,j;
    
    /*********** ELIMINATION OF REDUNDANT FACES ***********/
    /* if all these vertices also verify a previous constraint */
    /* which is NOT already selected, we eliminate this face */
    /* This keeps the lexicographically greatest selection */
	for(i=0;i<pos-1;i++)
	{
		if(egalite[i])
			continue;       /* already selected */
      
		/* if Sat[i] & mf == mf then it's redundant */
		for (j=0;j<nr;j++)
		{
	
#ifdef DEBUGPP4
			fprintf(stderr,"mf=%08X Sat[%d]=%08X &=%08X\n",mf[j],i,Sat->p[i][j],
				(mf[j] & Sat->p[i][j]) );
#endif

			if (((mf[j]) & (Sat->p[i][j])) != mf[j])
				break;	/* it's not redundant */
		}
      
#ifdef DEBUGPP4
		if (j==nr) fprintf(stderr, "Redundant with constraint %d\n", i);
#endif
    
		if (j==nr) return;	/* it is redundant */
	}
    /********* END OF ELIMINATION OF DEGENERATE FACES *********/    
    /* if we haven't found a constraint verified by all */
    /* the rays, its OK, it's a new face. */
    traite_m_face(D,mf);
    return;
  }

  /* See if there are enough constraints left to finish */
  if((pos+nb_un)>D->NbConstraints) return;
  
  /* Recurring part of the procedure */
  /* Add the pos'th constraint, compute new saturation vector */
  {	
    int k;
    new_mf  = (unsigned int *)malloc(nr*sizeof(unsigned int));
    for (k=0; k<nr; k++)
      new_mf[k] = mf[k] & Sat->p[pos][k];
  }
#ifdef DEBUGPP4
fprintf(stderr,"new_mf = ");
 { 
   int i;
   for(i=0;i<nr;i++) {
     fprintf(stderr,"%08X", new_mf[i]);
   }
   fprintf(stderr,"\ncount(new_mf) = %d\n",count_sat(new_mf));
 }
#endif

  {
 	int c;
	c = count_sat(new_mf);
	/* optimization : at least m_dim+1 rays must be saturated to add this constraint */
	if (c>m_dim )
	{
		int redundant = 0;

                egalite[pos]=1;		/* Try it with the pos-th constraint */

		/* If this constraint does not change anything,
		 * it is redundant with respect to the selected
		 * equalities and the remaining inequalities.
		 * Check whether it is redundant with respect
		 * to just the selected equalities.
		 */
		if( c==count_sat(mf) ) {
		    int i, c, j;

		    for (i = 0, c = 0; i < D->NbConstraints; ++i) {
			if (egalite[i] == 0 || egalite[i] == -1)
			    continue;
			for (j = 0; j < D->Dimension+1; ++j)
			    value_assign(CTest->p[j][c], 
					 D->Constraint[i][j+1]);
			++c;
		    }
		    CTest->NbColumns = c;
#ifdef DEBUGPP41
		    Matrix_Print(stderr,P_VALUE_FMT,CTest);
#endif
		    redundant = !TestRank(CTest);
		}

		/* Do not decrement nb_un if equality is redundant. */
		if( redundant )
		{
		   egalite[pos]=-1;	/* Don't use in further redundance test
					 */
		   scan_m_face(pos+1,nb_un,D,new_mf);
		}
		else
		{
		   scan_m_face(pos+1,nb_un-1,D,new_mf);
		}
	}
  }
 free(new_mf);
 egalite[pos]=0;		/* Try it without the pos-th constraint */
 if ((pos+nb_un)>=D->NbConstraints) return;
 scan_m_face(pos+1,nb_un,D,mf);
 return;
} /* scan_m_face */

/*
 * Create a saturation matrix with rows correspond to the constraints and 
 * columns correspond to the rays of the polyhedron 'Pol'. Global variable 
 * 'nr' is set in the function.             
 */
static SatMatrix *Poly2Sat(Polyhedron *Pol,unsigned int **L) { 
  
  SatMatrix *Sat;
  int i, j, k, kx;
  unsigned int *Temp;
  Value *p1, *p2, p3,tmp;
  unsigned Dimension, NbRay, NbCon, bx;
  
  /* Initialize all the 'Value' variables */
  value_init(p3); value_init(tmp);
  
  NbRay = Pol->NbRays;
  NbCon = Pol->NbConstraints;
  Dimension = Pol->Dimension+1;                /* Homogeneous Dimension */
  
  /* Build the Sat matrix */
  nr      = (NbRay - 1)/(sizeof(int)*8) + 1;   /* Set globally */
  Sat     = SMAlloc(NbCon,nr);
  Temp     = (unsigned int *)malloc(nr*sizeof(unsigned int));
  memset(Sat->p_init,0,nr*NbCon*sizeof(int));
  memset(Temp,0,nr*sizeof(unsigned int));
  kx=0; bx=MSB;
  for (k=0; k<NbRay; k++) { 
    for (i=0; i<NbCon; i++) {
      p1 = &Pol->Constraint[i][1];
      p2 = &Pol->Ray[k][1];
      value_set_si(p3,0);
      for (j=0;j<Dimension;j++) {
	value_multiply(tmp,*p1,*p2);
	value_addto(p3,p3,tmp);
	p1++; p2++;
      }
      if (value_zero_p(p3))
	Sat->p[i][kx]|=bx;
    }
    Temp[kx] |= bx;
    NEXT(kx, bx);
  }
  
  /* Set 'L' to an array containing ones in every bit position of its */
  /* elements.                                                        */
  *L = Temp;
  
  /* Clear all the 'Value' variables */
  value_clear(p3); value_clear(tmp);
  
  return Sat;
} /* Poly2Sat */

/*
 * Create a parametrized polyhedron with zero parameters. This function was 
 * first written by Xavier Redon, and was later modified by others.
 */
Param_Polyhedron *GenParamPolyhedron(Polyhedron *Pol) {
  
  Param_Polyhedron *result;
  int nbRows, nbColumns;
  int i, size, rays;
  
  nbRows=Pol->NbRays;
  nbColumns=Pol->Dimension+2;
  
  /* Count the number of rays */
  for(i=0, rays=0; i<nbRows; i++)
    if(value_notone_p(Pol->Ray[i][0]) ||
       value_zero_p(Pol->Ray[i][nbColumns-1]))
      ++rays;
  
  /* Initialize the result */
  result=(Param_Polyhedron *)malloc(sizeof(Param_Polyhedron));
  result->nbV=nbRows-rays;
  result->V=NULL;
  
  /* Build the parametric vertices */
  for(i=0;i<nbRows;i++) {
    Matrix *vertex;
    Param_Vertices *paramVertex;
    int j;

    if (value_notone_p(Pol->Ray[i][0]) ||
        value_zero_p(Pol->Ray[i][nbColumns-1]))
      continue;

    vertex=Matrix_Alloc(nbColumns-2,2);
    for(j=1;j<nbColumns-1;j++) {
      value_assign(vertex->p[j-1][0],Pol->Ray[i][j]);
      value_assign(vertex->p[j-1][1],Pol->Ray[i][nbColumns-1]);
    }
    paramVertex=(Param_Vertices *)malloc(sizeof(Param_Vertices));
    paramVertex->Vertex=vertex;
    
    /* There is one validity domain : universe of dimension 0 */
    paramVertex->Domain=Matrix_Alloc(1,2);
    value_set_si(paramVertex->Domain->p[0][0],1);
    value_set_si(paramVertex->Domain->p[0][1],1);    
    paramVertex->next=result->V;
    result->V=paramVertex;
  }
  
  /* Build the parametric domains (only one here) */
  if (nbRows > 1)
    size=(nbRows-1)/(8*sizeof(int))+1;
  else
    size = 1;
  result->D=(Param_Domain *)malloc(sizeof(Param_Domain));
  result->D->next=NULL;
  result->D->Domain=Universe_Polyhedron(0);
  result->D->F=(unsigned int *)malloc(size*sizeof(int));
  memset(&result->D->F[0],0xFF,size*sizeof(int));
  
  return result;
} /* GenParamPolyhedron */


/*----------------------------------------------------------------------*/
/* PreElim_Columns                                                      */
/* function being called before Elim_Columns                            */
/* Equalities in E are analysed to initialize ref and p.                */
/* These two vectors are used to construct the new constraint matrix    */
/* PreElim_Columns returns the transformation matrix to re-convert the  */
/* resulting domains in the same format (by adding empty columns)       */
/* in the parameter space                                               */
/*----------------------------------------------------------------------*/
Matrix *PreElim_Columns(Polyhedron *E,int *p,int *ref,int m) {
	
  int i,j,l;
  Matrix *T;
  
  /* find which columns to eliminate */
  /* p contains, for each line in E, the column to eliminate */
  /* (i.e. the corresponding parameter index to eliminate) */
  /* 0 <= i <= E->NbEq, and  1 <= p[i] <= E->Dimension */
  memset(p,0,sizeof(int)*(E->NbEq));

#ifdef DEBUGPP32
  fprintf(stderr,"\n\nPreElim_Columns starting\n");
  fprintf(stderr,"E =\n");
  Polyhedron_Print(stderr,P_VALUE_FMT,E);
#endif
  
  for(l=0;l<E->NbEq;++l) {
    if(value_notzero_p(E->Constraint[l][0])) {
      fprintf(stderr,"Internal error: Elim_Columns (polyparam.c), expecting equalities first in E.\n");
      free(p);
      return(NULL);
    }
    for(i=1;value_zero_p(E->Constraint[l][i]);++i) {
      if(i==E->Dimension+1) {
	fprintf(stderr,"Internal error: Elim_Columns (polyparam.c), expecting non-empty constraint in E.\n");
	free(p);
	return( NULL );
      }
    }
    p[l] = i;
    
#ifdef DEBUGPP32
    fprintf(stderr,"p[%d] = %d,",l,p[l]);
#endif
  }

  /* Reference vector: column ref[i] in A corresponds to column i in M */
  for(i=0;i<E->Dimension+2-E->NbEq;++i) {
    ref[i]=i;
    for(j=0;j<E->NbEq;++j)
      if(p[j]<=ref[i])
	ref[i]++;
    
#ifdef DEBUGPP32
    fprintf(stderr,"ref[%d] = %d,",i,ref[i]);
#endif
  }
  
  /* Size of T : phdim-nbEq * phdim */
  T = Matrix_Alloc(m+1-E->NbEq,m+1);
  for(i=0;i<T->NbColumns;i++)
    for(j=0;j<T->NbRows;j++)
      if(ref[E->Dimension-m+j+1] == E->Dimension-m+i+1)
	value_set_si(T->p[j][i],1);
      else
	value_set_si(T->p[j][i],0);
  
#ifdef DEBUGPP32
  fprintf(stderr,"\nTransMatrix =\n");
  Matrix_Print(stderr,P_VALUE_FMT,T);
#endif
  
  return(T);

} /* PreElim_Columns */

/*----------------------------------------------------------------------*/
/* Elim_Columns                                                         */
/* Eliminate columns from A, using the equalities in E.                 */
/* ref and p are precalculated by PreElim_Columns, using E;             */
/* these two vectors are used to construct the new constraint matrix    */
/*----------------------------------------------------------------------*/
Polyhedron *Elim_Columns(Polyhedron *A,Polyhedron *E,int *p,int *ref) {
  
  int i,l,c;
  Matrix *M, *C;
  Polyhedron *R;
  Value tmp1,tmp2;
  
  /* Initialize all the 'Value' variables */
  value_init(tmp1); value_init(tmp2);
  
#ifdef DEBUGPP32
  fprintf(stderr,"\nElim_Columns starting\n");
  fprintf(stderr,"A =\n" );
  Polyhedron_Print(stderr,P_VALUE_FMT,A);
#endif
  
  /* Builds M = constraint matrix of A, useless columns zeroed */
  M = Polyhedron2Constraints(A);
  for(l=0;l<E->NbEq;++l) {    
    for(c=0;c<M->NbRows;++c) {
      if(value_notzero_p(M->p[c][p[l]])) {
	
	/* A parameter to eliminate here ! */
	for(i=1;i<M->NbColumns;++i) {
	  value_multiply(tmp1,M->p[c][i],E->Constraint[l][p[l]]);
	  value_multiply(tmp2,E->Constraint[l][i],M->p[c][p[l]]);
	  value_subtract(M->p[c][i],tmp1,tmp2);
	}
      }
    }
  } 
  
#ifdef DEBUGPP32
  fprintf(stderr,"\nElim_Columns after zeroing columns of A.\n");
  fprintf(stderr,"M =\n");
  Matrix_Print(stderr,P_VALUE_FMT,M);
#endif
  
  /* Builds C = constraint matrix, useless columns eliminated */
  C = Matrix_Alloc(M->NbRows,M->NbColumns - E->NbEq);
  for(l=0;l<C->NbRows;++l)
    for(c=0;c<C->NbColumns;++c) {
      value_assign(C->p[l][c],M->p[l][ref[c]]);
    }
    
#ifdef DEBUGPP32
  fprintf(stderr,"\nElim_Columns after eliminating columns of A.\n");
  fprintf(stderr,"C =\n");
  Matrix_Print(stderr,P_VALUE_FMT,C);
#endif
    
  R = Constraints2Polyhedron(C,ws);
  Matrix_Free(C);
  Matrix_Free(M);
  value_clear(tmp1); value_clear(tmp2);
  return(R);
} /* Elim_Columns */

/* 
 * Given a polyhedron 'Di' in combined data and parameter space and a context 
 * polyhedron 'C' representing the constraints on the parameter space, create
 * a list of parameterized vertices and assign values to global variables: 
 * m,n,ws,Sat,egalite,mf,Xi,Pi,PiInv,RaysDi,CEqualities. 
 */
Param_Polyhedron *Find_m_faces(Polyhedron **Di,Polyhedron *C,int keep_dom,int working_space,Polyhedron **CEq,Matrix **CT) {	
  
  unsigned int *mf;
  int i, j;
  Polyhedron *D=*Di,
             *C1,         /* true context */
             *D1;         /* the combined polyhedron, including context C */
  Matrix *M;
  Param_Polyhedron *res;
  int *p, *ref;

  if(CT) {
    *CEq = NULL;
    *CT = NULL;
  }
  
  if(!D || !C) 
    return (Param_Polyhedron *) 0;

  ws = working_space;
  m = C->Dimension;
  n = D->Dimension - m;
  if(n<0) {
    fprintf(stderr,
	    "Find_m_faces: ?%d parameters of a %d-polyhedron !\n",m,n);
    return (Param_Polyhedron *) 0;
  }
  if (m==0)
    return GenParamPolyhedron(D);
  
  /* Add constraints from Context to D -> result in D1 */
  C1 = align_context(C,D->Dimension,ws);

#ifdef DEBUGPP31
  fprintf(stderr,"m = %d\n",m);
  fprintf(stderr, "D = ");
  Polyhedron_Print(stderr,P_VALUE_FMT,D);
  fprintf(stderr,"C1 = ");
  Polyhedron_Print(stderr,P_VALUE_FMT,C1);
#endif
  
  D1 = DomainIntersection(D,C1,ws);

#ifdef DEBUGPP31
  fprintf(stderr,"D1 = ");
  Polyhedron_Print(stderr,P_VALUE_FMT,D1);
#endif
  
  Domain_Free(C1);

  if(!D1 || emptyQ(D1))
    return(NULL);
  
  /* Compute the true context C1 */
  /* M : lines in the direction of the first n indices (index space) */
  M   = Matrix_Alloc(n, D1->Dimension+2);
  for (i=0; i<n; i++)
    value_set_si(M->p[i][i+1],1);
  C1 = DomainAddRays(D1,M,ws);
  Matrix_Free(M);
  
#ifdef DEBUGPP31
  fprintf(stderr,"True context C1 = ");
  Polyhedron_Print(stderr,P_VALUE_FMT,C1);
#endif
  
  /* CEqualities contains the constraints (to be added again later) */
  /* *CT is the transformation matrix to add the removed parameters */
  if(!CT) {
    if (C1->NbEq == 0) {
      Polyhedron_Free(C1);
      CEqualities = NULL;
    } else {
      Polyhedron *CEq1,	/* CEqualities, in homogeneous dim */
	         *C2,	/* C1 (temporary) simplified */
	         *D2;	/* D1, (temporary) simplified */
      
      /* Remove equalities from true context C1 and from D1             */     
      /* Compute CEqualities = matrix of equalities in C1, projected in */
      /* the parameter space                                            */
      M = Matrix_Alloc(C1->NbEq,m+2);
      for(j=0,i=0;i<C1->NbEq;++i,++j) {
	while(value_notzero_p(C1->Constraint[j][0]))
	  ++j;
	value_assign(M->p[i][0],C1->Constraint[j][0]);
	Vector_Copy(&C1->Constraint[j][D->Dimension-m+1],&M->p[i][1],(m+1));
      }
      CEqualities = Constraints2Polyhedron(M,ws);
      Matrix_Free(M);
      CEq1 = align_context(CEqualities,D->Dimension,ws);

      /* Simplify D1 and C1 (remove the equalities) */
      D2 = DomainSimplify(D1,CEq1,ws);
      Polyhedron_Free(D1);
      Polyhedron_Free(C1);
      Polyhedron_Free(CEq1);
      D1 = D2;
      C1 = NULL;
    }
  }
  else { /* if( CT  ) */
    Polyhedron *CEq1,	/* CEqualities */
               *C2,	/* C1 (temporary) simplified */
               *D2;	/* D1, (temporary) simplified */

    /* Suppress all useless constraints in parameter domain */
    /* when CT is not NULL (ehrhart) */
    /* Vin100, march 01 */
    CEq1 = C1;
    M = Matrix_Alloc(C1->NbConstraints,m+2);
    for(i=0;i<C1->NbConstraints;++i) {
      value_assign(M->p[i][0],C1->Constraint[i][0]);
      Vector_Copy(&C1->Constraint[i][D->Dimension-m+1],&M->p[i][1],(m+1));
    }
    CEqualities = Constraints2Polyhedron( M, ws );
    Matrix_Free(M);

    D2 = DomainSimplify(D1,CEq1,ws);
    Polyhedron_Free(D1);
    D1 = D2;
    C1 = Universe_Polyhedron(D2->Dimension);
    
    /* if CT is not NULL, the constraints are eliminated                */
    /* *CT will contain the transformation matrix to come back to the   */
    /* original dimension (for a polyhedron, in the parameter space)    */
    if( CEq1->NbEq )
    {
      m -= CEq1->NbEq;
      p = (int *)malloc(sizeof(int)*(CEq1->NbEq));
    }
    else
      p = NULL;
    ref = (int*) malloc(sizeof(int)*
			(CEq1->Dimension+2-CEq1->NbEq));
    *CT = PreElim_Columns(CEq1,p,ref,CEqualities->Dimension);
    D2 = Elim_Columns(D1,CEq1,p,ref);
    if (p)
      free(p);
    free(ref);
    
#ifdef DEBUGPP3
    fprintf(stderr,"D2\t Dim = %3d\tNbEq = %3d\tLines = %3d\n",
	    D2->Dimension,D2->NbEq,D2->NbBid);
    C2 = Elim_Columns(C1,CEq1,p,ref);
    fprintf(stderr,"C2\t Dim = %3d\tNbEq = %3d\tLines = %3d\n",
	    C2->Dimension,C2->NbEq,C2->NbBid);
    Polyhedron_Free(C2);
#endif
    
    Polyhedron_Free(D1);
    Polyhedron_Free(C1);
    D1 = D2;
    C1 = NULL;
    *CEq = CEqualities;
    
#ifdef DEBUGPP3
    fprintf(stderr,"Polyhedron CEq = ");
    Polyhedron_Print(stderr,P_VALUE_FMT,*CEq);
    fprintf(stderr,"Matrix CT = ");
    Matrix_Print(stderr,P_VALUE_FMT,*CT);
#endif
    
    Polyhedron_Free(CEq1);
    CEqualities = NULL;	/* don't simplify ! */

    /* m changed !!! */
    if(m==0) {
      /* return the new D1 too */
      *Di = D1;

      return GenParamPolyhedron(D1);
    }
  }

#ifdef DEBUGPP3
  fprintf(stderr,"Polyhedron D1 (D AND C) = ");
  Polyhedron_Print(stderr,P_VALUE_FMT, D1);
  fprintf(stderr,"Polyhedron CEqualities = ");
  if(CEqualities) Polyhedron_Print(stderr,P_VALUE_FMT, CEqualities);
  else fprintf(stderr,"NULL\n");
#endif
  
  KD = keep_dom;
  PDomains = NULL;
  PV_Result = NULL;
  nbPV = 0;
  
  if (D1->NbRays==0) return 0;
  
  /* mf : a bit array indicating which rays are part of the m-face */
  /* Poly2Sat initializes mf to all ones */
  /* set global variable nr to size (number of words) of mf */
  Sat = Poly2Sat(D1,&mf);

#ifdef DEBUGPP4
    fprintf(stderr,"Sat = ");
    SMPrint(Sat);

  fprintf(stderr,"mf = ");
  for (i=0; i<nr; i++)
    fprintf(stderr,"%08X", mf[i]);
  fprintf(stderr, "\n");
#endif
  
  /* A boolean array saying which constraints are part of the m-face */
  egalite = (int *)malloc(sizeof(int)*D1->NbConstraints );
  memset(egalite,0, sizeof(int)*D1->NbConstraints);

  for (i=0; i<D1->NbEq; i++)
    egalite[i] = 1;

  Xi     = Matrix_Alloc(n+1,m+1);
  Pi     = Matrix_Alloc(m+1,m+1);
  PiTest = Matrix_Alloc(m+1,m+1);
  CTest  = Matrix_Alloc(D->Dimension+1,D->NbConstraints);
  PiInv  = Matrix_Alloc(m+1,m+2);
  RaysDi = Matrix_Alloc(D1->NbRays,m+2);
  m_dim = m;
  
#ifdef DEBUGPP
  nbfaces=0;
#endif
#ifdef DEBUGPP3
  fprintf(stderr,
	  "Target: find faces that saturate %d constraints and %d rays/lines\n",
	  D1->Dimension - m_dim,m_dim+1);
#endif
	
  /* D1->NbEq constraints already saturated ! */
  scan_m_face(D1->NbEq,(D1->Dimension - m_dim - D1->NbEq),D1,mf);

  /* pos, number of constraints needed     */

#ifdef DEBUGPP
  fprintf( stderr, "Number of m-faces: %d\n", nbfaces );
#endif
  
  Matrix_Free(RaysDi);
  Matrix_Free(PiInv);
  Matrix_Free(PiTest);
  Matrix_Free(CTest);
  Matrix_Free(Pi);
  Matrix_Free(Xi);
  free(egalite);
  free(mf);
  SMFree(Sat);
  
  /*	if(CEqualities && keep_dom==0) {
	   Domain_Free(CEqualities);
	}
  */

  if(CT)		/* return the new D1 too ! */
    *Di = D1;
  else
    Domain_Free(D1);

  res = (Param_Polyhedron *) malloc (sizeof(Param_Polyhedron));
  res->nbV = nbPV;
  res->V = PV_Result;
  res->D = PDomains;
  return(res);
} /* Find_m_faces */

/*
 * Given parametric domain 'PD' and number of parametric vertices 'nb_domains',
 * find the vertices that belong to distinct sub-domains. 
 */
void Compute_PDomains(Param_Domain *PD,int nb_domains,int working_space) {
  
  unsigned bx;
  int i, ix, nv;
  Polyhedron *dx, *d1, *d2;
  Param_Domain *p1, *p2, *p2prev, *PDNew;
  
  if (nb_domains==0) {
    
#ifdef DEBUGPP5
    fprintf(stderr,"No domains\n");
#endif
    
    return;
  }

  /* Already filled out by GenParamPolyhedron */
  if (!PD->next && PD->F)
    return;

   /* Initialization */
   nv = (nb_domains - 1)/(8*sizeof(int)) + 1;

#ifdef DEBUGPP5
   fprintf(stderr,"nv = %d\n",nv);
#endif

   for(p1=PD,i=0,ix=0,bx=MSB;p1;p1=p1->next,i++) {
     
     /* Assign a bit array 'p1->F' of suitable size to include the vertices */
     p1->F = (unsigned *) malloc (nv * sizeof(unsigned));
     
     /* Set the bit array to zeros */
     memset(p1->F,0,nv * sizeof(unsigned));
     p1->F[ix] |= bx;      /* Set i'th bit to one */
     NEXT(ix, bx);
   }

#ifdef DEBUGPP5
   fprintf(stderr,"nb of vertices=%d\n",i);
#endif

   /* Walk the PD list with two pointers */
   ix = 0; bx=MSB;
   for (p1=PD;p1;p1=p1->next) {
     for (p2prev=p1,p2=p1->next;p2;p2prev=p2,p2=p2->next) {
	
       /* Find intersection */
       dx = PDomainIntersection(p1->Domain,p2->Domain,working_space);
       
       if (!dx || emptyQ(dx)) {
	 
#ifdef DEBUGPP5
	 fprintf( stderr, "Empty dx (p1 inter p2). Continuing\n");
#endif
	 if(dx)
	   Domain_Free(dx);
	 continue;
       }

#ifdef DEBUGPP5      
       fprintf(stderr,"Begin PDomainDifference\n");
       fprintf(stderr, "p1=");
       Polyhedron_Print(stderr,P_VALUE_FMT,p1->Domain);
       fprintf(stderr,"p2=");
       Polyhedron_Print(stderr,P_VALUE_FMT,p2->Domain);
#endif
       d1 = PDomainDifference(p1->Domain,p2->Domain,working_space);
       d2 = PDomainDifference(p2->Domain,p1->Domain,working_space);       

#ifdef DEBUGPP5       
       fprintf(stderr,"p1\\p2=");
       Polyhedron_Print(stderr,P_VALUE_FMT,d1);
       fprintf(stderr,"p2\\p1=");
       Polyhedron_Print(stderr,P_VALUE_FMT,d2);
       fprintf(stderr,"END PDomainDifference\n\n");       
#endif
       if (!d1 || emptyQ(d1) || d1->NbEq!=0) {

#ifdef DEBUGPP5
	 fprintf(stderr,"Empty d1\n");
#endif
	 if (d1) 
	   Domain_Free(d1);
	 Domain_Free(dx);
	 
	 if (!d2 || emptyQ(d2) || d2->NbEq!=0) {
	   
#ifdef DEBUGPP5
	   fprintf( stderr, "Empty d2 (deleting)\n");
#endif
	   /* dx = p1->Domain = p2->Domain */
	   if (d2) Domain_Free(d2);
	   
	   /* Update p1 */
	   for (i=0;i<nv;i++)
	     p1->F[i] |= p2->F[i];
	   
	   /* Delete p2 */
	   p2prev->next = p2->next;
	   Domain_Free(p2->Domain);
	   free(p2->F);
	   free(p2);
	   p2 = p2prev;
	 }
	 else {  /* d2 is not empty --> dx==p1->domain */
	   
#ifdef DEBUGPP5
	   fprintf( stderr, "p2 replaced by d2\n");
#endif
	   /* Update p1 */
	   for(i=0;i<nv;i++)
	     p1->F[i] |= p2->F[i];
	   
	   /* Replace p2 with d2 */
	   Domain_Free( p2->Domain );
	   p2->Domain = d2;
	 }
       }
       else { /* d1 is not empty */         
	 if (!d2 || emptyQ(d2) || d2->NbEq!=0) {
	   
#ifdef DEBUGPP5
	   fprintf( stderr, "p1 replaced by d1\n");
#endif
	   if (d2) Domain_Free(d2);
	   
	   /* dx = p2->domain */
	   Domain_Free(dx);
	   
	   /* Update p2 */
	   for(i=0;i<nv;i++)
	     p2->F[i] |= p1->F[i];

	   /* Replace p1 with d1 */
	   Domain_Free(p1->Domain);
	   p1->Domain = d1;
	 }
	 else { /*d2 is not empty-->d1,d2,dx are distinct */
	   
#ifdef DEBUGPP5
	   fprintf(stderr,"Non-empty d1 and d2\nNew node created\n");
#endif
	   /* Create a new node for dx */
	   PDNew = (Param_Domain *) malloc( sizeof(Param_Domain) );
	   PDNew->F = (unsigned int *)malloc( nv*sizeof(int) );
	   memset(PDNew->F,0,nv*sizeof(int));
	   PDNew->Domain = dx;
	   
	   for (i=0;i<nv;i++)
	     PDNew->F[i] = p1->F[i] | p2->F[i];
	   
	   /* Replace p1 with d1 */
	   Domain_Free( p1->Domain );
	   p1->Domain = d1;
	   
	   /* Replace p2 with d2 */
	   Domain_Free( p2->Domain );
	   p2->Domain = d2;
	   
	   /* Insert new node after p1 */
	   PDNew->next = p1->next;
	   p1->next = PDNew;
	 }
       }
     }  /* end of p2 scan */
     if (p1->Domain->next) {
	Polyhedron *C = DomainConvex(p1->Domain, working_space);
	Domain_Free(p1->Domain);
	p1->Domain = C;
     }
   } /* end of p1 scan */
} /* Compute_PDomains */
					
/* 
 * Given a polyhedron 'Din' in combined data and parametre space, a context
 * polyhedron 'Cin' representing the constraints on the parameter space and 
 * a working space size 'working_space', return a parametric polyhedron with
 * a list of parametric vertices and their defining domains. 
 */
Param_Polyhedron *Polyhedron2Param_Vertices(Polyhedron *Din,Polyhedron *Cin,int working_space) {
  
  Param_Polyhedron *result;
  
  POL_ENSURE_FACETS(Din);
  POL_ENSURE_VERTICES(Din);
  POL_ENSURE_FACETS(Cin);
  POL_ENSURE_VERTICES(Cin);
 
#ifdef DEBUGPP
  fprintf(stderr,"Polyhedron2Param_Vertices algorithm starting at : %.2fs\n",
	  (float)clock()/CLOCKS_PER_SEC);
#endif
  
  /***************** Scan the m-faces ****************/
  result = Find_m_faces(&Din,Cin,0,working_space,NULL,NULL);
  
#ifdef DEBUGPP
  fprintf(stderr, "nb of points : %d\n",result->nbV);
#endif
  
#ifdef DEBUGPP
  fprintf(stderr, "end main loop : %.2fs\n", (float)clock()/CLOCKS_PER_SEC);
#endif
  
  return(result);
} /* Polyhedron2Param_Vertices */

/*
 * Free the memory allocated to a list of parametrized vertices  
 */
void Param_Vertices_Free(Param_Vertices *PV) {
  
  Param_Vertices *next_pv;
  
  while(PV) {
    next_pv = PV->next;
    if (PV->Vertex) Matrix_Free(PV->Vertex);
    if (PV->Domain) Matrix_Free(PV->Domain);
    free(PV);  
    PV = next_pv;
  }
} /* Param_Vertices_Free */

/*
 * Print a list of parametrized vertices *
 */
void Print_Vertex(FILE *DST,Matrix *V,char **param_names){
  
  int l, v;
  int first;
  Value gcd,tmp;
  
  value_init(gcd); value_init(tmp);
  
  fprintf(DST, "[" );
  for(l=0;l<V->NbRows;++l){
    
    /* Variables */
    first=1;
    fprintf(DST, " " );
    for(v=0;v < V->NbColumns-2;++v) {
      if(value_notzero_p(V->p[l][v])) {
	Gcd(V->p[l][v],V->p[l][V->NbColumns-1],&gcd);
	value_absolute(gcd,gcd);
	value_division(tmp,V->p[l][v],gcd);
	if(value_posz_p(tmp)) {
	  if(!first) 
	    fprintf(DST, "+");
	  value_division(tmp,V->p[l][v],gcd);
	  if(value_notone_p(tmp)) { 
	    value_print(DST,VALUE_FMT,tmp);
	  }  
	}
	else { /* V->p[l][v]/gcd<0 */
	  value_division(tmp,V->p[l][v],gcd);
	  if(value_mone_p(tmp))
	    fprintf(DST, "-" );
	  else {
	    value_print(DST,VALUE_FMT,tmp);
	  }
	}
	value_division(tmp,V->p[l][V->NbColumns-1],gcd);
	if(value_notone_p(tmp)) {
	  fprintf(DST, "%s/", param_names[v]);
	  value_print(DST,VALUE_FMT,tmp);
	}
	else
	  fprintf(DST, "%s", param_names[v]);
	first=0;
      }
    }

    /* Constant */
    if(value_notzero_p(V->p[l][v]) || first) {
      if(value_posz_p(V->p[l][v]) && !first)
	fprintf(DST,"+");
	Gcd(V->p[l][v],V->p[l][V->NbColumns-1],&gcd);
      value_absolute(gcd,gcd);
      value_division(tmp,V->p[l][v],gcd);
      value_print(DST,VALUE_FMT,tmp);
      value_division(tmp,V->p[l][V->NbColumns-1],gcd);
      if(value_notone_p(tmp)) {
	fprintf(DST,"/");
	value_print(DST,VALUE_FMT,tmp);
	fprintf(DST," ");
      }
    }
    if (l<V->NbRows-1) 
      fprintf(DST, ", ");
  }
  fprintf(DST, " ]");
  value_clear(gcd); value_clear(tmp);
  return;
} /* Print_Vertex */

/*----------------------------------------------------------------------*/
/* VertexCT                                                             */
/* convert a paramvertex from reduced space to normal m-space           */
/*----------------------------------------------------------------------*/
Matrix *VertexCT(Matrix *V,Matrix *CT) {
  
  Matrix *Vt;
  int i,j,k;
  
  if(CT) {
    
    /* Have to transform the vertices to original dimension */
    Vt = Matrix_Alloc(V->NbRows,CT->NbColumns+1);
    for(i=0;i<V->NbRows;++i) {
      value_assign(Vt->p[i][CT->NbColumns],V->p[i][V->NbColumns-1]);
      for(j=0;j<CT->NbColumns;j++) {
	for(k=0;k<CT->NbRows;k++)
	  if(value_notzero_p(CT->p[k][j]))
	    break;
	if(k<CT->NbRows)
	  value_assign(Vt->p[i][j],V->p[i][k]);
	else
	  value_set_si(Vt->p[i][j],0);
      }
    }
    return(Vt);
  }
  else
    return(NULL);
} /* VertexCT */

/*
 * Print the validity Domain 'D' of a parametric polyhedron 
 */ 
void Print_Domain(FILE *DST,Polyhedron *D,char **pname) {
  
  int l, v;
  int first;
  
  POL_ENSURE_FACETS(D);
  POL_ENSURE_VERTICES(D);
 
  for(l=0;l<D->NbConstraints;++l) {
    fprintf(DST, "         ");
    first = 1;
    for(v=1;v<=D->Dimension;++v) {
      if(value_notzero_p(D->Constraint[l][v])) {
	if(value_one_p(D->Constraint[l][v])) {
	  if(first)
	    fprintf(DST, "%s ", pname[v-1]);
	  else
	    fprintf(DST, "+ %s ", pname[v-1] );
	}
	else if(value_mone_p(D->Constraint[l][v]))
	  fprintf(DST, "- %s ", pname[v-1] );
	else {
	  if(value_pos_p(D->Constraint[l][v]) && !first )
	    fprintf(DST, "+ " );
	  value_print(DST,VALUE_FMT,D->Constraint[l][v]);
	  fprintf(DST,"%s ",pname[v-1]);
	}
	first = 0;
      }
    }
    if(value_notzero_p(D->Constraint[l][v])) {
      if(value_pos_p(D->Constraint[l][v]) && !first)
	fprintf(DST,"+");
      fprintf(DST," ");
      value_print(DST,VALUE_FMT,D->Constraint[l][v]);
    }
    fprintf(DST,(value_notzero_p(D->Constraint[l][0])) ?" >= 0":" = 0");
    fprintf(DST, "\n" );
  }
  fprintf(DST, "\n");

	if( D->next )
	{
		fprintf( DST, "UNION\n" );
		Print_Domain( DST, D->next, pname );
	}
  return;
} /* Print_Domain */

/* 
 * Given a list of parametrized vertices and an array of parameter names, Print
 * a list of parametrized vertices in a comprehensible format. 
 */
void Param_Vertices_Print(FILE *DST,Param_Vertices *PV,char **param_names) {
   
  Polyhedron *poly;
  
  while(PV) {
    fprintf(DST, "Vertex :\n" );
    Print_Vertex(DST,PV->Vertex,param_names);
    
    /* Pour le domaine : */
    fprintf(DST, "   If :\n" );
    poly = Constraints2Polyhedron(PV->Domain,200);
    Print_Domain(DST,poly,param_names);
    Domain_Free(poly);   
    PV = PV->next;
  }
  return;
} /* Param_Vertices_Print */

/* 
 * Given a polyhedron 'Din' in combined data and parametre space, a context
 * polyhedron 'Cin' representing the constraints on the parameter space and 
 * a working space size 'working_space', return a parametric polyhedron with
 * a list of distinct validity domains and a complete list of valid vertices 
 * associated to each validity domain. 
 */
Param_Polyhedron *Polyhedron2Param_Domain(Polyhedron *Din,Polyhedron *Cin,int working_space) {
		
  Param_Polyhedron *result;
  Param_Domain *D;

  POL_ENSURE_FACETS(Din);
  POL_ENSURE_VERTICES(Din);
  POL_ENSURE_FACETS(Cin);
  POL_ENSURE_VERTICES(Cin);
 
#ifdef DEBUGPP
  fprintf(stderr,"Polyhedron2Param_Polyhedron algorithm starting at : %.2fs\n",
	  (float)clock()/CLOCKS_PER_SEC);
#endif
  
  /* Find the m-faces, keeping the corresponding domains */
  /* in the linked list PDomains */
  result = Find_m_faces(&Din,Cin,1,working_space,NULL,NULL);
  
#ifdef DEBUGPP
  if(result) fprintf(stderr, "Number of vertices : %d\n",result->nbV);
  fprintf(stderr,"Vertices found at : %.2fs\n",(float)clock()/CLOCKS_PER_SEC);
#endif
  
  /* Processing of PVResult and PDomains */
  if(result && Cin->Dimension>0)		/* at least 1 parameter */
    Compute_PDomains(result->D,result->nbV,working_space);
  if(result && CEqualities)
    for(D=result->D;D;D=D->next)
      D->Domain = Add_CEqualities(D->Domain);
  Polyhedron_Free(CEqualities); 
  
#ifdef DEBUGPP
  fprintf(stderr, "domains found at : %.2fs\n", (float)clock()/CLOCKS_PER_SEC);
#endif

  return(result);
} /* Polyhedon2Param_Domain */

/*
 *
 */
Param_Polyhedron *Polyhedron2Param_SimplifiedDomain(Polyhedron **Din,Polyhedron *Cin,int working_space,Polyhedron **CEq,Matrix **CT) {
						     
  Param_Polyhedron *result;

  assert(CEq != NULL);
  assert(CT != NULL);

  POL_ENSURE_FACETS(*Din);
  POL_ENSURE_VERTICES(*Din);
  POL_ENSURE_FACETS(Cin);
  POL_ENSURE_VERTICES(Cin);
 
#ifdef DEBUGPP
  fprintf(stderr,"Polyhedron2Param_Polyhedron algorithm starting at : %.2fs\n",
	  (float)clock()/CLOCKS_PER_SEC);
#endif
  
  /* Find the m-faces, keeping the corresponding domains */
  /* in the linked list PDomains */
  result = Find_m_faces(Din,Cin,1,working_space,CEq,CT);
  
#ifdef DEBUGPP
  if(result) fprintf(stderr, "Number of vertices : %d\n",result->nbV);
  fprintf(stderr,"Vertices found at : %.2fs\n",(float)clock()/CLOCKS_PER_SEC);
#endif

  /* Processing of PVResult and PDomains */
  if(result && Cin->Dimension>0)		/* at least 1 parameter */
    Compute_PDomains(result->D,result->nbV,working_space);
  
  /* Removed this, Vin100, March 01 */
  /*	if(result && CEqualities )
	for(D=result->D;D;D=D->next)
	D->Domain = Add_CEqualities(D->Domain);
  */
  
#ifdef DEBUGPP
  fprintf(stderr, "domains found at : %.2fs\n", (float)clock()/CLOCKS_PER_SEC);
#endif
  
  return(result);
} /* Polyhedron2Param_SimplifiedDomain */

/*
 * Free the memory allocated to a list of validity domain of a parametrized 
 * polyhedron.
 */
void Param_Domain_Free(Param_Domain *PD) {
	
  Param_Domain *next_pd;
  
  while(PD) {
    free(PD->F);
    Domain_Free(PD->Domain);
    next_pd = PD->next;
    free(PD);
    PD = next_pd;
  }
  return;
} /* Param_Domain_Free */

/*
 * Free the memory allocated to a parametric polyhedron 'P' 
 */
void Param_Polyhedron_Free(Param_Polyhedron *P) {
  
  if (!P) return;
  Param_Vertices_Free(P->V);
  Param_Domain_Free(P->D);
  free(P);
  return;
} /* Param_Polyhedron_Free */




