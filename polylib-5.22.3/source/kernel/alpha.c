/* polytest.c */
#include <stdio.h>
#include <polylib/polylib.h>

/* alpha.c
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*---------------------------------------------------------------------*/
/* int exist_points(pos,P,context)                                     */
/*    pos : index position of current loop index (0..hdim-1)           */
/*    P: loop domain                                                   */
/*    context : context values for fixed indices                       */
/* recursive procedure, recurs for each imbriquation                   */
/* returns 1 if there exists any integer points in this polyhedron     */
/* returns 0 if there are none                                         */
/*---------------------------------------------------------------------*/
static int exist_points(int pos,Polyhedron *Pol,Value *context) {
  
  Value LB, UB, k,tmp;
  
  value_init(LB); value_init(UB); 
  value_init(k);  value_init(tmp);
  value_set_si(LB,0);
  value_set_si(UB,0);
  
  /* Problem if UB or LB is INFINITY */
  if (lower_upper_bounds(pos,Pol,context,&LB,&UB) !=0) {
    errormsg1("exist_points", "infdom", "infinite domain");
    value_clear(LB);
    value_clear(UB);
    value_clear(k);
    value_clear(tmp);
    return -1;
  }
  value_set_si(context[pos],0);
  if(value_lt(UB,LB)) {
    value_clear(LB); 
    value_clear(UB);
    value_clear(k);
    value_clear(tmp);
    return 0;
  }  
  if (!Pol->next) {
    value_subtract(tmp,UB,LB);
    value_increment(tmp,tmp);
    value_clear(UB);
    value_clear(LB);
    value_clear(k);
    return (value_pos_p(tmp));
  }
  
  for (value_assign(k,LB);value_le(k,UB);value_increment(k,k)) {
    
    /* insert k in context */
    value_assign(context[pos],k);    
    if (exist_points(pos+1,Pol->next,context) > 0 ) {
      value_clear(LB); value_clear(UB);
      value_clear(k); value_clear(tmp);
      return 1;
    }
  }   
  /* Reset context */
  value_set_si(context[pos],0);
  value_clear(UB); value_clear(LB);
  value_clear(k); value_clear(tmp);
  return 0;
}
    
/*--------------------------------------------------------------*/
/* Test to see if there are any integral points in a polyhedron */
/*--------------------------------------------------------------*/
int Polyhedron_Not_Empty(Polyhedron *P,Polyhedron *C,int MAXRAYS) {

  int res,i;
  Value *context;
  Polyhedron *L;
  
  POL_ENSURE_FACETS(P);
  POL_ENSURE_VERTICES(P);
  POL_ENSURE_FACETS(C);
  POL_ENSURE_VERTICES(C);

  /* Create a context vector size dim+2 and set it to all zeros */
  context = (Value *) malloc((P->Dimension+2)*sizeof(Value));
  
  /* Initialize array 'context' */
  for (i=0;i<(P->Dimension+2);i++) 
    value_init(context[i]);
  
  Vector_Set(context,0,(P->Dimension+2));
  
  /* Set context[P->Dimension+1] = 1  (the constant) */
  value_set_si(context[P->Dimension+1],1);
  
  L = Polyhedron_Scan(P,C,MAXRAYS);
  res = exist_points(1,L,context);
  Domain_Free(L);
  
  /* Clear array 'context' */
  for (i=0;i<(P->Dimension+2);i++) 
    value_clear(context[i]);
  free(context);
  return res;
}

/* PolyhedronLTQ ( P1, P2 ) */
/* P1 and P2 must be simple polyhedra */
/* result =  0 --> not comparable */
/* result = -1 --> P1 < P2        */
/* result =  1 --> P1 > P2        */
/* INDEX  = 1 .... Dimension      */
int PolyhedronLTQ (Polyhedron *Pol1,Polyhedron *Pol2,int INDEX, int PDIM, int NbMaxConstrs) { 
  
  int res, dim, i, j, k;
  Polyhedron *Q1, *Q2, *Q3, *Q4, *Q;
  Matrix *Mat;

  if (Pol1->next || Pol2->next) {
    errormsg1("PolyhedronLTQ", "compoly", "Can only compare polyhedra");
    return 0;
  }
  if (Pol1->Dimension != Pol2->Dimension) {
    errormsg1("PolyhedronLTQ","diffdim","Polyhedra are not same dimension");
    return 0;
  }
  dim = Pol1->Dimension+2;

  POL_ENSURE_FACETS(Pol1);
  POL_ENSURE_VERTICES(Pol1);
  POL_ENSURE_FACETS(Pol2);
  POL_ENSURE_VERTICES(Pol2);
  
#ifdef DEBUG
  fprintf(stdout, "P1\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Pol1);
  fprintf(stdout, "P2\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Pol2);
#endif
  
  /* Create the Line to add */
  k = Pol1->Dimension-INDEX+1-PDIM;
  Mat = Matrix_Alloc(k,dim);
  Vector_Set(Mat->p_Init,0,dim*k);
  for(j=0,i=INDEX;j<k;i++,j++)
    value_set_si(Mat->p[j][i],1);
  
  Q1 = AddRays(Mat->p[0],k,Pol1,NbMaxConstrs);
  Q2 = AddRays(Mat->p[0],k,Pol2,NbMaxConstrs);

#ifdef DEBUG
  fprintf(stdout, "Q1\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Q1);
  fprintf(stdout, "Q2\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Q2);
#endif
  
  Matrix_Free(Mat);
  Q  = DomainIntersection(Q1,Q2,NbMaxConstrs);
  
#ifdef DEBUG
  fprintf(stdout, "Q\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Q);
#endif
  
  Domain_Free(Q1);
  Domain_Free(Q2);
  
  if (emptyQ(Q)) res = 0;	/* not comparable */
  else {
    Q1 = DomainIntersection(Pol1,Q,NbMaxConstrs);
    Q2 = DomainIntersection(Pol2,Q,NbMaxConstrs);
    
#ifdef DEBUG
    fprintf(stdout, "Q1\n");
    Polyhedron_Print(stdout,P_VALUE_FMT,Q1);
    fprintf(stdout, "Q2\n");
    Polyhedron_Print(stdout,P_VALUE_FMT,Q2);
#endif

    k = Q1->NbConstraints + Q2->NbConstraints;
    Mat = Matrix_Alloc(k, dim);
    Vector_Set(Mat->p_Init,0,k*dim);
    
    /* First compute surrounding polyhedron */    
    j=0;
    for (i=0; i<Q1->NbConstraints; i++) {
      if ((value_one_p(Q1->Constraint[i][0])) && (value_pos_p(Q1->Constraint[i][INDEX]))) {
	
	/* keep Q1's lower bounds */
	for (k=0; k<dim; k++) 
	  value_assign(Mat->p[j][k],Q1->Constraint[i][k]);
	j++;
      }
    }
    for (i=0; i<Q2->NbConstraints; i++) {
      if ((value_one_p(Q2->Constraint[i][0])) && (value_neg_p(Q2->Constraint[i][INDEX]))) {
	
	/* and keep Q2's upper bounds */
	for (k=0; k<dim; k++) 
	  value_assign(Mat->p[j][k],Q2->Constraint[i][k]);
	j++;
      }
    }
    Q4 = AddConstraints(Mat->p[0], j, Q, NbMaxConstrs);
    Matrix_Free(Mat);
    
#ifdef debug
    fprintf(stderr, "Q4 surrounding polyhedron\n");
    Polyhderon_Print(stderr,P_VALUE_FMT, Q4);
#endif

    /* if surrounding polyhedron is empty, D1>D2 */
    if (emptyQ(Q4)) {
      res = 1;
      
#ifdef debug
      fprintf(stderr, "Surrounding polyhedron is empty\n");
#endif
      goto LTQdone2; 
    }
    
    /* Test if Q1 < Q2 */      
    /* Build a constraint array for >= Q1 and <= Q2 */
    Mat = Matrix_Alloc(2,dim);
    Vector_Set(Mat->p_Init,0,2*dim);
    
    /* Choose a contraint from Q1 */
    for (i=0; i<Q1->NbConstraints; i++) {
      if (value_zero_p(Q1->Constraint[i][0])) {
	
	/* Equality */
	if (value_zero_p(Q1->Constraint[i][INDEX])) {
	  
	  /* Ignore side constraint (they are in Q) */
	  continue;
	}
	else if (value_neg_p(Q1->Constraint[i][INDEX])) {
	  
	  /* copy -constraint to Mat */
	  value_set_si(Mat->p[0][0],1);
	  for (k=1; k<dim; k++)
	    value_oppose(Mat->p[0][k],Q1->Constraint[i][k]);
	}
	else {
	  
	  /* Copy constraint to Mat */
	  
	  value_set_si(Mat->p[0][0],1);
	  for (k=1; k<dim; k++)
	    value_assign(Mat->p[0][k],Q1->Constraint[i][k]);
	}
      }
      else if(value_neg_p(Q1->Constraint[i][INDEX])) {
	
	/* Upper bound -- make a lower bound from it */
	value_set_si(Mat->p[0][0],1);
	for (k=1; k<dim; k++)
	  value_oppose(Mat->p[0][k],Q1->Constraint[i][k]);
      }
      else {	
	
	/* Lower or side bound -- ignore it */
	continue;
      }
      
      /* Choose a constraint from Q2 */
      for (j=0; j<Q2->NbConstraints; j++) {
	if (value_zero_p(Q2->Constraint[j][0])) {   /* equality */
	  if (value_zero_p(Q2->Constraint[j][INDEX])) {
	    
	    /* Ignore side constraint (they are in Q) */
	    continue;
	  }
	  else if (value_pos_p(Q2->Constraint[j][INDEX])) {
	    
	    /* Copy -constraint to Mat */
	    value_set_si(Mat->p[1][0],1);
	    for (k=1; k<dim; k++)
	      value_oppose(Mat->p[1][k],Q2->Constraint[j][k]);
	  }
	  else {
	    
	    /* Copy constraint to Mat */
	    value_set_si(Mat->p[1][0],1);
	    for (k=1; k<dim; k++)
	      value_assign(Mat->p[1][k],Q2->Constraint[j][k]);
	  };
	}
	else if (value_pos_p(Q2->Constraint[j][INDEX])) {
	  
	  /* Lower bound -- make an upper bound from it */
	  value_set_si(Mat->p[1][0],1);
	  for(k=1;k<dim;k++)
	    value_oppose(Mat->p[1][k],Q2->Constraint[j][k]);
	}
	else {
	  
	  /* Upper or side bound -- ignore it */
	  continue;
	};
	
#ifdef DEBUG
	fprintf(stdout, "i=%d j=%d M=\n", i+1, j+1);
	Matrix_Print(stdout,P_VALUE_FMT,Mat);
#endif
	
	/* Add Mat to Q and see if anything is made */
	Q3 = AddConstraints(Mat->p[0],2,Q,NbMaxConstrs);

#ifdef DEBUG
	fprintf(stdout, "Q3\n");
	Polyhedron_Print(stdout,P_VALUE_FMT,Q3);
#endif
	
	if (!emptyQ(Q3)) { 
	  Domain_Free(Q3);
	  
#ifdef DEBUG
	  fprintf(stdout, "not empty\n");
#endif
	  res = -1;
	  goto LTQdone;
	}
#ifdef DEBUG
	fprintf(stdout,"empty\n");	
#endif
	Domain_Free(Q3);
      } /* end for j */
    } /* end for i */
    res = 1;
LTQdone:
    Matrix_Free(Mat);
LTQdone2: 
    Domain_Free(Q4);
    Domain_Free(Q1);
    Domain_Free(Q2);
  }
  Domain_Free(Q);
  
#ifdef DEBUG
  fprintf(stdout, "res = %d\n", res);
#endif
  
  return res;
} /* PolyhedronLTQ */

/* GaussSimplify --
   Given Mat1, a matrix of equalities, performs Gaussian elimination.
   Find a minimum basis, Returns the rank.
   Mat1 is context, Mat2 is reduced in context of Mat1
*/
int GaussSimplify(Matrix *Mat1,Matrix *Mat2) {
  
  int NbRows = Mat1->NbRows;
  int NbCols = Mat1->NbColumns;
  int *column_index;
  int i, j, k, n, t, pivot, Rank; 
  Value gcd, tmp, *cp; 
  
  column_index=(int *)malloc(NbCols * sizeof(int));
  if (!column_index) {
    errormsg1("GaussSimplify", "outofmem", "out of memory space\n");
    Pol_status = 1;
    return 0;
  }
  
  /* Initialize all the 'Value' variables */
  value_init(gcd); value_init(tmp);
  
  Rank=0;
  for (j=0; j<NbCols; j++) {		  /* for each column starting at */ 
    for (i=Rank; i<NbRows; i++)		  /* diagonal, look down to find */
      if (value_notzero_p(Mat1->p[i][j])) /* the first non-zero entry    */
	break;	                         
    if (i!=NbRows) {			  /* was one found ? */
      if (i!=Rank)			  /* was it found below the diagonal?*/
	Vector_Exchange(Mat1->p[Rank],Mat1->p[i],NbCols);
      
      /* Normalize the pivot row */
      Vector_Gcd(Mat1->p[Rank],NbCols,&gcd);
      
      /* If (gcd >= 2) */
      value_set_si(tmp,2);
      if (value_ge(gcd,tmp)) {
	cp = Mat1->p[Rank];
        for (k=0; k<NbCols; k++,cp++)
          value_division(*cp,*cp,gcd);		
      }
      if (value_neg_p(Mat1->p[Rank][j])) {
	cp = Mat1->p[Rank];
	for (k=0; k<NbCols; k++,cp++)
	  value_oppose(*cp,*cp);
      }
      /* End of normalize */
      pivot=i;
      for (i=0;i<NbRows;i++)	/* Zero out the rest of the column */
	if (i!=Rank) {
	  if (value_notzero_p(Mat1->p[i][j])) {
	    Value a, a1, a2, a1abs, a2abs;
	    value_init(a); value_init(a1); value_init(a2);
            value_init(a1abs); value_init(a2abs);
            value_assign(a1,Mat1->p[i][j]);
            value_absolute(a1abs,a1);
            value_assign(a2,Mat1->p[Rank][j]); 
            value_absolute(a2abs,a2);
            Gcd(a1abs,a2abs,&a);
	    value_division(a1,a1,a);
	    value_division(a2,a2,a);
	    value_oppose(a1,a1);
	    Vector_Combine(Mat1->p[i],Mat1->p[Rank],Mat1->p[i],a2, 
			   a1,NbCols);
	    Vector_Normalize(Mat1->p[i],NbCols);
	    value_clear(a); value_clear(a1); value_clear(a2);
            value_clear(a1abs); value_clear(a2abs);
          }
	}
      column_index[Rank]=j;
      Rank++;
    }
  } /* end of Gauss elimination */


  if (Mat2) {  /* Mat2 is a transformation matrix  (i,j->f(i,j))....
		  can't scale it because can't scale both sides of -> */
    /* normalizes an affine transformation        */
    /* priority of forms                          */
    /*    1. i' -> i                (identity)    */
    /*    2. i' -> i + constant     (uniform)     */
    /*    3. i' -> constant         (broadcast)   */
    /*    4. i' -> j                (permutation) */
    /*    5. i' -> j + constant     (      )      */
    /*    6. i' -> i + j + constant (non-uniform) */
    for (k=0; k<Rank; k++) {
      j = column_index[k];
      for (i=0; i<(Mat2->NbRows-1);i++) {   /* all but the last row 0...0 1 */
	if ((i!=j) && value_notzero_p(Mat2->p[i][j])) {
	  
	  /* Remove dependency of i' on j */
          Value a, a1, a1abs, a2, a2abs;
	  value_init(a); value_init(a1); value_init(a2);
          value_init(a1abs); value_init(a2abs);
	  value_assign(a1,Mat2->p[i][j]);
	  value_absolute(a1abs,a1);
	  value_assign(a2,Mat1->p[k][j]);
	  value_absolute(a2abs,a2);
	  Gcd(a1abs,a2abs,&a);
	  value_division(a1,a1,a);
	  value_division(a2,a2,a);
	  value_oppose(a1,a1);
	  if (value_one_p(a2)) {
            Vector_Combine(Mat2->p[i],Mat1->p[k],Mat2->p[i],a2,
			   a1,NbCols);
	    
	    /* Vector_Normalize(Mat2->p[i],NbCols); -- can't do T        */
	  } /* otherwise, can't do it without mult lhs prod (2i,3j->...) */
	  value_clear(a); value_clear(a1); value_clear(a2);
          value_clear(a1abs); value_clear(a2abs);
                
	}
        else if ((i==j) && value_zero_p(Mat2->p[i][j])) {
	  
	  /* 'i' does not depend on j */
	  for (n=j+1; n < (NbCols-1); n++) {
	    if (value_notzero_p(Mat2->p[i][n])) { /* i' depends on some n */
	      value_set_si(tmp,1);
              Vector_Combine(Mat2->p[i],Mat1->p[k],Mat2->p[i],tmp,
			     tmp,NbCols);
	      break;
	    }  /* if 'i' depends on just a constant, then leave it alone.*/
	  }
        }
      }
    }
    
    /* Check last row of transformation Mat2 */
    for (j=0; j<(NbCols-1); j++)
      if (value_notzero_p(Mat2->p[Mat2->NbRows-1][j])) {
	errormsg1("GaussSimplify", "corrtrans", "Corrupted transformation\n");
	break;
      }
    
    if (value_notone_p(Mat2->p[Mat2->NbRows-1][NbCols-1])) {
      errormsg1("GaussSimplify", "corrtrans", "Corrupted transformation\n");
    }
  }
  value_clear(gcd); value_clear(tmp);
  free(column_index);
  return Rank;
} /* GaussSimplify */

/* 
 * Topologically sort 'n' polyhdera and return 0 on failure, otherwise return 
 * 1 on success. Here 'L' is a an array of pointers to polyhedra, 'n' is the 
 * number of polyhedra, 'index' is the level to consider for partial ordering
 * 'pdim' is the parameter space dimension, 'time' is an array of 'n' integers
 * to store logical time values, 'pvect', if not NULL, is an array of 'n' 
 * integers that contains a permutation specification after call and MAXRAYS is
 * the workspace size for polyhedra operations. 
 */
int PolyhedronTSort (Polyhedron **L,unsigned int n,unsigned int index,unsigned int pdim,int *time,int *pvect,unsigned int MAXRAYS) {
 
  unsigned int const nbcells = ((n*(n-1))>>1);    /* Number of memory cells 
						     to allocate, see below */
  int *dag;  /* The upper triangular matrix */
  int **p;   /* Array of matrix row addresses */
  unsigned int i, j, k;
  unsigned int t, nb, isroot;

  if (n<2) return 0;

  /* we need an upper triangular matrix (example with n=4):

     . o o o
     . . o o     . are unuseful cells, o are useful cells
     . . . o
     . . . .
     
     so we need to allocate (n)(n-1)/2 cells
     - dag will point to this memory.
     - p[i] will point to row i of the matrix
     p[0] = dag - 1 (such that p[0][1] == dag[0])
     p[1] = dag - 1 + (n-1)
     p[2] = dag - 1 + (n-1) + (n-2)
     ...
     p[i] = p[i-1] + (n-1-i)
  */

  /* malloc n(n-1)/2 for dag and n-1 for p (node n does not have any row) */
  dag = (int *) malloc(nbcells*sizeof(int));
  if (!dag) return 0;
  p = (int **) malloc ((n-1) * sizeof(int *));
  if (!p) { 
    free(dag); return 0; 
  }

  /* Initialize 'p' and 'dag' */
  p[0] = dag-1;
  for (i=1; i<n-1; i++)
    p[i] = p[i-1] + (n-1)-i;
  for (i=0; i<nbcells; i++)
    dag[i] = -2;      /* -2 means 'not computed yet' */
  for (i=0; i<n; i++) time[i] = -1;

  /* Compute the dag using transitivity to reduce the number of */
  /*   PolyhedronLTQ calls.                                     */
  for (i=0; i<n-1; i++) {
    POL_ENSURE_FACETS(L[i]);
    POL_ENSURE_VERTICES(L[i]);
    for (j=i+1; j<n; j++) {
      if (p[i][j] == -2) /* not computed yes */
	p[i][j] = PolyhedronLTQ(L[i], L[j], index, pdim, MAXRAYS);
      if (p[i][j] != 0) {
	
	/* if p[i][j] is 1 or -1, look for -p[i][j] on the same row:
	   p[i][j] == -p[i][k] ==> p[j][k] = p[i][k] (transitivity)
	   note: p[r][c] == -p[c][r], use this to avoid reading or writing
	   under the matrix diagonal
	*/  
	   
	/* first, k<i so look for p[i][j] == p[k][i] (i.e. -p[i][k]) */
	for (k=0; k<i; k++)
	  if (p[k][i] == p[i][j]) p[k][j] = p[k][i];
	
	/* then, i<k<j so look for p[i][j] == -p[i][k] */
	for (k=i+1; k<j; k++)
	  if (p[i][k] == -p[i][j]) p[k][j] = -p[i][k];
	
	/* last, k>j same search but */
	for (k=j+1; k<n; k++)
	  if (p[i][k] == -p[i][j]) p[j][k] = p[i][k];
      }
    }
  }
  
  /* walk thru the dag to compute the partial order (and optionally
     the permutation)
     Note: this is not the fastest way to do it but it takes
     negligible time compared to a single call of PolyhedronLTQ !
     Each macro-step (while loop) gives the same logical time to all
     found roots and optionally add these nodes in the permutation vector
  */ 
  
  t = 0; /* current logical time, assigned to current roots and
	    increased by 1 at the end of each step */
  nb = 0; /* number of processed nodes (have been given a time) */
  while (nb<n) {
    for (i=0; i<n; i++) { /* search for roots */
      /* for any node, if it's not already been given a logical time
	 then walk thru the node row; if we find a 1 at some column j,
	 it means that node j preceeds the current node, so it is not
	 a root */
      if (time[i]<0) {
	isroot = 1; /* assume that it is until it is definitely not */
	/* first search on a column */
	for (j=0; j<i; j++) {
	  if (p[j][i]==-1) { /* found a node lower than it */
	    isroot = 0; break;
	  }
	}
	if /*still*/ (isroot)
	  for (j=i+1; j<n; j++) {
	    if (p[i][j]==1) { /* greater than this node */
	      isroot = 0; break;
	    }
	  }
	if (isroot) { /* then it definitely is */
	  time[i] = t; /* this node gets current logical time */
	  if (pvect)
	    pvect[nb] = i+1; /* nodes will be numbered from 1 to n */
	  nb++; /* one more node processed */
	}
      }
    }
    /* now make roots become neutral, i.e. non comparable with all other nodes */
    for (i=0; i<n; i++) {
      if (time[i]==t) {
	for (j=0; j<i; j++)
	  p[j][i] = 0;
	for (j=i+1; j<n; j++)
	  p[i][j] = 0;
      }
    }
    t++; /* ready for next set of root nodes */
  }
  
  free (p);   /* let's be clean, collect the garbage */
  free (dag);
  return 1;
} /* PolyhedronTSort */




