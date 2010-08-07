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
int PolyhedronLTQ (Polyhedron *Pol1,Polyhedron *Pol2,int NbMaxConstrs) { 
  
  int res, dim, i, j, k;
  Polyhedron *Q1, *Q2, *Q3, *Q;
  Matrix *Mat;

#define INDEX 1

  if (Pol1->next || Pol2->next) {
    errormsg1("PolyhedronLTQ", "compoly", "Can only compare polyhedra");
    return 0;
  }
  if (Pol1->Dimension != Pol2->Dimension) {
    errormsg1("PolyhedronLTQ","diffdim","Polyhedra are not same dimension");
    return 0;
  }
  dim = Pol1->Dimension+2;
  
#ifdef DEBUG
  fprintf(stdout, "P1\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Pol1);
  fprintf(stdout, "P2\n");
  Polyhedron_Print(stdout,P_VALUE_FMT,Pol2);
#endif
  
  /* Create the Line to add */
  Mat = Matrix_Alloc(1,dim);
  Vector_Set(Mat->p_Init,0,dim);
  value_set_si(Mat->p[0][INDEX],1);  /* line in INDEX dimension */
  
  Q1 = AddRays(Mat->p[0],1,Pol1,NbMaxConstrs);
  Q2 = AddRays(Mat->p[0],1,Pol2,NbMaxConstrs);

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
  int i, j, k, n, pivot, Rank; 
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
	    Value a, a1, a2;
	    value_init(a); value_init(a1); value_init(a2);
	    value_absolute(a1,Mat1->p[i][j]);
	    value_absolute(a2,Mat1->p[Rank][j]);
	    Gcd(a1,a2,&a));
	    value_division(a1,a1,a);
	    value_division(a2,a2,a);
	    value_oppose(a1,a1);
	    Vector_Combine(Mat1->p[i],Mat1->p[Rank],Mat1->p[i],a2, 
			   a1,NbCols);
	    Vector_Normalize(Mat1->p[i],NbCols);
	    value_clear(a); value_clear(a1); value_clear(a2);
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
	  Value a, a1, a2;
	  value_init(a); value_init(a1); value_init(a2);
	  value_absolute(a1,Mat2->p[i][j]);
	  value_absolute(a2,Mat2->p[k][j]);
	  Gcd(a1,a2,&a);
	  value_division(a1,a1,a);
	  value_division(a2,a2,a);
	  value_oppose(a1,a1);
	  if (value_one_p(a2)) {
	    Vector_Combine(Mat2->p[i],Mat1->p[k],Mat2->p[i],a2,
			   a1,NbCols);
	    
	    /* Vector_Normalize(Mat2->p[i],NbCols); -- can't do T        */
	  } /* otherwise, can't do it without mult lhs prod (2i,3j->...) */
	  value_clear(a); value_clear(a1); value_clear(a2);
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
    tmp = Constraints2Polyhedron(a,600);
    Matrix_Free(a);
    if (!last) A = last = tmp;
    else {
      last->next = tmp;
      last = tmp;
    }
    }

    if (nbMat) {
      a = Matrix_Read(); 
    }
    fgets(s,128,stdin);
    nbPol = nbMat = 0;
    while ( (*s=='#') ||
        ((sscanf(s, "D %d", &nbPol)<1) && (sscanf(s, "M %d", &nbMat)<1)) )
      fgets(s, 128, stdin);

    for (i=0, B=last=(Polyhedron *)0; i<nbPol; i++) {
      b = Matrix_Read();
      tmp = Constraints2Polyhedron(b,200);
      Matrix_Free(b);
      if (!last) B = last = tmp;
      else {
	last->next = tmp;
	last = tmp;
      }
    }

    if (nbMat) {
      b = Matrix_Read();
    }    
    fgets(s, 128, stdin);
    while ((*s=='#') || (sscanf(s, "F %d", &func)<1))
      fgets(s, 128, stdin);
    
    switch (func) {
    case 1:
      C = DomainUnion(A, B, 200);
      D = DomainConvex(C, 200);
      d = Polyhedron2Constraints(D);
      Matrix_Print(stdout,P_VALUE_FMT,d);
      Matrix_Free(d);
      Domain_Free(C);
      Domain_Free(D);
      break;
    case 2:
      D = DomainSimplify(A, B, 200);
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
      
      C = Polyhedron_Image(A, b, 400);
      Polyhedron_Print(stdout,P_VALUE_FMT,C);
      break;
    case 8:
      
      printf("%s\n",
	     Polyhedron_Not_Empty(A,B,600) ? "Not Empty" : "Empty");
      break;
    case 9:
      
      i = PolyhedronLTQ(A,B,600);
      printf("%s\n",
	     i==-1 ? "A<B" : i==1 ? "A>B" : i==0 ? "A><B" : "error");
      i = PolyhedronLTQ(B,A,600);
      printf("%s\n",
	     i==-1 ? "A<B" : i==1 ? "A>B" : i==0 ? "A><B" : "error");
      break;
    default:
      printf("? unknown function\n");
    }
    
    Domain_Free(A);
    Domain_Free(B);
    
    return 0;
}


