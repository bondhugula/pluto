#include <stdlib.h>
#include <polylib/polylib.h>

/* nota bene: on stocke les matrices par lignes */
/* nota bene: matrices are stored in row major order */

/*------------------------------------------------------------------------
	change les signes de la ligne i de la matrice A
        change the sign of row i of matrix A
------------------------------------------------------------------------*/

static void moins_l(Value *a,int i,int n,int p) {
  
  int k;
  Value *c;
  
  c=a+(i-1)*p;
  
  for(k=1; k<=p; k++) {
    value_oppose(*c,*c);
    c++;
  }
  return;
} /* moins_l */

/*------------------------------------------------------------------------
	change les signes de la colonne i de la matrice A
        change the sign of column i of matrix A
------------------------------------------------------------------------*/

static void moins_c(Value *a,int i,int n,int p) {
  
  int k;
  Value *c;
  
  c=a+(i-1);
  
  for(k=1;k<=n;k++) {
    value_oppose(*c,*c);
    c=c+p;
  }
  return;
} /* moins_c */

/*------------------------------------------------------------------------
	echange les lignes i et j de la matrice A
        exchange row i and j of matrix A
------------------------------------------------------------------------*/

static void echange_l(Value *a,int i,int j,int n,int p) {	
  
  int k;
  Value s;
  Value *c1,*c2;
  
  value_init(s);
  c1=a+(i-1)*p;
  c2=a+(j-1)*p;
  
  for(k=1;k<=p;k++) {
    value_assign(s,*c1);
    value_assign(*c1,*c2);
    value_assign(*c2,s);
    c1++;
    c2++;
  }
  value_clear(s);
  return;
} /* echange_l */ 

/*------------------------------------------------------------------------
	echange les colonnes i et j de la matrice A
        exchange columns i and j of matrix A
------------------------------------------------------------------------*/

static void echange_c(Value *a,int i,int j,int n,int p) {	
  
  int k;
  Value s;
  Value *c1,*c2;

  value_init(s);
  c1=a+(i-1);
  c2=a+(j-1);
  
  for(k=1;k<=n;k++) {	
    value_assign(s,*c1);
    value_assign(*c1,*c2);
    value_assign(*c2,s);
    c1=c1+p;
    c2=c2+p;
  }
  value_clear(s);
  return;
} /* echange_c */

/*------------------------------------------------------------------------
	A est une matrice de taille n*p
	on ajoute a la ligne i, x fois la ligne j
        A is a matrix of size n*p
        we add x times the row j to the row i
------------------------------------------------------------------------*/

static void ligne(Value *a,int i,int j,Value x,int n,int p) {	
 
  int k;
  Value *c1,*c2;

  c1=a+(i-1)*p;
  c2=a+(j-1)*p;
  
  for(k=1;k<=p;k++) {

    value_addmul(*c1, x, *c2);
    c1++;
    c2++;
  }
  return;
} /* ligne */

/*------------------------------------------------------------------------
	A est une matrice de taille n*p
	on ajoute a la colonne i, x fois la colonne j
        A is a matrix of size n*p
        we add x times the column j to the column i

------------------------------------------------------------------------*/

static void colonne(Value *a,int i,int j,Value x,int n,int p) {	
  
  int k;
  Value *c1,*c2;

  c1=a+(i-1);
  c2=a+(j-1);
  
  for(k=1;k<=n;k++) {
    value_addmul(*c1, x, *c2);
    c1=c1+p;
    c2=c2+p;
  }
  return;
} /* colonne */
	
/*----------------------------------------------------------------------
	trouve le numero de colonne du plus petit element non nul de 
	la ligne q, d'indice superieur a q, de la matrice A,  mais
	renvoie zero si tous les elements sont nuls sauf le qieme.

        find the column number (greater than q) of the smallest non
        null element of row q . it returns 
        0 if all the last elements of row q are null except the qth

----------------------------------------------------------------------*/

static int petit_l(Value *a,int n,int p,int q) {
  
  int numero=0, k, tousnuls;
  Value minus, comp;
  Value *c;

  value_init(minus); value_init(comp);
  c=a+(q-1)*p+(p-1);
  tousnuls=1;
  for(k=p;k>q;k--) {
    value_absolute(comp,*c);
    if (value_notzero_p(comp)) {	
      if (tousnuls==1) {	
	value_assign(minus,comp);
	numero=k;
	tousnuls=0;
      }
      else if (value_ge(minus,comp)) {
	value_assign(minus,comp);
	numero=k;
      }
    }
    c--;
  }  
  if (tousnuls==1) {
    value_clear(minus); value_clear(comp);
    return(0);
  }  
  else  {
    value_absolute(comp,*c);
    if ((value_notzero_p(comp))&&(value_ge(minus,comp))) {
      value_clear(minus); value_clear(comp);
      return(q);
    }  
    else {
      value_clear(minus); value_clear(comp);
      return(numero);
    }  
  }
} /* petit_l */ 
	
/*----------------------------------------------------------------------
	trouve le numero de ligne du plus petit element non nul de 
	la colonne q, d'indice superieur a q, de la matrice A,  mais
	renvoie zero si tous les elements sont nuls sauf le qieme.

        find the row number (greater than q) of the smallest non
        null element of column q . it returns 
        0 if all the last elements of column q are null except the qth

----------------------------------------------------------------------*/

static int petit_c(Value *a,int n,int p,int q) {
  
  int numero=0, k, tousnuls;  
  Value minus, comp;
  Value *c;
  
  value_init(minus); value_init(comp);
  c = a+q-1+p*(n-1);
  tousnuls=1;
  for(k=n;k>q;k--) {
    value_absolute(comp,*c);
    if (value_notzero_p(comp)) {	
      if (tousnuls==1) {	
	value_assign(minus,comp);
	numero=k;
	tousnuls=0;
      }
      else if (value_ge(minus,comp)) {
	value_assign(minus,comp);
	numero=k;
      }
    }
    c=c-p;
  }
  
  if (tousnuls==1) {
    value_clear(minus); value_clear(comp);
    return(0);
  }  
  else {
    value_absolute(comp,*c);
    if ((value_notzero_p(comp)) && (value_ge(minus,comp))) {
      value_clear(minus); value_clear(comp);
      return(q);
    }  
    else {
      value_clear(minus); value_clear(comp);
      return(numero);
    }  
  }
} /* petit_c */
	
/*-----------------------------------------------------------------------
	A est initialisee a une matrice "identite" de taille n*p
	A is initialized to an "identity" matrix of size n*p
-----------------------------------------------------------------------*/

static void identite(Value *a,int n,int p) {
  
  int i,j;
  Value *b;
  
  b = a;
  for(i=1;i<=n;i++) {
    for(j=1;j<=p;j++) {
      if (i==j) 
	value_set_si(*b,1);
      else 
	value_set_si(*b,0);
      b++;
    }
  }
  return;
} /* identite */

/*-----------------------------------------------------------------------
  transpose la sous-matrice de A dont les indices sont
  superieurs ou egaux a q
  transposition of the sub-matrix of A composed of the elements 
  whose indices are greater or equal to q
-----------------------------------------------------------------------*/

static void transpose(Value *a,int n,int q) {
  
  int i;
  Value val;
  Value *b,*c;

  value_init(val);
  if (q<n) {
    b=a+(q-1)+(q-1)*n;
    c=b;
    for(i=q+1;i<=n;i++) {	
      b++;
      c=c+n;
      value_assign(val,*b);
      value_assign(*b,*c);
      value_assign(*c,val);
    }
    transpose(a,n,q+1);
  }
  value_clear(val);
  return;
} /* transpose */

/*----------------------------------------------------------------------
	trouve le numero de colonne du premier element non divisible 
	par val dans la sous-matrice d'indices > q mais renvoie 0 sinon
        
        find the column number of the first non-divisible by val element
        in the sub-matrix of a composed of the elements of indices >q.
        Return 0 is there is no such element
----------------------------------------------------------------------*/

static int encore(Value *a,int n,int p,int q,Value val) {
  
  int k,l;
  Value comp, tmp;
  Value *c;

  if ((q>=n)||(q>=p)) return(0);
  
  value_init(comp); value_init(tmp); 
  c=a+q*p+q;
  k=q+1;
  l=q+1;
  while (k<=n) {
    value_absolute(comp,*c);
    if (value_zero_p(val)) {
      if (value_notzero_p(comp)) { 
	value_clear(comp);
	value_clear(tmp);
	return(l);
      }	
    }
    else {
      value_modulus(tmp,comp,val);
      if (value_notzero_p(tmp)) {
	value_clear(comp);
	value_clear(tmp);
	return(l);
      }
    }  
    if (l==p) {
      k=k+1;
      l=q+1;
      c=c+q+1;
    }
    else {
      l++;
      c++;
    }
  }
  value_clear(comp);
  value_clear(tmp);
  return(0);
} /* encore */

/*----------------------------------------------------------------------
	transforme la matrice A (de taille n*p), a partir de la position 
	q*q en sa forme de smith, les matrices b et c subissant les memes 
	transformations respectivement sur les lignes et colonnes

        transform the matrix A (size n*p), from position (q,q) into 
        its smith normal form. the same transformations are applied to 
        matrices b and c, respectively on rows and on columns
----------------------------------------------------------------------*/

static void smith(Value *a,Value *b,Value *c,Value *b_inverse,Value *c_inverse,int n,int p,int q) {
  
  int i,k,fini;
  Value x, pivot, tmp, x_inv;
  Value *f;
  
  value_init(pivot); value_init(tmp); 
  value_init(x); value_init(x_inv);
  
  if ((q<=n)&&(q<=p)) {
    fini = 1;
    
    while (fini!=0) {
      i=petit_c(a,n,p,q);
      while (i!=0) {
	if (i!=q) {
	  echange_l(a,i,q,n,p);
	  echange_l(b,i,q,n,n);
	  echange_c(b_inverse,i,q,n,n);
	}	
	f=a+(q-1)+(q-1)*p;
	value_assign(pivot,*f);
	if (value_neg_p(pivot)) {
	  moins_l(a,q,n,p);
	  moins_l(b,q,n,n);
	  moins_c(b_inverse,q,n,n);
	  value_oppose(pivot,pivot);
	}	
	for(k=q+1;k<=n;k++) {	
	  f=f+p;
	  if (value_notzero_p(*f)) {
	    value_division(x,*f,pivot);
	    value_modulus(tmp,*f,pivot);
	    if (value_neg_p(tmp)) 
	      value_decrement(x,x);
	    value_oppose(x_inv,x);
	    ligne(a,k,q,x_inv,n,p);
	    ligne(b,k,q,x_inv,n,n);
	    colonne(b_inverse,q,k,x,n,n);	
	  }		
	}	
	i=petit_c(a,n,p,q);
      }	      
      fini=0;
      i=petit_l(a,n,p,q);
      while (i!=0) {
	if (i!=q) {
	  echange_c(a,i,q,n,p);
	  echange_c(c,i,q,p,p);
	  echange_l(c_inverse,i,q,p,p);
	  fini=1;
	}	
	f=a+(q-1)+(q-1)*p;
	value_assign(pivot,*f);
	if (value_neg_p(pivot)) {
	  moins_c(a,q,n,p);
	  moins_c(c,q,p,p);
	  moins_l(c_inverse,q,p,p);
	  value_oppose(pivot,pivot);
	}	
	for(k=q+1;k<=p;k++) {	
	  f++;
	  if (value_notzero_p(*f)) {
	    value_division(x,*f,pivot);
	    value_modulus(tmp,*f,pivot);
	    if (value_neg_p(tmp)) 
	      value_decrement(x,x);
	    value_oppose(x_inv,x);
	    colonne(a,k,q,x_inv,n,p);
	    colonne(c,k,q,x_inv,p,p);	
	    ligne(c_inverse,q,k,x,p,p);
	  }		
	}	
	i=petit_l(a,n,p,q);
      }	
    }
    value_assign(pivot,*(a+(q-1)+(q-1)*p));
    if (value_neg_p(pivot)) {
      moins_l(a,q,n,p);
      moins_l(b,q,n,n);
      moins_c(b_inverse,q,n,n);
      value_oppose(pivot,pivot);
    }
    
    i=encore(a,n,p,q,pivot);
    if (i!=0) {
      value_set_si(x,1);
      value_set_si(x_inv,-1);
      colonne(a,q,i,x,n,p);
      colonne(c,q,i,x,p,p);
      ligne(c_inverse,i,q,x_inv,p,p);
      smith(a,b,c,b_inverse,c_inverse,n,p,q);
    }	
    else smith(a,b,c,b_inverse,c_inverse,n,p,q+1);
  }
  value_clear(pivot); value_clear(tmp);
  value_clear(x); value_clear(x_inv);
  
  return;
} /* smith */ 

/*----------------------------------------------------------------------
	decompose la matrice A en le produit d'une matrice B
	unimodulaire et d'une matrice triangulaire superieure
	D est l'inverse de B
        Decompose the matrix A in the product of a matrix B unimodular
        and a matrix upper triangular, D is the inverse of B
----------------------------------------------------------------------*/

static void hermite(Value *a,Value *b,Value *d,int n,int p,int q) {
  
  int i,k;
  Value x, pivot, x_inv, tmp;
  Value *c1;

  value_init(pivot);  value_init(tmp);
  value_init(x); value_init(x_inv);

  if ((q<=p)&&(q<=n)) {
    i=petit_c(a,n,p,q);    
    while (i!=0) {
      if (i!=q) {
	echange_l(a,i,q,n,p);
	echange_c(b,i,q,n,n);
	echange_l(d,i,q,n,n);
      }
      
      c1=a+(q-1)+(q-1)*p;
      value_assign(pivot,*c1);
      if (value_neg_p(pivot)) {
	moins_l(a,q,n,p);
	moins_l(d,q,n,n);
	moins_c(b,q,n,n);
	value_oppose(pivot,pivot);
      }      
      for(k=q+1;k<=n;k++) {	
	c1=c1+p;
	if (value_notzero_p(*c1)) {
	  value_division(x,*c1,pivot);
	  value_modulus(tmp,*c1,pivot);
	  if (value_neg_p(tmp)) 
	    value_decrement(x,x);
	  value_oppose(x_inv,x);
	  ligne(a,k,q,x_inv,n,p);
	  colonne(b,q,k,x,n,n);	
	  ligne(d,k,q,x_inv,n,n);
	}		
      }      
      i=petit_c(a,n,p,q);	
    }    
    c1=a+(q-1);
    value_assign(pivot,*(c1+(q-1)*p));
    if (value_neg_p(pivot)) {
      moins_l(a,q,n,p);
      moins_l(d,q,n,n);
      moins_c(b,q,n,n);
      value_oppose(pivot,pivot);
    }    
    if (value_notzero_p(pivot)) {
      for(k=1;k<q;k++) {
	if (value_notzero_p(*c1)) {
	  value_division(x,*c1,pivot);
	  value_modulus(tmp,*c1,pivot);
	  if (value_neg_p(tmp)) 
	    value_decrement(x,x);
	  value_oppose(x_inv,x);
	  ligne(a,k,q,x_inv,n,p);
	  colonne(b,q,k,x,n,n);	
	  ligne(d,k,q,x_inv,n,n);
	}
	c1=c1+p;		
      }
    }    
    hermite(a,b,d,n,p,q+1);
  }	
  value_clear(pivot); value_clear(tmp);
  value_clear(x); value_clear(x_inv);
  return;
} /* hermite */
	
/*----------------------------------------------------------------------
  B est une g_base de dimension n dont les p premiers vecteurs
  colonnes sont colineaires aux vecteurs colonnes de A
  C est l'invers de B
  
  B is an integral (g_base ?) of dimension n whose p first columns
  vectors are proportionnal to the column vectors of A. C is the
  inverse of B.
  ----------------------------------------------------------------------*/

/** Convert PolmattoDarmat :
***  This function converts the matrix of a Polylib to a int * as necessary
***    for the functions in Darte's implementation.
***
***  Input  : A Polylib Matrix ( a pointer to it );
***  Output : An int * with the matrix copied into it
**/
static Value *ConvertPolMattoDarMat(Matrix *A ) {
  
  int i,j;
  Value *result;
  
  result = (Value *)malloc(sizeof(Value) * A->NbRows * A->NbColumns);
  for(i=0;i<A->NbRows * A->NbColumns;i++)
    value_init(result[i]);

  for (i = 0; i < A->NbRows; i++)
    for (j = 0 ;  j < A->NbColumns; j++)
      value_assign(result[i*A->NbColumns + j],A->p[i][j]);
  return result;
} /* ConvertPolMattoDarMat */

/** Convert DarmattoPolmat
***   This function converts the matrix from Darte representation to a matrix
***    in PolyLib.
***
***   Input  : The matrix (a pointer to it),  number of Rows, number of columns
***   Output : The matrix of the PolyLib 
***
**/

static Matrix *ConvertDarMattoPolMat (Value *A, int NbRows, int NbCols) {
  
  int i,j;
  Matrix *result;
  
  result = Matrix_Alloc(NbRows, NbCols);
  
  for (i = 0; i < NbRows; i ++)
    for (j = 0; j < NbCols; j ++)
      value_assign(result->p[i][j],A[i*NbCols + j]);  
  return result;
} /* ConvertDarMattoPolMat */

/**
*** Smith : This function takes a Matrix A of dim n * l as its input 
***         and returns the three matrices U, V and Product such that
***         A = U * Product * V, where U is an unimodular matrix of dimension
***         n * n and V is an unimodular matrix of dimension l * l.
***         Product is a diagonal matrix of dimension n * l.
***     
***         We use Alan Darte's implementation of Smith for computing
***         the Smith Normal Form 
**/
void Smith(Matrix *A, Matrix **U, Matrix **V, Matrix **Product)
{ 
  int i;
  Matrix *u, *v;
  
  u = Identity(A->NbRows);
  v = Identity(A->NbColumns);

  *U = Identity(A->NbRows);
  *V = Identity(A->NbColumns);
  
  *Product = Matrix_Copy(A);
  smith((*Product)->p_Init, u->p_Init, v->p_Init, (*U)->p_Init, (*V)->p_Init,
	A->NbRows, A->NbColumns, 1);

  Matrix_Free(u);
  Matrix_Free(v);
} /* Smith */

/** Hermite :
***   This function takes a Matrix as its input and finds its HNF 
***    ( Left form )
***
***   Input  : A Matrix A (The Matrix A is not necessarily a Polylib matrix.
***             It is just a matrix as far as Hermite is concerned. It
***             does not even check if the matrix is a Polylib matrix or not) 
***   Output : The Hnf matrix H and the Unimodular matrix U such that
***               A = H * U.
***
***   We use Alan Darte's implementation of Hermite to compute the HNF.
***   Alan Darte's implementation computes the Upper Triangular HNF.
***   So We work on the fact that if A = H * U then 
***       A (transpose) = U(transpose) * H (transpose)
***       There are a set of interface functions written in Interface.c
***        which convert a matrix from Polylib representationt to that of
***        Alan Darte's and vice versa.
***     
***   This Function Does the Following
***    Step 1 : Given the matrix A it finds its Transpose.
***    Step 2 : Finds the HNF (Right Form) using Alan Darte's Algorithm.
***    Step 3 : The H1 and U1 obtained in Step2 are both Transposed to get
***              the actual H and U such that A = HU.
**/
void Hermite (Matrix *A, Matrix **H, Matrix **U) {
  
  int i;
  Matrix *transpose, *tempH, *tempU ;
  Value *darte_matA, *darte_identite, *darte_id_inv;
  
  transpose = Transpose(A);
  darte_matA = ConvertPolMattoDarMat(transpose);
  
  darte_identite = (Value *)malloc(sizeof(Value) * A->NbColumns* A->NbColumns);
  darte_id_inv = (Value *) malloc(sizeof(Value) * A->NbColumns*A->NbColumns);
  for (i=0; i< A->NbColumns * A->NbColumns; i++)
    value_init(darte_identite[i]);
  for (i=0; i< A->NbColumns * A->NbColumns; i++)
    value_init(darte_id_inv[i]);
  
  identite(darte_identite, transpose->NbRows, transpose->NbRows);
  identite(darte_id_inv, transpose->NbRows, transpose->NbRows);  
  hermite(darte_matA, darte_identite, darte_id_inv,A->NbColumns,A->NbRows, 1); 
  Matrix_Free (transpose);
  transpose = ConvertDarMattoPolMat(darte_matA, A->NbColumns, A->NbRows);
  tempU = ConvertDarMattoPolMat(darte_identite, A->NbColumns, A->NbColumns);
  
  /* Finding H Transpose */
  tempH = Transpose(transpose);  
  Matrix_Free(transpose);
  
  /* Finding U Transpose */
  transpose = Transpose(tempU);
  
  H[0] = Matrix_Copy(tempH);
  U[0] = Matrix_Copy(transpose);
    
  Matrix_Free (tempH);
  Matrix_Free (transpose);
  Matrix_Free (tempU);
  
  for (i=0; i< A->NbRows * A->NbColumns; i++)
    value_clear(darte_matA[i]);
  for (i=0; i< A->NbColumns * A->NbColumns; i++)
    value_clear(darte_identite[i]);
  for (i=0; i< A->NbColumns * A->NbColumns; i++)
    value_clear(darte_id_inv[i]);
  free (darte_matA);
  free (darte_identite);
  free (darte_id_inv);
  return;
} /* Hermite */


