
#include <stdlib.h>
#include <polylib/polylib.h>
#define WS 200 

/* #define UE_DEBUG */


typedef struct _Polyhedron_union { 
	Enumeration *pt;
        struct _Polyhedron_union *next;} Polyhedron_union;

static int ppcm1 (int a, int b);
Matrix *CalcBase( Matrix *R);


/**** ker.c ****/

static void Soustraire_ligne(Matrix *R, int l1, int l2, int piv );
static int existepivot( Matrix *R, int l );
static void swap_line(Matrix *R, int l1, int l2);


/* Caclule une base du noyau de l'application definie
par la matrice _carree_ R (completee par des 0..0) */

/* /!\ R est alteree par cette fonction. */

Matrix *CalcBase( Matrix *R )
{
	Matrix *B/*,*BTran*/;	/* matrice des vect. de la base */
	int i, j;
	Value p;
	int l;		/* ligne en cours de `pivotage' */
	int lnn;
	int dimbase;
	int u;		/* vecteur util en cours de calcul */
	Value  som;
	int c;
	/* diagonalisation de R : algo du pivot */
	/* avec conservation des pivots nuls */
	for( l=0 ; l<R->NbRows ; ++l )
	{
		/* cherche le prochain pivot non nul */
		if( (lnn=existepivot(R,l)) != -1 )
		{
			swap_line( R, l, lnn );

			/* met la colonne dessous a 0 */
			for( j=l+1 ; j<R->NbRows ; ++j )
				if(value_notzero_p( R->p[j][l]) )
					Soustraire_ligne( R, l, j, l );

			/* et celle dessus aussi */
			for( j=0 ; j<l ; ++j )
				if( value_notzero_p(R->p[j][l]) )
					Soustraire_ligne( R, l, j, l );

		}
		/* s'ils sont tous nuls on passe au suivant */
		/* en laissant cette ligne inchangee. */
	}

	dimbase = 0;
	for( l=0 ; l<R->NbRows ; ++l )
		if( value_zero_p(R->p[l][l]) )
			++dimbase;     /* = nb de 0 dans la diagonale */

	B = Matrix_Alloc( dimbase, R->NbRows );

	l=0;
	for( u=0 ; u<dimbase ; ++u )
	{
		while(value_notzero_p( R->p[l][l]) )
			++l;

	/* calcul effectif de ce vecteur : chaque coord en commencant par le bas */
		for( i=R->NbRows-1 ; i>l ; --i )
			value_set_si(B->p[u][i], 0);
		value_set_si(B->p[u][l], 1);
		for( i=l-1 ; i>=0 ; --i )	/* i=chaque coord du vecteur */
		{
			if(value_zero_p( R->p[i][i]) )
			  /* on a deux 0... un seul dans */
			  /* ce vect util suffit */
				value_set_si(B->p[u][i],0); 
			else
			{
				/* somme des coef deja calcules * coef dans la matrice */
				value_set_si(som,0);
				for( c=l ; c>i ; --c )
				{
				  value_addmul(som, R->p[i][c], B->p[u][c]);
				  value_multiply(B->p[u][c] ,B->p[u][c] , R->p[i][i]);	
				}
				value_oppose(B->p[u][i] , som );
			}
		}
		/* reste a faire le pgcd du vecteur et a l'orienter */
		value_set_si(p,0);
		for( i=0 ; i<R->NbRows ; ++i )
		  Gcd( p, B->p[u][i], &p );
		if( value_zero_p(p))
			value_set_si(p,1);
		for( i=0 ; i<R->NbRows && value_zero_p(B->p[u][i]); ++i )
			;
		if( i<R->NbRows )
			if( value_neg_p(B->p[u][i]) )
		      value_oppose( p,p );

		for( i=0 ; i<R->NbRows ; ++i )
			value_division(B->p[u][i],B->p[u][i], p);

		/* incrementer le compteur de lignes */
		++l;
	}
	return B;
}


/* fonction qui calcule les vect generateurs de l'espace vectoriel
   contenant le polyhedre P */
/*
Matrix *CalcPolyhedronBase( Polyhedron *P )
{
  Matrix *R;
  Matrix *B;
  int n, lines;
  
  if( emptyQ(P) )
    return( Matrix_Alloc( 0, P->Dimension ) ); */ /* pas de vect. gen */ /*
  
  R = Matrix_Alloc( P->Dimension, P->Dimension );
  
  */ /* recopie le 'lineality space' du polyedre dans la matrice R */ /*
  for( lines=0,n=0 ; n<P->NbConstraints ; ++n )
    {
      if( P->Constraint[n][0]==0 )
	*/ /* c'est une direction definissant un ss-espace */ /*
	{
	  memcpy( &R->p[lines][0], &P->Constraint[n][1],
		  sizeof(int)*P->Dimension );
	  ++lines;
	}
    }
  */ /* remplit le reste de 0..0 */ /*
  for( ; lines<R->NbRows ; ++lines )
    memset( &R->p[lines][0], 0, sizeof(int)*P->Dimension );
  
  B = CalcBase( R );
  Matrix_Free( R );
  
  return( B );
}*/


/* fonction qui calcule l'espace vectoriel non parametre de dimension dim
	contenant le polyhedre parametre P */
/* les _egalites_ sont stockees par la polylib sous forme triangulaire
superieure donc il suffit de prendre les premieres. */

/*Matrix *CalcEVPolyhedronNP( Polyhedron *P, int dim )
{
	Matrix *R;
	Matrix *B;
	int n, lines;

	if( emptyQ(P) ) */ /* pas de vect. gen */ /*
	{
	  B = Matrix_Alloc( 1, dim );	*/ /* on ne peut pas allouer 0 lignes ! */ /*
		B->NbRows = 0;
		return( B );
	}

	R = Matrix_Alloc( dim, dim );

*/ /* recopie le 'lineality space' du polyedre dans la matrice R */ /*
	for( lines=0,n=0 ; n<P->NbConstraints && lines<dim ; ++n )
	{
		if( P->Constraint[n][0]==0 )
		  */ /* c'est une direction definissant un ss-espace */ /*
		{
			memcpy( &R->p[lines][0], &P->Constraint[n][1],
					sizeof(int)*P->Dimension );
			++lines;
		}
	}
*/ /* remplit le reste de 0..0 */ /*
	for( ; lines<R->NbRows ; ++lines )
		memset( &R->p[lines][0], 0, sizeof(int)*dim );

	B = CalcBase( R );
	Matrix_Free( R );

	return( B );
}*/


/* renvoie la ligne sur laquelle on a trouve un coef non nul */
/* pareil mais cherche dans toutes les lignes */
/* et -1 s'il n'y en a pas. */
static int existepivot( Matrix *R, int l )
{
	int j, c;

	for( j=l ; j<R->NbRows ; ++j )
		if(value_notzero_p( R->p[j][l]) )
			return( j );

	/* on ne l'a pas trouve pour l'instant... on cherche au dessus */
	/* les lignes ayant que des 0 jusqu'a la position l */
	for( j=0 ; j<l ; ++j )
	{
		for( c=0 ; c<l && value_zero_p(R->p[j][c]) ; c++ )
			;
		if( c==l && value_notzero_p(R->p[j][l]) )
			return( j );
	}

	return( -1 );
}


/* echange les lignes l1 et l2 dans la matrice R */
static void swap_line(Matrix *R, int l1, int l2)
{
	int i;
	Value  tmp;

	if( l1 != l2 )
		for(i = 0;i < R->NbColumns;i ++)
		{
			value_assign(tmp , R->p[l1][i]);
			value_assign(R->p[l1][i] , R->p[l2][i]);
			value_assign(R->p[l2][i] , tmp);
		}
}



int pgcd1( int a, int b)
{
 int r;
 if( a== 0 )
 return( abs(b) );
  if(b==0 )
  return(abs(a) );
  do
  {
    r= a % b; 
    a= b;
   b = r;
   }
   while ( r!=0 );
 return(abs(a));
}




/* Soustraire la ligne l1 de l2 */
/* On effectue l2 = (l1[piv]/pgcd)*l2 - l1 * (l2[piv]/pgcd) */
static void Soustraire_ligne(Matrix *R, int l1, int l2, int piv )
{
	int i;
	Value  a, b, p, t;
	/* l2 = a*l2 - b*l1 */

	if (value_zero_p(R->p[l2][piv] ))		/* c'est deja fait ! */
		return;

	value_init(a);
	value_init(b);
	value_init(p);
	value_init(t);

	Gcd( R->p[l1][piv], R->p[l2][piv], &p );
	value_division(a, R->p[l1][piv] , p);
	value_division(b , R->p[l2][piv] , p);

	value_set_si(R->p[l2][piv] , 0);
	value_set_si(p,0);
	for(i = piv + 1;i < R->NbColumns;i ++)
	{
		value_multiply(t,b,R->p[l1][i]);
		value_multiply(R->p[l2][i],a,R->p[l2][i]);
		value_subtract(R->p[l2][i],R->p[l2][i],t);
		Gcd(p, R->p[l2][i], &p );
	}
	/* Simplification par le pgcd de toute la ligne */
	for( i=piv+1 ; i<R->NbColumns && p!=0 ; i++ )
		value_division(R->p[l2][i],R->p[l2][i], p);

	value_clear(a);
	value_clear(b);
	value_clear(p);
	value_clear(t);
}


/*** ext_ehrhart.c ****/


void new_eadd(evalue *e1,evalue *res) {
int i, p, x, y;

evalue *ne;
Value g,m1,m2;

  value_init(g);
  value_init(m1);
  value_init(m2);
    if (value_notzero_p(e1->d) && value_notzero_p(res->d)) {
         /* Add two rational numbers*/
         value_multiply(m1,e1->x.n,res->d);
         value_multiply(m2,res->x.n,e1->d);
         value_addto(res->x.n,m1,m2);
         value_multiply(res->d,e1->d,res->d);
         Gcd(res->x.n,res->d,&g);
         if (value_notone_p(g)) {
              value_division(res->d,res->d,g);
              value_division(res->x.n,res->x.n,g);
         }
         value_clear(g); value_clear(m1); value_clear(m2);
         return ;
     }
     else if (value_notzero_p(e1->d) && value_zero_p(res->d)) {
              if (res->x.p->type==polynomial) {
                  /* Add the constant to the constant term of a polynomial*/
                   new_eadd(e1, &res->x.p->arr[0]);
                   value_clear(g); value_clear(m1); value_clear(m2);
                   return ;
              }
              else if (res->x.p->type==periodic) {
                          /* Add the constant to all elements of a periodic number */
                          for (i=0; i<res->x.p->size; i++) {
                              new_eadd(e1, &res->x.p->arr[i]);
                          }
                          value_clear(g); value_clear(m1); value_clear(m2);
			     
                          return ;
                    } 
                    else {
                            fprintf(stderr, "eadd: cannot add const with vector\n");
                            value_clear(g); value_clear(m1); value_clear(m2);
                        
                            return;
                    }
     }
    /* ######### add polynomial or periodic to constant #############
       you have to exchange e1 and res, before doing addition */
     
     else if (value_zero_p(e1->d) && value_notzero_p(res->d)) {
              enode *tmp;
	      evalue x;
	      x=*res;
	      tmp= ecopy(e1->x.p);
	      value_init(res->d);
	      value_set_si( res->d, 0 );
	      res->x.p=tmp;

              new_eadd(&x,res);
	      value_clear(g); value_clear(m1); value_clear(m2);
              return ;
     }
    else {   /* ((e1->d==0) && (res->d==0)) */
                 if ((e1->x.p->type != res->x.p->type) ) {
		   /* ##### adding to evalues of different type. two cases are possible  #### 
		      
		   #### res is periodic and e1 is polynomial, you have to exchange
		   e1 and res then to add e1 to the constant term of res  #### */
		     if ((res->x.p->type == periodic)&&(e1->x.p->type == polynomial)) {
	               
		          evalue eval;
		          value_set_si( eval.d, 0 );
		          eval.x.p=ecopy(res->x.p);
	                  res->x.p= ecopy(e1->x.p);
                          new_eadd(&eval,&res->x.p->arr[0]);
			 		         	     
	             }
                     else if ((res->x.p->type == polynomial)&&(e1->x.p->type == periodic)) {
		       /* #### res is polynomial and e1 is periodic,
		          add e1 to the constant term of res  #### */
			 
			  new_eadd(e1,&res->x.p->arr[0]);
		     }
	                	 
	             value_clear(g); value_clear(m1); value_clear(m2);
		     return;
	         }
	         else if (e1->x.p->pos  != res->x.p->pos ) { 
		   /* ### adding evalues of different position (i.e function of different unknowns
		      to case are possible  ### */
			   
			 if (res->x.p->type == polynomial) {/*  ### res and e1 are polynomials
								add e1 to the constant term of res */
			       
		               new_eadd(e1,&res->x.p->arr[0]);
		               value_clear(g); value_clear(m1); value_clear(m2);
		               return;
		          }
		          else {  /* ### res and e1 are pointers to periodic numbers
				     add e1 to all elements of res */
				   
			          for (i=0;i<res->x.p->size;i++) {
			               new_eadd(e1,&res->x.p->arr[i]);
			          }
			          value_clear(g); value_clear(m1); value_clear(m2);
			          return;
		          }
                 
				          
	         }  /* ### */
                 
                
		 /* same type , same pos  and same size */
                 if (e1->x.p->size == res->x.p->size) {
		   /* add any element in e1 to the corresponding element in res */
	              for (i=0; i<res->x.p->size; i++) {
                            new_eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
                      }
                      value_clear(g); value_clear(m1); value_clear(m2);
                      return ;
                }
                
		/* Sizes are different */
                if (res->x.p->type==polynomial) {
                    /* VIN100: if e1-size > res-size you have to copy e1 in a   */
                    /* new enode and add res to that new node. If you do not do */
                    /* that, you lose the the upper weight part of e1 !         */

                     if(e1->x.p->size > res->x.p->size) {
                          enode *tmp;
                          tmp = ecopy(e1->x.p);
                          for(i=0;i<res->x.p->size;++i) {
                              new_eadd(&res->x.p->arr[i], &tmp->arr[i]);
                              /*  free_evalue_refs(&res->x.p->arr[i]); */
                          }
                          res->x.p = tmp;
    	  
                     }
                     else {
	  	     
                        for (i=0; i<e1->x.p->size ; i++) {
                             new_eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
                        } 
                        value_clear(g); value_clear(m1); value_clear(m2);
                      		   
                        return ;
                     } 
                }
                
		/* ### add two periodics of the same pos (unknown) but whith different sizes (periods) ### */
                else if (res->x.p->type==periodic) {
		      /* you have to create a new evalue 'ne' in whitch size equals to the scm
		       of the sizes of e1 and res, then to copy res periodicaly in ne, after
		       to add periodicaly elements of e1 to elements of ne, and finaly to 
		       return ne. */
		                                    
	               x=e1->x.p->size;
	               y= res->x.p->size;
	               p=ppcm1(x,y);
	               ne= (evalue *) malloc (sizeof(evalue)); 
	               value_init(ne->d);
	               value_set_si( ne->d,0);
	    	           	    
	               ne->x.p=new_enode(res->x.p->type,p, res->x.p->pos);
	               for(i=0;i<p;i++)  {
		              value_assign(ne->x.p->arr[i].d, res->x.p->arr[i%y].d);
		              if (value_notzero_p(ne->x.p->arr[i].d))   {
			    value_init(ne->x.p->arr[i].x.n);
			    value_assign(ne->x.p->arr[i].x.n, res->x.p->arr[i%y].x.n);
		              }
		              else { 
			          ne->x.p->arr[i].x.p =ecopy(res->x.p->arr[i%y].x.p);
		              }
	               }
	               for(i=0;i<p;i++)  {
	                    new_eadd(&e1->x.p->arr[i%x], &ne->x.p->arr[i]);
	               }
      
	               res=ne;
	    	    
                        value_clear(g); value_clear(m1); value_clear(m2);
                        return ;
                }
                else { /* evector */
                     fprintf(stderr, "eadd: ?cannot add vectors of different length\n");
                     value_clear(g); value_clear(m1); value_clear(m2);
                     return ;
                }
     }
     value_clear(g); value_clear(m1); value_clear(m2);
     return ;
 } /* new_eadd */
   
/* remove the last row and the last column of a matrix Mat */
Matrix *Reduce_Matrix (Matrix *Mat) {
	   int i; 
	   Value *p;
	   p=*(Mat->p+(Mat->NbRows-1));
	   for (i=0;i<Mat->NbColumns; i++)  {
		   value_clear(*p++);
	   }
	   for (i=0; i<Mat->NbRows-1; i++)  { 
		   p=*(Mat->p+i);
		   value_clear(*(p+(Mat->NbColumns-1)));
	   }
	   Mat->NbRows--;
	   Mat->NbColumns--;
	   return Mat;
   } /* Reduce_Matrix */
    
  
  /* Computes the scalar product (in euclidien space) of two vectors */ 
  
void Scalar_product(Value *p1,Value *p2,unsigned length, Value *r) {
        Value *cp1, *cp2;
	
        int i;
          cp1=p1;
       	  cp2=p2;
	  value_set_si(*r,0);
               for (i=0;i<length;i++) {
		 value_addmul(*r, *cp1, *cp2);
		  cp1++; cp2++;
	       }
  } /* Scalar_product */

    /* computes the scm of two integrals  */
     
      int ppcm1 (int a, int b) {
	   int t;
	   t = (a*b)/ pgcd1(a,b);
	  return t;
       } /* ppcm1 */
  
/* computes the scm of two values */

void ppcm(Value a, Value b, Value *r) {
	Value g;
	value_init(g); 
	Gcd(a,b,&g);
	value_multiply(*r,a,b);
	value_division(*r,*r,g);
} /* ppcm  */


Matrix *Orthogonal_Base(Matrix *Mat)  {
  Matrix *OrthMat;
  Value a,b,c,d;
  Vector *q,*p,*f;
  unsigned length;
  int i,j,k;
  value_init(a);
  value_init(b);
  value_init(c);
  value_init(d);
  OrthMat= Matrix_Alloc(Mat->NbRows,Mat->NbColumns);
  length=Mat->NbColumns;
  for(k=0; k<length; k++)  {
    value_assign(OrthMat->p[0][k],Mat->p[0][k]);		
  }
  f=Vector_Alloc(length);
  p=Vector_Alloc(length);
  q=Vector_Alloc(length);
  for(i=1; i<Mat->NbRows; i++)  {
    for(k=0;k<length;k++)  {
      value_assign(f->p[k],Mat->p[i][k]);
      value_assign(q->p[k],Mat->p[i][k]);
    }
    value_set_si(d,1);
    for(j=0; j<i; j++) {
      for(k=0;k<length;k++)  {
	value_assign(p->p[k],OrthMat->p[j][k]);
      }
      
      Scalar_product(p->p,f->p,length,&a);
      Scalar_product(p->p,p->p,length,&b);	
      Gcd(a,b,&c);
      value_division(a,a,c);
      value_division(b,b,c);
      for(k=0;k<length;k++) {
	value_multiply(p->p[k],p->p[k],a);
      }
      
      if(value_notone_p(d)|value_notone_p(b))  {
	ppcm(d,b,&c);
	value_division(a,c,b);
	value_division(b,c,d);
	value_assign(d,c);
	for(k=0;k<length;k++) {
	  value_multiply(p->p[k],p->p[k],a);
	  value_multiply(q->p[k],q->p[k],b);
	}
	
      }
      
      for(k=0;k<length;k++) {
	value_subtract(q->p[k],q->p[k],p->p[k]);
      }
      
    }
    Vector_Gcd(q->p,length,&c); 
    for(k=0;k<length;k++) {
      value_division(OrthMat->p[i][k],q->p[k],c);
    }
    
  }
  value_clear(a);
  value_clear(b);
  value_clear(c);
  value_clear(d);
  return OrthMat;
}  /* Orthogonal_Base  */
     

/* Remove an element of a list */
void  Remove_Element(Enumeration *en,Enumeration **re, Enumeration *prev)  {
  if (en== *re)  {
    *re= (*re)->next;
  }
  else {
    prev->next=en->next;
  }
} /* Remove_Element */
  
 /* Remove validite domains and correspending ehrhart polynomials whitch are redundant after 
  the enumeration of a polyhedron */

 void Remove_RedundantDomains (Enumeration  **Ures)  {
	Enumeration *ren1, *ren2, *previous=NULL;
	int red;
	for (ren1=*Ures; ren1; ren1=ren1->next)  {
		red=0;
		for (ren2=*Ures; ren2; ren2=ren2->next)  {
			if (ren1!=ren2)  {
                        	  
			    if (PolyhedronIncludes(ren2->ValidityDomain, ren1->ValidityDomain)) {
					red= 1;
					break;
		
				}
				
			}
		}
		if (red)  {
			Remove_Element(ren1,Ures,previous);
		}
		previous=ren1;
	}
 }/*Remove_RedendentDomains */
						      
				  
int IncludeInRes (Polyhedron *p, Enumeration *e, unsigned MR) {
	Enumeration *en;
	   for (en=e; en; en=en->next)  {
		  
		if (PolyhedronIncludes(e->ValidityDomain,p)) 
			 return 1;
	   }
	return 0; 	
}

Polyhedron *DMUnion(Enumeration *en, unsigned MR)  {
	Enumeration *e1;
	Polyhedron *d;
	e1=en;
        d=e1->ValidityDomain; 
	 for (e1=en->next; e1; e1=e1->next) {
	   d= DomainUnion( d, e1->ValidityDomain, MR);
	 }
	 return d;
  }


void AffConstraints(Polyhedron *Poldisj)
{
	Polyhedron *p;

	for(p=Poldisj;p;p=p->next)
	{
		Polyhedron_PrintConstraints( stdout, P_VALUE_FMT, p);
		printf("\n");
	}
}
int Degenerate (Enumeration *en) {
	if(value_notzero_p(en->EP.d)) {
		
           if(value_mone_p(en->EP.x.n )) {
	       return 1;    
           }
	}
   return 0;	
}

/* Enumeration of a domain D */

Enumeration *Domain_Enumerate(Polyhedron *D, Polyhedron *C, unsigned MAXRAYS,char **pn)
{     Polyhedron_union *Polun,*pu;
      Polyhedron  *lp, *lp1, *lp1next;
      Polyhedron *d1,*d2,*d;  
      Enumeration *e,*pr,*en,*en1, *en2,*tmp, *res, *sen;
      Polun=NULL;
     
      for (d = D; d; d = d->next) {
	  POL_ENSURE_FACETS(d);
	  POL_ENSURE_VERTICES(d);
      }
      POL_ENSURE_FACETS(C);
      POL_ENSURE_VERTICES(C);

     lp = Disjoint_Domain( D, 0, MAXRAYS );

#ifdef UE_DEBUG
     printf("##############################################################\n");
     printf("\n###### DISJOINT UNION  ######\n\n");
     AffConstraints(lp); 
     printf("##############################################################\n");
#endif
     
	for (lp1=lp ; lp1; lp1=lp1->next)
	{
		Enumeration *enext;
		lp1next = lp1->next;
		lp1->next = NULL;
		en= Polyhedron_Enumerate(lp1, C, MAXRAYS,NULL);
		lp1->next = lp1next;
		sen= NULL;
		for (e=en;e;e=enext) {
			enext = e->next;
			if (!Degenerate(e)) {
				e->next = sen;
				sen=e;
			} else {
				free_evalue_refs(&e->EP);
				Domain_Free(e->ValidityDomain);
				free(e);
			}
		}

		if(sen!= NULL)
		{
			pu = (Polyhedron_union  *)malloc(sizeof(Polyhedron_union));
			pu->pt=sen;
			pu->next = Polun;
			Polun = pu;
		}
	}
	if(!Polun)
	{
#ifdef UE_DEBUG
		fprintf(stdout,"Empty Polun\n");
#endif
		return ((Enumeration *) 0);
	}
      
	while(Polun->next != NULL)  {
		Enumeration *enext;
		res=NULL;
		en1=Polun->pt;
		en2=(Polun->next)->pt;

		d1=DMUnion(en1, MAXRAYS);
		d2=DMUnion(en2, MAXRAYS);

		for (en1=Polun->pt;en1;en1=enext) {
			enext = en1->next;
			for(en2=(Polun->next)->pt;en2;en2=en2->next)
			{
				d = DomainIntersection(en1->ValidityDomain,en2->ValidityDomain,MAXRAYS);
				if( d && !emptyQ(d)&&!IncludeInRes(d,res,MAXRAYS))  {
					tmp = (Enumeration  *)malloc(sizeof(Enumeration));
					value_init(tmp->EP.d);
					value_assign( tmp->EP.d, en2->EP.d );
					if(value_zero_p(tmp->EP.d))
						tmp->EP.x.p=ecopy(en2->EP.x.p);
					else
					{
						value_init(tmp->EP.x.n);
						value_assign( tmp->EP.x.n, en2->EP.x.n );
					}

					new_eadd(&en1->EP,&tmp->EP);
					tmp->ValidityDomain =d;
					tmp->next= res;
					res=tmp;
				}
			}
			d=DomainDifference(en1->ValidityDomain,d2 ,MAXRAYS);
			if (d && !emptyQ(d) && !IncludeInRes(d,res,MAXRAYS)) {
				en1->ValidityDomain = d;
				en1->next= res;
				res=en1;
			} else {
				free_evalue_refs(&en1->EP);
				free(en1);
			}
		}
		for (en2=(Polun->next)->pt; en2; en2 = enext) {
			enext = en2->next;
			d= DomainDifference(en2->ValidityDomain,d1,MAXRAYS);
			if (d && !emptyQ(d)&&!IncludeInRes(d,res,MAXRAYS)) {
				en2->ValidityDomain = d;
				en2->next = res;
				res = en2;
			} else {
				free_evalue_refs(&en2->EP);
				free(en2);
			}
		}
		Domain_Free(d1);
		Domain_Free(d2);
	    
		Polun->pt=res;
	        		     
		Polun->next= (Polun->next)->next;
	}
	res=Polun->pt;
		
	Remove_RedundantDomains(&res); 
	return(res);
}


/**********
 DO NOT USE THE FOLLOWING FUNCTION IT'S NOT WORKING PROPERLY YET.
 **********/

/* Enumeration of the image by T of domain D */
Enumeration *Polyhedron_Image_Enumerate(Polyhedron *D,  Polyhedron *C, Matrix *T, unsigned MAXRAYS, char **par_name)
{   Polyhedron *polun,*pol;
    Enumeration *ee;
    Matrix *TCopy,*Tred, *d1,*d;
    Vector *v1,*v2;
    Value h;
    int i,j,k;

  POL_ENSURE_FACETS(D);
  POL_ENSURE_VERTICES(D);
  POL_ENSURE_FACETS(C);
  POL_ENSURE_VERTICES(C);

   value_init(h);
    if(!D) {
	 fprintf(stdout,"             Error: in reading input domain \n");   
	 value_clear(h);
	   return ((Enumeration *) 0);
    }
    else {
      printf("\n ################   INPUT  POLYHEDRON  #######################\n\n");
      AffConstraints(D);
     }

#ifdef DOMAIN_IMAGE
       fpol=DomainImage(D,T,MAXRAYS);
       printf("\n $$$$$$$$$$$$$  THE  DOMAIN IMAGE    $$$$$$$$$$$$$\n\n");
         AffConstraints(fpol);
          if(emptyQ(fpol)) {
		  value_clear(h);
 		  return ((Enumeration *) 0);
	  } 
          ee = Domain_Enumerate(fpol,C,MAXRAYS,par_name);
	  value_clear(h);
       return  (ee);
#endif
 
     TCopy= Matrix_Copy(T);
     Tred= Reduce_Matrix(TCopy);
     printf("\n ##################  INPUT REDUCED TRANSFORMATION MATRIX ##################\n" );
     Matrix_Print(stdout,P_VALUE_FMT,Tred);

         if (Tred->NbRows <Tred->NbColumns) {
               d1=(Matrix *) Matrix_Alloc(Tred->NbColumns,Tred->NbColumns);
               for (i=0;i<Tred->NbRows;i++) {
	 	    for (j=0; j<Tred->NbColumns;j++) {
			 	value_assign( d1->p[i][j], Tred->p[i][j] );
		    }
	       }
	       for(i=Tred->NbRows;i<Tred->NbColumns;i++) {
		     for (j=0;j<Tred->NbColumns;j++) {
		          value_set_si( d1->p[i][j], 0 );
	             }
	       }
	       d= (Matrix *) CalcBase(d1);
	       Matrix_Free(Tred);
	       Matrix_Free(d1);
	 
         }
         else {
              d=(Matrix *) CalcBase(Tred);
              Matrix_Free(Tred);
         }
    if(d->NbRows==0) {
	if(emptyQ(D)) {
		value_clear(h);
		return ((Enumeration *) 0);
	}
	else {
	   printf( "        Ker(A)=0  implys directly Enumeration on input polyhedron\n\n");    	
           ee=Domain_Enumerate(D,C,MAXRAYS,par_name);
	   value_clear(h);
           return ee;
	}
   }
   
   d1=Transpose(d);
   Matrix_Free(d);
   
   if(d1->NbRows!=D->Dimension) {
      fprintf(stdout,"      \n Error: incompatible dimension \n");
      value_clear(h);
      return ((Enumeration *) 0);
   }
   if(d1->NbColumns > 1) {
 fprintf(stdout,"   \n Error: Can not compute intégral points : More then vector in ker(A)! \n");
    value_clear(h);
    return ((Enumeration *) 0);
	 
    }
   printf( "           \n Ker(A)=1  implys adding constraints befor Enumeration\n");
   v1=Vector_Alloc(d1->NbRows);
   v2=Vector_Alloc(d1->NbRows);
   
   	  polun=(Polyhedron *) NULL; 
            for (k=0;k<d1->NbRows;k++) {
	          value_assign(v1->p[k],d1->p[k][0]) ;
	    }
	  /* adding a new constraint for all constraints of D in which the scalar product of the*/
	  /* normal whith vector v1 is greter then zero*/
	    
	  for (j=0;j<D->NbConstraints;j++)  {
                 for (k=0;k<=D->Dimension-1;k++) {
			value_assign(v2->p[k],D->Constraint[j][k+1]) ;
	          }
                  Scalar_product(v1->p,v2->p,D->Dimension,&h);
		
		  if(value_pos_p(h)&&!value_zero_p(D->Constraint[j][0])) {
	               Vector *NCont;
		       Value val;
		       value_init( val );
		       /* Create a new contraint whitch is added to the polyhedron*/
		       
		       NCont=Vector_Alloc(d1->NbRows+2);
		       value_set_si( NCont->p[0],1); /* the constraint is an inequality */
				       
		       for (k=1;k<=D->Dimension;k++) {
		            value_oppose( NCont->p[k], D->Constraint[j][k]);
					}
		       value_decrement(val,h);
		       value_subtract(val,val,D->Constraint[j][D->Dimension+1]);
		       value_assign (NCont->p[D->Dimension+1],val);
		       value_clear(val);
		       /* add the new constraint to polyhedron D */
		       pol=AddConstraints(NCont->p,1,D,MAXRAYS);
		       polun=AddPolyToDomain(Polyhedron_Copy(pol),polun);
		       Polyhedron_Free(pol);
		       Vector_Free(NCont);
				 value_clear( val );
		    }
	   }
	  if(polun==NULL) { /*  No constraint is added to input polyhedron */
	      if(emptyQ(D)) {
		      value_clear(h);
		      return ((Enumeration *) 0);
	      }
	      else {
	         ee= Domain_Enumerate(D,C,MAXRAYS,par_name);
	      }
	  }
	  else { /* some constraintes are added to input polyhedron */
	      if(emptyQ(polun)){
		      value_clear(h);
         	      return ((Enumeration *) 0);
              }
              else {
	printf("\n ##################################################################");       
        printf("\n ****** THE RESULT OF ADDING CONSTRAINTS TO THE INPUT POLYHEDRON  ****** \n");  	 
	       AffConstraints(polun);
       	       ee= Domain_Enumerate(polun,C,MAXRAYS,par_name);
	       value_clear(h);
	       return (ee );
	      }
      }
	      

	return( NULL );	          
}


