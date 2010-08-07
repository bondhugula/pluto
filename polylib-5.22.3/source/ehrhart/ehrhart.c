/***********************************************************************/
/*                Ehrhart V4.20                                        */
/*                copyright 1997, Doran Wilde                          */
/*                copyright 1997-2000, Vincent Loechner                */
/*       Permission is granted to copy, use, and distribute            */
/*       for any commercial or noncommercial purpose under the terms   */
/*       of the GNU General Public license, version 2, June 1991       */
/*       (see file : LICENSING).                                       */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include <polylib/polylib.h>
#include <polylib/homogenization.h>


/*! \class Ehrhart

The following are mainly for debug purposes. You shouldn't need to
change anything for daily usage...

<p>

you may define each macro independently 
<ol>
<li> #define EDEBUG minimal debug 
<li> #define EDEBUG1 prints enumeration points
<li> #define EDEBUG11 prints number of points
<li> #define EDEBUG2 prints domains
<li> #define EDEBUG21 prints more domains
<li> #define EDEBUG3 prints systems of equations that are solved
<li> #define EDEBUG4 prints message for degree reduction
<li> #define EDEBUG5 prints result before simplification 
<li> #define EDEBUG6 prints domains in Preprocess 
<li> #define EDEBUG61 prints even more in Preprocess
<li> #define EDEBUG62 prints domains in Preprocess2
</ol>
*/


/** 
    
define this to print all constraints on the validity domains if not
defined, only new constraints (not in validity domain given by the
user) are printed

*/
#define EPRINT_ALL_VALIDITY_CONSTRAINTS

/* #define EDEBUG 	*/		/* minimal debug */
/* #define EDEBUG1	*/		/* prints enumeration points */
/* #define EDEBUG11	*/		/* prints number of points */
/* #define EDEBUG2	*/		/* prints domains */
/* #define EDEBUG21	*/		/* prints more domains */
/* #define EDEBUG3	*/		/* prints systems of equations that are solved */
/* #define EDEBUG4	*/		/* prints message for degree reduction */
/* #define EDEBUG5	*/		/* prints result before simplification */
/* #define EDEBUG6	*/		/* prints domains in Preprocess */
/* #define EDEBUG61	*/		/* prints even more in Preprocess */
/* #define EDEBUG62	*/		/* prints domains in Preprocess2 */


/**
 Reduce the degree of resulting polynomials
*/
#define REDUCE_DEGREE

/** 
define this to print one warning message per domain overflow these
overflows should no longer happen since version 4.20
*/
#define ALL_OVERFLOW_WARNINGS

/******************* -----------END USER #DEFS-------- *********************/

int overflow_warning_flag = 1;

/*-------------------------------------------------------------------*/
/* EHRHART POLYNOMIAL SYMBOLIC ALGEBRA SYSTEM                        */
/*-------------------------------------------------------------------*/ 
/** 

EHRHART POLYNOMIAL SYMBOLIC ALGEBRA SYSTEM. The newly allocated enode
can be freed with a simple free(x)

@param type : enode type
@param size : degree+1 for polynomial, period for periodic
@param pos  : 1..nb_param, position of parameter            
@return a newly allocated enode 

*/
enode *new_enode(enode_type type,int size,int pos) {
  
  enode *res;
  int i;
  
  if(size == 0) {
    fprintf(stderr, "Allocating enode of size 0 !\n" );
    return NULL;
  }
  res = (enode *) malloc(sizeof(enode) + (size-1)*sizeof(evalue));
  res->type = type;
  res->size = size;
  res->pos = pos;
  for(i=0; i<size; i++) {
    value_init(res->arr[i].d);
    value_set_si(res->arr[i].d,0);
    res->arr[i].x.p = 0;
  }
  return res;
} /* new_enode */

/**
releases all memory referenced by e.  (recursive)        
@param e pointer to an evalue
*/
void free_evalue_refs(evalue *e) {
  
  enode *p;
  int i;
  
  if (value_notzero_p(e->d)) {
    
    /* 'e' stores a constant */
    value_clear(e->d);
    value_clear(e->x.n);
    return; 
  }  
  value_clear(e->d);
  p = e->x.p;
  if (!p) return;	/* null pointer */
  for (i=0; i<p->size; i++) {
    free_evalue_refs(&(p->arr[i]));
  }
  free(p);
  return;
} /* free_evalue_refs */

/**

@param e pointer to an evalue 
@return description

*/
enode *ecopy(enode *e) {
  
  enode *res;
  int i;
  
  res = new_enode(e->type,e->size,e->pos);
  for(i=0;i<e->size;++i) {
    value_assign(res->arr[i].d,e->arr[i].d);
    if(value_zero_p(res->arr[i].d))
      res->arr[i].x.p = ecopy(e->arr[i].x.p);
    else {
      value_init(res->arr[i].x.n);
      value_assign(res->arr[i].x.n,e->arr[i].x.n);
    }
  }
  return(res);
} /* ecopy */

/**

@param DST destination file
@param e pointer to evalue to be printed
@param pname array of strings, name of the parameters

*/
void print_evalue(FILE *DST,evalue *e,char **pname) {
  
  if(value_notzero_p(e->d)) {    
    if(value_notone_p(e->d)) {
      value_print(DST,VALUE_FMT,e->x.n);
      fprintf(DST,"/");
      value_print(DST,VALUE_FMT,e->d);
    }  
    else {
      value_print(DST,VALUE_FMT,e->x.n);
    }
  }  
  else
    print_enode(DST,e->x.p,pname);
  return;
} /* print_evalue */

/** prints the enode  to DST  

@param DST destination file 
@param p pointer to enode  to be printed 
@param pname array of strings, name of the parameters

*/
void print_enode(FILE *DST,enode *p,char **pname) {
  
  int i;
  
  if (!p) {
    fprintf(DST, "NULL");
    return;
  }
  if (p->type == evector) {
    fprintf(DST, "{ ");
    for (i=0; i<p->size; i++) {
      print_evalue(DST, &p->arr[i], pname);
      if (i!=(p->size-1))
	fprintf(DST, ", ");
    }
    fprintf(DST, " }\n");
  }
  else if (p->type == polynomial) {
    fprintf(DST, "( ");
    for (i=p->size-1; i>=0; i--) {
      print_evalue(DST, &p->arr[i], pname);
      if (i==1) fprintf(DST, " * %s + ", pname[p->pos-1]);
      else if (i>1) 
	fprintf(DST, " * %s^%d + ", pname[p->pos-1], i);
    }
    fprintf(DST, " )\n");
  }
  else if (p->type == periodic) {
    fprintf(DST, "[ ");
    for (i=0; i<p->size; i++) {
      print_evalue(DST, &p->arr[i], pname);
      if (i!=(p->size-1)) fprintf(DST, ", ");
    }
    fprintf(DST," ]_%s", pname[p->pos-1]);
  }
  return;
} /* print_enode */ 

/**

@param  e1 pointers to evalues
@param  e2 pointers to evalues
@return 1 (true) if they are equal, 0 (false) if not       

*/
static int eequal(evalue *e1,evalue *e2) { 
 
    int i;
    enode *p1, *p2;
  
    if (value_ne(e1->d,e2->d))
        return 0;
  
    /* e1->d == e2->d */
    if (value_notzero_p(e1->d)) {    
        if (value_ne(e1->x.n,e2->x.n))
            return 0;
    
        /* e1->d == e2->d != 0  AND e1->n == e2->n */
        return 1;
    }
  
    /* e1->d == e2->d == 0 */
    p1 = e1->x.p;
    p2 = e2->x.p;
    if (p1->type != p2->type) return 0;
    if (p1->size != p2->size) return 0;
    if (p1->pos  != p2->pos) return 0;
    for (i=0; i<p1->size; i++)
        if (!eequal(&p1->arr[i], &p2->arr[i]) ) 
            return 0;
    return 1;
} /* eequal */

/** 

@param e pointer to an evalue

*/
void reduce_evalue (evalue *e) {
  
    enode *p;
    int i, j, k;
  
    if (value_notzero_p(e->d))
        return;	/* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	/* hum... an overflow probably occured */
  
    /* First reduce the components of p */
    for (i=0; i<p->size; i++)
        reduce_evalue(&p->arr[i]);

    if (p->type==periodic) {
    
        /* Try to reduce the period */
        for (i=1; i<=(p->size)/2; i++) {
            if ((p->size % i)==0) {
	
                /* Can we reduce the size to i ? */
                for (j=0; j<i; j++)
                    for (k=j+i; k<e->x.p->size; k+=i)
                        if (!eequal(&p->arr[j], &p->arr[k])) goto you_lose;

                /* OK, lets do it */
                for (j=i; j<p->size; j++) free_evalue_refs(&p->arr[j]);
                p->size = i;
                break;

            you_lose:   /* OK, lets not do it */
                continue;
            }
        }

        /* Try to reduce its strength */
        if (p->size == 1) {
	    value_clear(e->d);
            memcpy(e,&p->arr[0],sizeof(evalue));
            free(p);
        }
    }
    else if (p->type==polynomial) {
	  
        /* Try to reduce the degree */
        for (i=p->size-1;i>=1;i--) {
            if (!(value_one_p(p->arr[i].d) && value_zero_p(p->arr[i].x.n)))
                break;
	    /* Zero coefficient */
	    free_evalue_refs(&p->arr[i]);
        }
        if (i+1<p->size)
	    p->size = i+1;

        /* Try to reduce its strength */
        if (p->size == 1) {
	    value_clear(e->d);
            memcpy(e,&p->arr[0],sizeof(evalue));
            free(p);
        }
    }
} /* reduce_evalue */


/**

multiplies two evalues and puts the result in res

@param e1 pointer to an evalue
@param e2 pointer to a constant evalue
@param res  pointer to result evalue = e1 * e2

*/
static void emul (evalue *e1,evalue *e2,evalue *res) {
  
    enode *p;
    int i;
    Value g;

    if (value_zero_p(e2->d)) {
        fprintf(stderr, "emul: ?expecting constant value\n");
        return;
    }  
    value_init(g);
    if (value_notzero_p(e1->d)) {
    
	value_init(res->x.n);
        /* Product of two rational numbers */
        value_multiply(res->d,e1->d,e2->d);
        value_multiply(res->x.n,e1->x.n,e2->x.n );
        Gcd(res->x.n, res->d,&g);
        if (value_notone_p(g)) {
            value_division(res->d,res->d,g);
            value_division(res->x.n,res->x.n,g);
        }
    }
    else { /* e1 is an expression */    
        value_set_si(res->d,0);
        p = e1->x.p;
        res->x.p = new_enode(p->type, p->size, p->pos);
        for (i=0; i<p->size; i++) {
            emul(&p->arr[i], e2, &(res->x.p->arr[i]) ); 
        }
    }
    value_clear(g);
    return;
} /* emul */

/**
adds one evalue to evalue 'res. result = res + e1 

@param e1 an evalue 
@param res 

*/
void eadd(evalue *e1,evalue *res) {

  int i;
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
    return;
  }
  else if (value_notzero_p(e1->d) && value_zero_p(res->d)) {
    if (res->x.p->type==polynomial) {
      
      /* Add the constant to the constant term */
      eadd(e1, &res->x.p->arr[0]);
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
    }
    else if (res->x.p->type==periodic) {
      
      /* Add the constant to all elements of periodic number */
      for (i=0; i<res->x.p->size; i++) { 
	eadd(e1, &res->x.p->arr[i]);
      }
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
    }
    else {
      fprintf(stderr, "eadd: cannot add const with vector\n");
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
    }
  }
  else if (value_zero_p(e1->d) && value_notzero_p(res->d)) {
    fprintf(stderr,"eadd: cannot add evalue to const\n");
    value_clear(g); value_clear(m1); value_clear(m2);
    return;
  }
  else {    /* ((e1->d==0) && (res->d==0)) */  
    if ((e1->x.p->type != res->x.p->type) ||
	(e1->x.p->pos  != res->x.p->pos )) {
      fprintf(stderr, "eadd: ?cannot add, incompatible types\n");
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
    }
    if (e1->x.p->size == res->x.p->size) {
      for (i=0; i<res->x.p->size; i++) {
	eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
      }
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
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
	  eadd(&res->x.p->arr[i], &tmp->arr[i]);
	  free_evalue_refs(&res->x.p->arr[i]);
	}
	res->x.p = tmp;
      }
      else {
	for (i=0; i<e1->x.p->size ; i++) {
	  eadd(&e1->x.p->arr[i], &res->x.p->arr[i]);
	}
	value_clear(g); value_clear(m1); value_clear(m2);
	return;
      }
    }
    else if (res->x.p->type==periodic) {
      fprintf(stderr, "eadd: ?addition of different sized periodic nos\n");
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
    }
    else { /* evector */
      fprintf(stderr, "eadd: ?cannot add vectors of different length\n");
      value_clear(g); value_clear(m1); value_clear(m2);
      return;
    }
  }
  value_clear(g); value_clear(m1); value_clear(m2);
  return;
} /* eadd */

/** 

computes the inner product of two vectors. Result = result (evalue) =
v1.v2 (dot product)

@param  v1 an enode (vector)
@param  v2 an enode (vector of constants)
@param res result (evalue)

*/
void edot(enode *v1,enode *v2,evalue *res) {
  
    int i;
    evalue tmp;
  
    if ((v1->type != evector) || (v2->type != evector)) {
        fprintf(stderr, "edot: ?expecting vectors\n");
        return;
    }
    if (v1->size != v2->size) {
        fprintf(stderr, "edot: ? vector lengths do not agree\n");
        return;
    }
    if (v1->size<=0) {
        value_set_si(res->d,1);	/* set result to 0/1 */
	value_init(res->x.n);
        value_set_si(res->x.n,0);
        return;
    }
  
    /* vector v2 is expected to have only rational numbers in */
    /* the array.  No pointers. */  
    emul(&v1->arr[0],&v2->arr[0],res);
    for (i=1; i<v1->size; i++) {
        value_init(tmp.d);
     
        /* res = res + v1[i]*v2[i] */
        emul(&v1->arr[i],&v2->arr[i],&tmp);
        eadd(&tmp,res);
        free_evalue_refs(&tmp);
    }
    return;
} /* edot */

/**

local recursive function used in the following ref contains the new
position for each old index position

@param e pointer to an evalue
@param ref transformation Matrix

*/
static void aep_evalue(evalue *e, int *ref) {
  
    enode *p;
    int i;
  
    if (value_notzero_p(e->d))
        return;	        /* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	        /* hum... an overflow probably occured */
  
    /* First check the components of p */
    for (i=0;i<p->size;i++)
        aep_evalue(&p->arr[i],ref);
  
    /* Then p itself */
    p->pos = ref[p->pos-1]+1;
    return;
} /* aep_evalue */

/** Comments */
static void addeliminatedparams_evalue(evalue *e,Matrix *CT) {
	
    enode *p;
    int i, j;
    int *ref;

    if (value_notzero_p(e->d))
        return;	         /* a rational number, its already reduced */
    if(!(p = e->x.p))
        return;	         /* hum... an overflow probably occured */
  
    /* Compute ref */
    ref = (int *)malloc(sizeof(int)*(CT->NbRows-1));
    for(i=0;i<CT->NbRows-1;i++)
        for(j=0;j<CT->NbColumns;j++)
            if(value_notzero_p(CT->p[i][j])) {
                ref[i] = j;
                break;
            }
  
    /* Transform the references in e, using ref */
    aep_evalue(e,ref);
    free( ref );
    return;
} /* addeliminatedparams_evalue */

/** 

This procedure finds an integer point contained in polyhedron D /
first checks for positive values, then for negative values returns
TRUE on success. Result is in min.  returns FALSE if no integer point
is found

<p>

This is the maximum number of iterations for a given parameter to find
a integer point inside the context. Kind of weird. cherche_min should

<p>

@param min
@param D
@param pos

*/
/* FIXME: needs to be rewritten ! */
#define MAXITER 100
int cherche_min(Value *min,Polyhedron *D,int pos) {
	
    Value binf, bsup;	/* upper & lower bound */
    Value i;
    int flag, maxiter;
  
    if(!D)
        return(1);
    if(pos > D->Dimension)
        return(1);
  
    value_init(binf); value_init(bsup);
    value_init(i);
  
#ifdef EDEBUG61
    fprintf(stderr,"Entering Cherche min --> \n");
    fprintf(stderr,"LowerUpperBounds :\n");
    fprintf(stderr,"pos = %d\n",pos);
    fprintf(stderr,"current min = (");
    value_print(stderr,P_VALUE_FMT,min[1]);
    {int j;
    for(j=2;j<=D->Dimension ; j++) {
        fprintf(stderr,", ");
        value_print(stderr,P_VALUE_FMT,min[j]);
    }  
    }
    fprintf(stderr,")\n");
#endif
	
    flag = lower_upper_bounds(pos,D,min,&binf,&bsup);
  
#ifdef EDEBUG61
    fprintf(stderr, "flag = %d\n", flag);
    fprintf(stderr,"binf = ");
    value_print(stderr,P_VALUE_FMT,binf);
    fprintf(stderr,"\n");
    fprintf(stderr,"bsup = ");
    value_print(stderr,P_VALUE_FMT,bsup);
    fprintf(stderr,"\n");
#endif

    if(flag&LB_INFINITY)
        value_set_si(binf,0);
  
    /* Loop from 0 (or binf if positive) to bsup */
    for(maxiter=0,(((flag&LB_INFINITY) || value_neg_p(binf)) ? 
            value_set_si(i,0) : value_assign(i,binf));
        ((flag&UB_INFINITY) || value_le(i,bsup)) && maxiter<MAXITER ;
        value_increment(i,i),maxiter++) {
    
        value_assign(min[pos],i);
        if(cherche_min(min,D->next,pos+1)) {
            value_clear(binf); value_clear(bsup);
            value_clear(i);
            return(1);
        }  
    }
  
    /* Descending loop from -1 (or bsup if negative) to binf */
    if((flag&LB_INFINITY) || value_neg_p(binf))
        for(maxiter=0,(((flag&UB_INFINITY) || value_pos_p(bsup))?
                value_set_si(i,-1)
                :value_assign(i,bsup));
            ((flag&LB_INFINITY) || value_ge(i,binf)) && maxiter<MAXITER  ;
            value_decrement(i,i),maxiter++) {

            value_assign(min[pos],i);
            if(cherche_min(min,D->next,pos+1)) {
                value_clear(binf); value_clear(bsup);
                value_clear(i);
                return(1);
            }	
        }
    value_clear(binf); value_clear(bsup);
    value_clear(i);

    value_set_si(min[pos],0);
    return(0);	     /* not found :-( */
} /* cherche_min */

/** 

This procedure finds the smallest parallelepiped of size
'<i>size[i]</i>' for every dimension i, contained in polyhedron D.
If this is not possible, NULL is returned

<p>

<pre>
written by vin100, 2000, for version 4.19
modified 2002, version 5.10
</pre>

<p>

It first finds the coordinates of the lexicographically smallest edge
of the hypercube, obtained by transforming the constraints of D (by
adding 'size' as many times as there are negative coeficients in each
constraint), and finding the lexicographical min of this
polyhedron. Then it builds the hypercube and returns it.

<p>
@param D
@param size
@param MAXRAYS

*/

Polyhedron *Polyhedron_Preprocess(Polyhedron *D,Value *size,unsigned MAXRAYS)
{
    Matrix *M;
    int i, j, d;
    Polyhedron *T, *S, *H, *C;
    Value *min;

    d = D->Dimension;
    if (MAXRAYS < 2*D->NbConstraints)
	MAXRAYS = 2*D->NbConstraints;
    M = Matrix_Alloc(MAXRAYS, D->Dimension+2);
    M->NbRows = D->NbConstraints;
  
    /* Original constraints */
    for(i=0;i<D->NbConstraints;i++)
        Vector_Copy(D->Constraint[i],M->p[i],(d+2));

#ifdef EDEBUG6
    fprintf(stderr,"M for PreProcess : ");
    Matrix_Print(stderr,P_VALUE_FMT,M);
    fprintf(stderr,"\nsize == ");
    for( i=0 ; i<d ; i++ )
        value_print(stderr,P_VALUE_FMT,size[i]);
    fprintf(stderr,"\n");
#endif

    /* Additionnal constraints */
    for(i=0;i<D->NbConstraints;i++) {
        if(value_zero_p(D->Constraint[i][0])) {
            fprintf(stderr,"Polyhedron_Preprocess: ");
            fprintf(stderr,
		    "an equality was found where I did expect an inequality.\n");
            fprintf(stderr,"Trying to continue...\n");
            continue;
        }
        Vector_Copy(D->Constraint[i],M->p[M->NbRows],(d+2));
        for(j=1;j<=d;j++)
            if(value_neg_p(D->Constraint[i][j])) {
	      value_addmul(M->p[M->NbRows][d+1], D->Constraint[i][j], size[j-1]);
            }
    
        /* If anything changed, add this new constraint */
        if(value_ne(M->p[M->NbRows][d+1],D->Constraint[i][d+1]))
            M->NbRows ++ ;
    }
  
#ifdef EDEBUG6
    fprintf(stderr,"M used to find min : ");
    Matrix_Print(stderr,P_VALUE_FMT,M);
#endif
  
    T = Constraints2Polyhedron(M,MAXRAYS);
    Matrix_Free(M);
    if (!T || emptyQ(T)) {
        if(T)
            Polyhedron_Free(T);
        return(NULL);
    }
  
    /* Ok, now find the lexicographical min of T */
    min = (Value *) malloc(sizeof(Value) * (d+2));
    for(i=0;i<=d;i++) {
        value_init(min[i]);
        value_set_si(min[i],0);
    }  
    value_init(min[i]);
    value_set_si(min[i],1);
    C = Universe_Polyhedron(0);
    S = Polyhedron_Scan(T,C,MAXRAYS);
    Polyhedron_Free(C);
    Polyhedron_Free(T);

#ifdef EDEBUG6
    for(i=0;i<=(d+1);i++) {
        value_print(stderr,P_VALUE_FMT,min[i]);
        fprintf(stderr," ,");
    }
    fprintf(stderr,"\n");
    Polyhedron_Print(stderr,P_VALUE_FMT,S);
    fprintf(stderr,"\n");
#endif
    
    if (!cherche_min(min,S,1))
        {
            for(i=0;i<=(d+1);i++)
                value_clear(min[i]);
            return(NULL);
        }
    Domain_Free(S);
  
#ifdef EDEBUG6
    fprintf(stderr,"min = ( ");
    value_print(stderr,P_VALUE_FMT,min[1]);
    for(i=2;i<=d;i++) {
        fprintf(stderr,", ");
        value_print(stderr,P_VALUE_FMT,min[i]);
    }  
    fprintf(stderr,")\n");
#endif
  
    /* Min is the point from which we can construct the hypercube */
    M = Matrix_Alloc(d*2,d+2);
    for(i=0;i<d;i++) {
    
        /* Creates inequality  1 0..0 1 0..0 -min[i+1] */
        value_set_si(M->p[2*i][0],1);
        for(j=1;j<=d;j++)
            value_set_si(M->p[2*i][j],0);
        value_set_si(M->p[2*i][i+1],1);
        value_oppose(M->p[2*i][d+1],min[i+1]);
    
        /* Creates inequality  1 0..0 -1 0..0 min[i+1]+size -1 */
        value_set_si(M->p[2*i+1][0],1);
        for(j=1;j<=d;j++)
            value_set_si(M->p[2*i+1][j],0);
        value_set_si(M->p[2*i+1][i+1],-1);
        value_addto(M->p[2*i+1][d+1],min[i+1],size[i]);
        value_sub_int(M->p[2*i+1][d+1],M->p[2*i+1][d+1],1);
    }
  
#ifdef EDEBUG6
    fprintf(stderr,"PolyhedronPreprocess: constraints H = ");
    Matrix_Print(stderr,P_VALUE_FMT,M);
#endif
  
    H = Constraints2Polyhedron(M,MAXRAYS);
  
#ifdef EDEBUG6
    Polyhedron_Print(stderr,P_VALUE_FMT,H);
    fprintf(stderr,"\n");
#endif
  
    Matrix_Free(M);
    for(i=0;i<=(d+1);i++)
        value_clear(min[i]);
    free(min);
    assert(!emptyQ(H));
    return(H);
} /* Polyhedron_Preprocess */


/** This procedure finds an hypercube of size 'size', containing
polyhedron D increases size and lcm if necessary (and not "too big")
If this is not possible, an empty polyhedron is returned

<p>

<pre> written by vin100, 2001, for version 4.19</pre>

@param D
@param size
@param lcm
@param MAXRAYS

*/
Polyhedron *Polyhedron_Preprocess2(Polyhedron *D,Value *size,
				   Value *lcm,unsigned MAXRAYS) {
  
    Matrix *c;
    Polyhedron *H;
    int i,j,r;
    Value n;      /* smallest/biggest value */
    Value s;	/* size in this dimension */
    Value tmp1,tmp2;

#ifdef EDEBUG62
    int np;
#endif

    value_init(n); value_init(s);
    value_init(tmp1); value_init(tmp2);
    c = Matrix_Alloc(D->Dimension*2,D->Dimension+2);
  
#ifdef EDEBUG62
    fprintf(stderr,"\nPreProcess2 : starting\n");
    fprintf(stderr,"lcm = ");
    for( np=0 ; np<D->Dimension; np++ )
        value_print(stderr,VALUE_FMT,lcm[np]);
    fprintf(stderr,", size = ");
    for( np=0 ; np<D->Dimension; np++ )
        value_print(stderr,VALUE_FMT,size[np]);
    fprintf(stderr,"\n");
#endif
  
    for(i=0;i<D->Dimension;i++) {
    
        /* Create constraint 1 0..0 1 0..0 -min */
        value_set_si(c->p[2*i][0],1);
        for(j=0;j<D->Dimension;j++)
            value_set_si(c->p[2*i][1+j],0);
        value_division(n,D->Ray[0][i+1],D->Ray[0][D->Dimension+1]);
        for(r=1;r<D->NbRays;r++) {
            value_division(tmp1,D->Ray[r][i+1],D->Ray[r][D->Dimension+1]);
            if(value_gt(n,tmp1)) {
	
                /* New min */
                value_division(n,D->Ray[r][i+1],D->Ray[r][D->Dimension+1]);
            }
        }
        value_set_si(c->p[2*i][i+1],1);
        value_oppose(c->p[2*i][D->Dimension+1],n);
    
        /* Create constraint 1 0..0 -1 0..0 max */
        value_set_si(c->p[2*i+1][0],1);
        for(j=0;j<D->Dimension;j++)
            value_set_si(c->p[2*i+1][1+j],0);
    
        /* n = (num+den-1)/den */
        value_addto(tmp1,D->Ray[0][i+1],D->Ray[0][D->Dimension+1]);
        value_sub_int(tmp1,tmp1,1);
        value_division(n,tmp1,D->Ray[0][D->Dimension+1]);
        for(r=1;r<D->NbRays;r++) {
            value_addto(tmp1,D->Ray[r][i+1],D->Ray[r][D->Dimension+1]);
            value_sub_int(tmp1,tmp1,1);
            value_division(tmp1,tmp1,D->Ray[r][D->Dimension+1]);
            if (value_lt(n,tmp1)) {
	
                /* New max */
                value_addto(tmp1,D->Ray[r][i+1],D->Ray[r][D->Dimension+1]);
                value_sub_int(tmp1,tmp1,1);
                value_division(n,tmp1,D->Ray[r][D->Dimension+1]);
            }
        }
        value_set_si(c->p[2*i+1][i+1],-1);
        value_assign(c->p[2*i+1][D->Dimension+1],n);
        value_addto(s,c->p[2*i+1][D->Dimension+1],c->p[2*i][D->Dimension+1]);
    
        /* Now test if the dimension of the cube is greater than the size */
        if(value_gt(s,size[i])) {
      
#ifdef EDEBUG62
            fprintf(stderr,"size on dimension %d\n",i);
            fprintf(stderr,"lcm = ");
            for( np=0 ; np<D->Dimension; np++ )
                value_print(stderr,VALUE_FMT,lcm[np]);
            fprintf(stderr,", size = ");
            for( np=0 ; np<D->Dimension; np++ )
                value_print(stderr,VALUE_FMT,size[np]);
            fprintf(stderr,"\n");
            fprintf(stderr,"required size (s) = ");
            value_print(stderr,VALUE_FMT,s);
            fprintf(stderr,"\n");
#endif
      
            /* If the needed size is "small enough"(<=20 or less than 
	       twice *size),
               then increase *size, and artificially increase lcm too !*/
            value_set_si(tmp1,20);
            value_addto(tmp2,size[i],size[i]);
            if(value_le(s,tmp1)
                    || value_le(s,tmp2)) {

                if( value_zero_p(lcm[i]) )
                    value_set_si(lcm[i],1);
                /* lcm divides size... */
                value_division(tmp1,size[i],lcm[i]);
	
                /* lcm = ceil(s/h) */
                value_addto(tmp2,s,tmp1);
                value_add_int(tmp2,tmp2,1);
                value_division(lcm[i],tmp2,tmp1);
	
                /* new size = lcm*h */
                value_multiply(size[i],lcm[i],tmp1);
	
#ifdef EDEBUG62
                fprintf(stderr,"new size = ");
                for( np=0 ; np<D->Dimension; np++ )
                    value_print(stderr,VALUE_FMT,size[np]);
                fprintf(stderr,", new lcm = ");
                for( np=0 ; np<D->Dimension; np++ )
                    value_print(stderr,VALUE_FMT,lcm[np]);
                fprintf(stderr,"\n");
#endif

            }
            else {
	
#ifdef EDEBUG62
                fprintf(stderr,"Failed on dimension %d.\n",i);
#endif
                break;
            }
        }
    }
    if(i!=D->Dimension) {
        Matrix_Free(c);
        value_clear(n); value_clear(s);
        value_clear(tmp1); value_clear(tmp2);
        return(NULL);
    }
    for(i=0;i<D->Dimension;i++) {
        value_subtract(c->p[2*i+1][D->Dimension+1],size[i], 
                c->p[2*i][D->Dimension+1]);
    }
  
#ifdef EDEBUG62
    fprintf(stderr,"PreProcess2 : c =");
    Matrix_Print(stderr,P_VALUE_FMT,c);
#endif

    H = Constraints2Polyhedron(c,MAXRAYS);
    Matrix_Free(c);
    value_clear(n); value_clear(s);
    value_clear(tmp1); value_clear(tmp2);
    return(H);
} /* Polyhedron_Preprocess2 */

/**  

This procedure adds additional constraints to D so that as each
parameter is scanned, it will have a minimum of 'size' points If this
is not possible, an empty polyhedron is returned

@param D
@param size
@param MAXRAYS

*/
Polyhedron *old_Polyhedron_Preprocess(Polyhedron *D,Value size,
				      unsigned MAXRAYS) {
  
    int p, p1, ub, lb;
    Value a, a1, b, b1, g, aa;
    Value abs_a, abs_b, size_copy;
    int dim, con, newi, needed;
    Value **C;
    Matrix *M;
    Polyhedron *D1;
  
    value_init(a); value_init(a1); value_init(b);
    value_init(b1); value_init(g); value_init(aa);
    value_init(abs_a); value_init(abs_b); value_init(size_copy);
  
    dim = D->Dimension;
    con = D->NbConstraints;
    M = Matrix_Alloc(MAXRAYS,dim+2);
    newi = 0;
    value_assign(size_copy,size);
    C = D->Constraint;
    for (p=1; p<=dim; p++) {
        for (ub=0; ub<con; ub++) {
            value_assign(a,C[ub][p]);
            if (value_posz_p(a))	/* a >= 0 */
                continue;		/* not an upper bound */
            for (lb=0;lb<con;lb++) {
                value_assign(b,C[lb][p]);
                if (value_negz_p(b))
                    continue;	/* not a lower bound */
	
                /* Check if a new constraint is needed for this (ub,lb) pair */
                /* a constraint is only needed if a previously scanned */
                /* parameter (1..p-1) constrains this parameter (p) */
                needed=0;
                for (p1=1; p1<p; p1++) {
                    if (value_notzero_p(C[ub][p1]) ||
			value_notzero_p(C[lb][p1])) {
                        needed=1; 
                        break;
                    }
                }
                if (!needed) continue;	
                value_absolute(abs_a,a);
                value_absolute(abs_b,b);

                /* Create new constraint: b*UB-a*LB >= a*b*size */
                Gcd(abs_a,abs_b,&g);
                value_division(a1,a,g);
                value_division(b1,b,g);
                value_set_si(M->p[newi][0],1);
                value_oppose(abs_a,a1);           /* abs_a = -a1 */
                Vector_Combine(&(C[ub][1]),&(C[lb][1]),&(M->p[newi][1]),
                        b1,abs_a,dim+1);
                value_multiply(aa,a1,b1);
		value_addmul(M->p[newi][dim+1], aa, size_copy);
                Vector_Normalize(&(M->p[newi][1]),(dim+1));
                newi++;
            }
        }
    }
    D1 = AddConstraints(M->p_Init,newi,D,MAXRAYS);
    Matrix_Free(M);
  
    value_clear(a); value_clear(a1); value_clear(b);
    value_clear(b1); value_clear(g); value_clear(aa);
    value_clear(abs_a); value_clear(abs_b); value_clear(size_copy);
    return D1;
} /* old_Polyhedron_Preprocess */

/** 

PROCEDURES TO COMPUTE ENUMERATION. recursive procedure, recurse for
each imbriquation

@param pos index position of current loop index (1..hdim-1)
@param P loop domain
@param context context values for fixed indices 
@param res the number of integer points in this
polyhedron 

*/
void count_points (int pos,Polyhedron *P,Value *context, Value *res) {
   
    Value LB, UB, k, c;

    POL_ENSURE_FACETS(P);
    POL_ENSURE_VERTICES(P);

    if (emptyQ(P)) {
	value_set_si(*res, 0);
	return;
    }

    value_init(LB); value_init(UB); value_init(k);
    value_set_si(LB,0);
    value_set_si(UB,0);

    if (lower_upper_bounds(pos,P,context,&LB,&UB) !=0) {
    
        /* Problem if UB or LB is INFINITY */
        fprintf(stderr, "count_points: ? infinite domain\n");
        value_clear(LB); value_clear(UB); value_clear(k);
	value_set_si(*res, -1);
	return;
    }
  
#ifdef EDEBUG1
    if (!P->next) {
        int i;
        for (value_assign(k,LB); value_le(k,UB); value_increment(k,k)) {
            fprintf(stderr, "(");
            for (i=1; i<pos; i++) {
                value_print(stderr,P_VALUE_FMT,context[i]);
                fprintf(stderr,",");
            }
            value_print(stderr,P_VALUE_FMT,k);
            fprintf(stderr,")\n");
        }
    }
#endif
  
    value_set_si(context[pos],0);
    if (value_lt(UB,LB)) {
        value_clear(LB); value_clear(UB); value_clear(k);
	value_set_si(*res, 0);
	return;  
    }  
    if (!P->next) {
        value_subtract(k,UB,LB);
        value_add_int(k,k,1);
	value_assign(*res, k);
        value_clear(LB); value_clear(UB); value_clear(k);
        return;
    } 

    /*-----------------------------------------------------------------*/
    /* Optimization idea                                               */
    /*   If inner loops are not a function of k (the current index)    */
    /*   i.e. if P->Constraint[i][pos]==0 for all P following this and */
    /*        for all i,                                               */
    /*   Then CNT = (UB-LB+1)*count_points(pos+1, P->next, context)    */
    /*   (skip the for loop)                                           */
    /*-----------------------------------------------------------------*/
  
    value_init(c);
    value_set_si(*res, 0);
    for (value_assign(k,LB);value_le(k,UB);value_increment(k,k)) {
        /* Insert k in context */
        value_assign(context[pos],k);
	count_points(pos+1,P->next,context,&c);
	if(value_notmone_p(c))
	    value_addto(*res, *res, c);
	else {
	    value_set_si(*res, -1);
	    break;
        }
    }
    value_clear(c);
  
#ifdef EDEBUG11
    fprintf(stderr,"%d\n",CNT);
#endif
  
    /* Reset context */
    value_set_si(context[pos],0);
    value_clear(LB); value_clear(UB); value_clear(k);
    return;
} /* count_points */

/*-------------------------------------------------------------------*/
/* enode *P_Enum(L, LQ, context, pos, nb_param, dim, lcm,param_name) */
/*     L : list of polyhedra for the loop nest                       */
/*     LQ : list of polyhedra for the parameter loop nest            */
/*     pos : 1..nb_param, position of the parameter                  */
/*     nb_param : number of parameters total                         */
/*     dim : total dimension of the polyhedron, param incl.          */
/*     lcm : denominator array [0..dim-1] of the polyhedron          */
/*     param_name : name of the parameters                           */
/* Returns an enode tree representing the pseudo polynomial          */
/*   expression for the enumeration of the polyhedron.               */
/* A recursive procedure.                                            */
/*-------------------------------------------------------------------*/
static enode *P_Enum(Polyhedron *L,Polyhedron *LQ,Value *context,int pos,
		     int nb_param,int dim,Value *lcm,char **param_name) {

  enode *res,*B,*C;
  int hdim,i,j,rank,flag;
  Value n,g,nLB,nUB,nlcm,noff,nexp,k1,nm,hdv,k,lcm_copy;
  Value tmp;
  Matrix *A;

#ifdef EDEBUG
  fprintf(stderr,"-------------------- begin P_Enum -------------------\n");
  fprintf(stderr,"Calling P_Enum with pos = %d\n",pos);
#endif
  
  /* Initialize all the 'Value' variables */
  value_init(n); value_init(g); value_init(nLB);
  value_init(nUB); value_init(nlcm); value_init(noff);
  value_init(nexp); value_init(k1); value_init(nm);
  value_init(hdv); value_init(k); value_init(tmp);
  value_init(lcm_copy);

  if( value_zero_p(lcm[pos-1]) )
  {
		hdim = 1;
		value_set_si( lcm_copy, 1 );
  }
  else
  {
	  /* hdim is the degree of the polynomial + 1 */
		hdim = dim-nb_param+1; /* homogenous dim w/o parameters */
		value_assign( lcm_copy, lcm[pos-1] );
  }  

  /* code to limit generation of equations to valid parameters only */
  /*----------------------------------------------------------------*/
  flag = lower_upper_bounds(pos,LQ,&context[dim-nb_param],&nLB,&nUB);
  if (flag & LB_INFINITY) {
    if (!(flag & UB_INFINITY)) {
      
      /* Only an upper limit: set lower limit */
      /* Compute nLB such that (nUB-nLB+1) >= (hdim*lcm) */
      value_sub_int(nLB,nUB,1);
      value_set_si(hdv,hdim);
      value_multiply(tmp,hdv,lcm_copy);
      value_subtract(nLB,nLB,tmp);
      if(value_pos_p(nLB))
	value_set_si(nLB,0);
    }
    else {
      value_set_si(nLB,0);
      
      /* No upper nor lower limit: set lower limit to 0 */
      value_set_si(hdv,hdim);
      value_multiply(nUB,hdv,lcm_copy);
      value_add_int(nUB,nUB,1);
    }
  }

  /* if (nUB-nLB+1) < (hdim*lcm) then we have more unknowns than equations */
  /* We can: 1. Find more equations by changing the context parameters, or */
  /* 2. Assign extra unknowns values in such a way as to simplify result.  */
  /* Think about ways to scan parameter space to get as much info out of it*/
  /* as possible.                                                          */
  
#ifdef REDUCE_DEGREE
  if (pos==1 && (flag & UB_INFINITY)==0) {
    /* only for finite upper bound on first parameter */
    /* NOTE: only first parameter because subsequent parameters may
       be artificially limited by the choice of the first parameter */
   
#ifdef EDEBUG
    fprintf(stderr,"*************** n **********\n");
    value_print(stderr,VALUE_FMT,n);
    fprintf(stderr,"\n");
#endif

    value_subtract(n,nUB,nLB);
    value_increment(n,n);
    
#ifdef EDEBUG 
    value_print(stderr,VALUE_FMT,n);
    fprintf(stderr,"\n*************** n ************\n");
#endif 

    /* Total number of samples>0 */
    if(value_neg_p(n))
      i=0;
    else {
      value_modulus(tmp,n,lcm_copy);
      if(value_notzero_p(tmp)) {
		  value_division(tmp,n,lcm_copy);
		  value_increment(tmp,tmp);
		  i = VALUE_TO_INT(tmp);
      }
      else {
		  value_division(tmp,n,lcm_copy);
		  i =   VALUE_TO_INT(tmp);	/* ceiling of n/lcm */
      }
    }

#ifdef EDEBUG 
    value_print(stderr,VALUE_FMT,n);
    fprintf(stderr,"\n*************** n ************\n");
#endif 
    
    /* Reduce degree of polynomial based on number of sample points */
    if (i < hdim){
      hdim=i;
      
#ifdef EDEBUG4
      fprintf(stdout,"Parameter #%d: LB=",pos);
      value_print(stdout,VALUE_FMT,nLB);
      fprintf(stdout," UB=");
      value_print(stdout,VALUE_FMT,nUB);
      fprintf(stdout," lcm=");
      value_print(stdout,VALUE_FMT,lcm_copy);
      fprintf(stdout," degree reduced to %d\n",hdim-1);
#endif
      
    }
  }
#endif /* REDUCE_DEGREE */
  
  /* hdim is now set */  
  /* allocate result structure */
  res = new_enode(polynomial,hdim,pos);
  for (i=0;i<hdim;i++) {
    int l;
    l = VALUE_TO_INT(lcm_copy);
    res->arr[i].x.p = new_enode(periodic,l,pos);
  }
  
  /* Utility arrays */
  A = Matrix_Alloc(hdim, 2*hdim+1);               /* used for Gauss */
  B = new_enode(evector, hdim, 0);
  C = new_enode(evector, hdim, 0);

  /* We'll create these again when we need them */
  for (j = 0; j < hdim; ++j)
    free_evalue_refs(&C->arr[j]);
  
  /*----------------------------------------------------------------*/
  /*                                                                */
  /*                                0<-----+---k--------->          */
  /* |---------noff----------------->-nlcm->-------lcm---->         */
  /* |--- . . . -----|--------------|------+-------|------+-------|-*/
  /* 0          (q-1)*lcm         q*lcm    |  (q+1)*lcm   |         */
  /*                                      nLB          nLB+lcm      */
  /*                                                                */
  /*----------------------------------------------------------------*/
  if(value_neg_p(nLB)) {
    value_modulus(nlcm,nLB,lcm_copy);
    value_addto(nlcm,nlcm,lcm_copy);
  }
  else {
    value_modulus(nlcm,nLB,lcm_copy);
  }
  
  /* noff is a multiple of lcm */
  value_subtract(noff,nLB,nlcm);
  value_addto(tmp,lcm_copy,nlcm);
  for (value_assign(k,nlcm);value_lt(k,tmp);) {
 
#ifdef EDEBUG
    fprintf(stderr,"Finding ");
    value_print(stderr,VALUE_FMT,k);
    fprintf(stderr,"-th elements of periodic coefficients\n");
#endif

    value_set_si(hdv,hdim);
    value_multiply(nm,hdv,lcm_copy);
    value_addto(nm,nm,nLB);
    i=0;
    for (value_addto(n,k,noff); value_lt(n,nm); 
	 value_addto(n,n,lcm_copy),i++) {
      
      /* n == i*lcm + k + noff;   */
      /* nlcm <= k < nlcm+lcm     */
      /* n mod lcm == k1 */
      
#ifdef ALL_OVERFLOW_WARNINGS
      if (((flag & UB_INFINITY)==0) && value_gt(n,nUB)) {
	fprintf(stdout,"Domain Overflow: Parameter #%d:",pos);
	fprintf(stdout,"nLB=");
	value_print(stdout,VALUE_FMT,nLB);
	fprintf(stdout," n=");
	value_print(stdout,VALUE_FMT,n);
	fprintf(stdout," nUB=");
	value_print(stdout,VALUE_FMT,nUB);
	fprintf(stdout,"\n");
      }
#else
      if (overflow_warning_flag && ((flag & UB_INFINITY)==0) &&
	  value_gt(n,nUB)) {
	fprintf(stdout,"\nWARNING: Parameter Domain Overflow.");
	fprintf(stdout," Result may be incorrect on this domain.\n");
	overflow_warning_flag = 0;
      }
#endif
      
      /* Set parameter to n */
      value_assign(context[dim-nb_param+pos],n);
      
#ifdef EDEBUG1
	if( param_name )
	{
      fprintf(stderr,"%s = ",param_name[pos-1]);
      value_print(stderr,VALUE_FMT,n);
      fprintf(stderr," (hdim=%d, lcm[%d]=",hdim,pos-1);
      value_print(stderr,VALUE_FMT,lcm_copy);
      fprintf(stderr,")\n");
	}
	else
	{
      fprintf(stderr,"P%d = ",pos);
      value_print(stderr,VALUE_FMT,n);
      fprintf(stderr," (hdim=%d, lcm[%d]=",hdim,pos-1);
      value_print(stderr,VALUE_FMT,lcm_copy);
      fprintf(stderr,")\n");      
	}
#endif
      
      /* Setup B vector */
      if (pos==nb_param) {
	
#ifdef EDEBUG	
	fprintf(stderr,"Yes\n");
#endif 
	
	/* call count */
	/* count can only be called when the context is fully specified */
	value_set_si(B->arr[i].d,1);
	value_init(B->arr[i].x.n);
	count_points(1,L,context,&B->arr[i].x.n);
	
#ifdef EDEBUG3
	for (j=1; j<pos; j++) fputs("   ",stdout);
	fprintf(stdout, "E(");
	for (j=1; j<nb_param; j++) {
	  value_print(stdout,VALUE_FMT,context[dim-nb_param+j]);
	  fprintf(stdout,",");
	}  
	value_print(stdout,VALUE_FMT,n);
	fprintf(stdout,") = ");
	value_print(stdout,VALUE_FMT,B->arr[i].x.n);
	fprintf(stdout," =");
#endif
	
      }
      else {    /* count is a function of other parameters */
	/* call P_Enum recursively */
	value_set_si(B->arr[i].d,0);
	B->arr[i].x.p = P_Enum(L,LQ->next,context,pos+1,
			       nb_param,dim,lcm,param_name);
	
#ifdef EDEBUG3
	if( param_name )
	{
		for (j=1; j<pos; j++)
		  fputs("   ", stdout);
		fprintf(stdout, "E(");
		for (j=1; j<=pos; j++) {
		  value_print(stdout,VALUE_FMT,context[dim-nb_param+j]);
		  fprintf(stdout,",");
		}	  
		for (j=pos+1; j<nb_param; j++)
		  fprintf(stdout,"%s,",param_name[j]);
		fprintf(stdout,"%s) = ",param_name[j]);
		print_enode(stdout,B->arr[i].x.p,param_name);
		fprintf(stdout," =");
	}
#endif
	
      }
      
      /* Create i-th equation*/
      /* K_0 + K_1 * n**1 + ... + K_dim * n**dim | Identity Matrix */
      
      value_set_si(A->p[i][0],0);  /* status bit = equality */
      value_set_si(nexp,1);
      for (j=1;j<=hdim;j++) {
	value_assign(A->p[i][j],nexp);
	value_set_si(A->p[i][j+hdim],0);
	
#ifdef EDEBUG3
	fprintf(stdout," + ");
	value_print(stdout,VALUE_FMT,nexp);
	fprintf(stdout," c%d",j);
#endif	
	value_multiply(nexp,nexp,n);
      }
      
#ifdef EDEBUG3
      fprintf(stdout, "\n");
#endif
      
      value_set_si(A->p[i][i+1+hdim],1);      
    }
    
    /* Assertion check */
    if (i!=hdim)
      fprintf(stderr, "P_Enum: ?expecting i==hdim\n");
    
#ifdef EDEBUG
	if( param_name )
	{
    fprintf(stderr,"B (enode) =\n");
    print_enode(stderr,B,param_name);
	}
    fprintf(stderr,"A (Before Gauss) =\n");
    Matrix_Print(stderr,P_VALUE_FMT,A);
#endif
    
    /* Solve hdim (=dim+1) equations with hdim unknowns, result in CNT */
    rank = Gauss(A,hdim,2*hdim);
    
#ifdef EDEBUG
    fprintf(stderr,"A (After Gauss) =\n");
    Matrix_Print(stderr,P_VALUE_FMT,A);
#endif
    
    /* Assertion check */
    if (rank!=hdim) {
      fprintf(stderr, "P_Enum: ?expecting rank==hdim\n");
    }
    
    /* if (rank < hdim) then set the last hdim-rank coefficients to ? */
    /* if (rank == 0) then set all coefficients to 0 */    
    /* copy result as k1-th element of periodic numbers */
    if(value_lt(k,lcm_copy))
      value_assign(k1,k);
    else
      value_subtract(k1,k,lcm_copy);
    
    for (i=0; i<rank; i++) {
      
      /* Set up coefficient vector C from i-th row of inverted matrix */
      for (j=0; j<rank; j++) {
	Gcd(A->p[i][i+1],A->p[i][j+1+hdim],&g);
	value_init(C->arr[j].d);
	value_division(C->arr[j].d,A->p[i][i+1],g);
	value_init(C->arr[j].x.n);
	value_division(C->arr[j].x.n,A->p[i][j+1+hdim],g);
      }
      
#ifdef EDEBUG
	if( param_name )
	{
      fprintf(stderr, "C (enode) =\n");
      print_enode(stderr, C, param_name);
	}
#endif
      
      /* The i-th enode is the lcm-periodic coefficient of term n**i */
      edot(B,C,&(res->arr[i].x.p->arr[VALUE_TO_INT(k1)]));
      
#ifdef EDEBUG
	if( param_name )
	{
      fprintf(stderr, "B.C (evalue)=\n");
      print_evalue(stderr,&(res->arr[i].x.p->arr[VALUE_TO_INT(k1)]),
		   param_name);
      fprintf(stderr,"\n");
	}
#endif
      
      for (j = 0; j < rank; ++j)
	free_evalue_refs(&C->arr[j]);
    }
    value_addto(tmp,lcm_copy,nlcm);

    value_increment(k,k);
    for (i = 0; i < hdim; ++i) {
      free_evalue_refs(&B->arr[i]);
      if (value_lt(k,tmp))
	value_init(B->arr[i].d);
    }
  }
  
#ifdef EDEBUG
	if( param_name )
	{
	  fprintf(stderr,"res (enode) =\n");
	  print_enode(stderr,res,param_name);
	}
  fprintf(stderr, "-------------------- end P_Enum -----------------------\n");
#endif
  
  /* Reset context */
  value_set_si(context[dim-nb_param+pos],0);
  
  /* Release memory */
  Matrix_Free(A);
  free(B);
  free(C);
  
  /* Clear all the 'Value' variables */
  value_clear(n); value_clear(g); value_clear(nLB);
  value_clear(nUB); value_clear(nlcm); value_clear(noff);
  value_clear(nexp); value_clear(k1); value_clear(nm);
  value_clear(hdv); value_clear(k); value_clear(tmp);
  value_clear(lcm_copy);
  return res;
} /* P_Enum */

/*----------------------------------------------------------------*/
/* Scan_Vertices(PP, Q, CT)                                       */
/*    PP : ParamPolyhedron                                        */
/*    Q : Domain                                                  */
/*    CT : Context transformation matrix                          */
/*    lcm : lcm array (output)                                    */
/*    nbp : number of parameters                                  */
/*    param_name : name of the parameters                         */
/*----------------------------------------------------------------*/
static void Scan_Vertices(Param_Polyhedron *PP,Param_Domain *Q,Matrix *CT,
   Value *lcm, int nbp, char **param_name )
{
  Param_Vertices *V;
  int i, j, ix, l, np;
  unsigned bx;
  Value k,m1;

  /* Compute the denominator of P */
  /* lcm = Least Common Multiple of the denominators of the vertices of P */
  /* and print the vertices */  

  value_init(k); value_init(m1);
  for( np=0 ; np<nbp ; np++ )
    value_set_si( lcm[np], 0 );

	if( param_name )
	  fprintf(stdout,"Vertices:\n");

  for(i=0,ix=0,bx=MSB,V=PP->V; V && i<PP->nbV; i++,V=V->next) {
    if (Q->F[ix] & bx) {
      if( param_name )
	{
	  if(CT) {
	    Matrix *v;
	    v = VertexCT(V->Vertex,CT);
	    Print_Vertex(stdout,v,param_name);
	    Matrix_Free(v);
	  }
	  else
	    Print_Vertex(stdout,V->Vertex,param_name);
	  fprintf(stdout,"\n");
	}
      
      for(j=0;j<V->Vertex->NbRows;j++) {
	/* A matrix */
	for( l=0 ; l<V->Vertex->NbColumns-1 ; l++ )
	  {
	    if( value_notzero_p(V->Vertex->p[j][l]) )
	      {
		Gcd(V->Vertex->p[j][V->Vertex->NbColumns-1],
		    V->Vertex->p[j][l],&m1);
		value_division(k,V->Vertex->p[j][V->Vertex->NbColumns-1],m1);
		if( value_notzero_p(lcm[l]) )
		  {
		    /* lcm[l] = lcm[l] * k / gcd(k,lcm[l]) */
		    if (value_notzero_p(k) && value_notone_p(k))
		      {
			Gcd(lcm[l],k,&m1);
			value_division(k,k,m1);
			value_multiply(lcm[l],lcm[l],k);
		      }
		  }
		else
		  {
		    value_assign(lcm[l],k);
		  }
	      }
	  }
      }
    }
    NEXT(ix,bx);
  }
  value_clear(k); value_clear(m1);
} /* Scan_Vertices */

/**

Procedure to count points in a non-parameterized polytope.

@param P       Polyhedron to count
@param C       Parameter Context domain
@param CT      Matrix to transform context to original
@param CEq     additionnal equalities in context
@param MAXRAYS workspace size
@param param_name parameter names

*/
Enumeration *Enumerate_NoParameters(Polyhedron *P,Polyhedron *C,
				    Matrix *CT,Polyhedron *CEq,
				    unsigned MAXRAYS,char **param_name) {
  
    Polyhedron *L;
    Enumeration *res;
    Value *context;
    int j;
    int hdim = P->Dimension + 1;
    int r,i;
  
    /* Create a context vector size dim+2 */
    context = (Value *) malloc((hdim+1)*sizeof(Value));
    for (j=0;j<= hdim;j++) 
        value_init(context[j]);
  
    res = (Enumeration *)malloc(sizeof(Enumeration));
    res->next = NULL;
    res->ValidityDomain = Universe_Polyhedron(0);  /* no parameters */
    value_init(res->EP.d);
    value_set_si(res->EP.d,0);  
    L = Polyhedron_Scan(P,res->ValidityDomain,MAXRAYS);
  
#ifdef EDEBUG2
    fprintf(stderr, "L = \n");
    Polyhedron_Print(stderr, P_VALUE_FMT, L);
#endif
  
    if(CT) {
        Polyhedron *Dt;
    
        /* Add parameters to validity domain */
        Dt = Polyhedron_Preimage(res->ValidityDomain,CT,MAXRAYS);
        Polyhedron_Free(res->ValidityDomain);
        res->ValidityDomain = DomainIntersection(Dt,CEq,MAXRAYS);
        Polyhedron_Free(Dt);
    }
  
	if( param_name )
	{
    fprintf(stdout,"---------------------------------------\n");
    fprintf(stdout,"Domain:\n");
    Print_Domain(stdout,res->ValidityDomain, param_name);
  
    /* Print the vertices */
    printf("Vertices:\n");
    for(r=0;r<P->NbRays;++r) {
        if(value_zero_p(P->Ray[r][0]))
            printf("(line) ");
        printf("[");
	if (P->Dimension > 0)
	    value_print(stdout,P_VALUE_FMT,P->Ray[r][1]);
        for(i=1;i<P->Dimension;i++) {
            printf(", ");
            value_print(stdout,P_VALUE_FMT,P->Ray[r][i+1]);
        }  
        printf("]");
        if(value_notone_p(P->Ray[r][P->Dimension+1])) {
            printf("/");
            value_print(stdout,P_VALUE_FMT, P->Ray[r][P->Dimension+1]);
        }
        printf("\n");
    }
	}

  res->EP.x.p = new_enode(polynomial,1,0);;
  value_set_si(res->EP.x.p->arr[0].d, 1);
  value_init(res->EP.x.p->arr[0].x.n);

  if (emptyQ(P)) {
    value_set_si(res->EP.x.p->arr[0].x.n, 0);
  } else if (!L) {
    /* Non-empty zero-dimensional domain */
    value_set_si(res->EP.x.p->arr[0].x.n, 1);
  } else {
    CATCH(overflow_error) {
        fprintf(stderr,"Enumerate: arithmetic overflow error.\n");
        fprintf(stderr,"You should rebuild PolyLib using GNU-MP"
		" or increasing the size of integers.\n");
        overflow_warning_flag = 0;
        assert(overflow_warning_flag);
    }
    TRY {
    
        Vector_Set(context,0,(hdim+1));
    
        /* Set context[hdim] = 1  (the constant) */
        value_set_si(context[hdim],1);
	count_points(1, L, context, &res->EP.x.p->arr[0].x.n);
        UNCATCH(overflow_error);
    }
  }
  
    Domain_Free(L);
  
    /*	**USELESS, there are no references to parameters in res**
	if( CT )
	addeliminatedparams_evalue(&res->EP, CT);
    */
  
	if( param_name )
	{
    fprintf(stdout,"\nEhrhart Polynomial:\n");
    print_evalue(stdout,&res->EP,param_name);
    fprintf(stdout, "\n");
	}

    for (j=0;j<= hdim;j++) 
        value_clear(context[j]);  
    free(context);
    return(res);
} /* Enumerate_NoParameters */

/**  Procedure to count points in a parameterized polytope.

@param Pi Polyhedron to enumerate
@param C Context Domain
@param MAXRAYS size of workspace
@param param_name parameter names (array of strings), may be NULL
@return a list of validity domains + evalues EP

*/
Enumeration *Polyhedron_Enumerate(Polyhedron *Pi,Polyhedron *C,
				  unsigned MAXRAYS,char **param_name)
{
  Polyhedron *L, *CQ, *CQ2, *LQ, *U, *CEq, *rVD, *P, *Ph = NULL;
  Matrix *CT;
  Param_Polyhedron *PP;
  Param_Domain   *Q;
  int i,hdim, dim, nb_param, np;
  Value *lcm, *m1, hdv;
  Value *context;
  Enumeration *en, *res;

  if (POL_ISSET(MAXRAYS, POL_NO_DUAL))
    MAXRAYS = 0;

  POL_ENSURE_FACETS(Pi);
  POL_ENSURE_VERTICES(Pi);
  POL_ENSURE_FACETS(C);
  POL_ENSURE_VERTICES(C);

  res = NULL;
  P = Pi;

#ifdef EDEBUG2
  fprintf(stderr,"C = \n");
  Polyhedron_Print(stderr,P_VALUE_FMT,C);
  fprintf(stderr,"P = \n");
  Polyhedron_Print(stderr,P_VALUE_FMT,P);
#endif

  hdim		= P->Dimension + 1;
  dim		= P->Dimension;
  nb_param	= C->Dimension;

  /* Don't call Polyhedron2Param_Domain if there are no parameters */
  if(nb_param == 0) {
    
    return(Enumerate_NoParameters(P,C,NULL,NULL,MAXRAYS,param_name));  
  }
  if(nb_param == dim) {
    res = (Enumeration *)malloc(sizeof(Enumeration));
    res->next = 0;
    res->ValidityDomain = DomainIntersection(P,C,MAXRAYS);
    value_init(res->EP.d);
    value_init(res->EP.x.n);
    value_set_si(res->EP.d,1);
    value_set_si(res->EP.x.n,1);
    if( param_name ) {
      fprintf(stdout,"---------------------------------------\n");
      fprintf(stdout,"Domain:\n");
      Print_Domain(stdout,res->ValidityDomain, param_name);
      fprintf(stdout,"\nEhrhart Polynomial:\n");
      print_evalue(stdout,&res->EP,param_name);
      fprintf(stdout, "\n");
    }
    return res;
  }
  PP = Polyhedron2Param_SimplifiedDomain(&P,C,MAXRAYS,&CEq,&CT);
  if(!PP) {
	if( param_name )
    fprintf(stdout, "\nEhrhart Polynomial:\nNULL\n");

    return(NULL);
  }
  
  /* CT : transformation matrix to eliminate useless ("false") parameters */
  if(CT) {
    nb_param -= CT->NbColumns-CT->NbRows;
    dim  -= CT->NbColumns-CT->NbRows;
    hdim -= CT->NbColumns-CT->NbRows;
    
    /* Don't call Polyhedron2Param_Domain if there are no parameters */
    if(nb_param == 0) {
	res = Enumerate_NoParameters(P,C,CT,CEq,MAXRAYS,param_name);
	Param_Polyhedron_Free(PP);
	goto out;
    }  
  }

  /* get memory for Values */
  lcm = (Value *)malloc((nb_param+1) * sizeof(Value));
  m1  = (Value *)malloc((nb_param+1) * sizeof(Value));
  /* Initialize all the 'Value' variables */
  for( np=0 ; np < nb_param+1; np++ )
  {
    value_init(lcm[np]); value_init(m1[np]);
  }
  value_init(hdv);

  for(Q=PP->D;Q;Q=Q->next) {
    int hom = 0;
    if(CT) {
      Polyhedron *Dt;
      CQ = Q->Domain;      
      Dt = Polyhedron_Preimage(Q->Domain,CT,MAXRAYS);
      rVD = DomainIntersection(Dt,CEq,MAXRAYS);
      
      /* if rVD is empty or too small in geometric dimension */
      if(!rVD || emptyQ(rVD) ||
	 (rVD->Dimension-rVD->NbEq < Dt->Dimension-Dt->NbEq-CEq->NbEq)) {
	if(rVD)
	  Polyhedron_Free(rVD);
	Polyhedron_Free(Dt);
	Polyhedron_Free(CQ);
	continue;		/* empty validity domain */
      }
      Polyhedron_Free(Dt);
    }
    else
      rVD = CQ = Q->Domain;    
    en = (Enumeration *)malloc(sizeof(Enumeration));
    en->next = res;
    res = en;
    res->ValidityDomain = rVD;
    
    if( param_name ) {
      fprintf(stdout,"---------------------------------------\n");
      fprintf(stdout,"Domain:\n");
    
#ifdef EPRINT_ALL_VALIDITY_CONSTRAINTS
      Print_Domain(stdout,res->ValidityDomain,param_name);
#else    
      {
	Polyhedron *VD;
	VD = DomainSimplify(res->ValidityDomain,C,MAXRAYS);
	Print_Domain(stdout,VD,param_name);
	Domain_Free(VD);
      }
#endif /* EPRINT_ALL_VALIDITY_CONSTRAINTS */
    }
    
    overflow_warning_flag = 1;
    
    /* Scan the vertices and compute lcm */
    Scan_Vertices(PP,Q,CT,lcm,nb_param,param_name);
    
#ifdef EDEBUG2
    fprintf(stderr,"Denominator = ");
    for( np=0;np<nb_param;np++)
      value_print(stderr,P_VALUE_FMT,lcm[np]);
    fprintf(stderr," and hdim == %d \n",hdim);
#endif
    
#ifdef EDEBUG2
    fprintf(stderr,"CQ = \n");
    Polyhedron_Print(stderr,P_VALUE_FMT,CQ);
#endif

    /* Before scanning, add constraints to ensure at least hdim*lcm */
    /* points in every dimension */
    value_set_si(hdv,hdim-nb_param);

    for( np=0;np<nb_param+1;np++) {
	if( value_notzero_p(lcm[np]) )
	    value_multiply(m1[np],hdv,lcm[np]);
	else
	    value_set_si(m1[np],1);
    }

#ifdef EDEBUG2 
    fprintf(stderr,"m1 == ");
    for( np=0;np<nb_param;np++)
      value_print(stderr,P_VALUE_FMT,m1[np]);
    fprintf(stderr,"\n");
#endif 

    CATCH(overflow_error) {
      fprintf(stderr,"Enumerate: arithmetic overflow error.\n");
      CQ2 = NULL;
    }
    TRY {      
      CQ2 = Polyhedron_Preprocess(CQ,m1,MAXRAYS);

#ifdef EDEBUG2
      fprintf(stderr,"After preprocess, CQ2 = ");
      Polyhedron_Print(stderr,P_VALUE_FMT,CQ2);
#endif
      
      UNCATCH(overflow_error);
    }

    /* Vin100, Feb 2001 */
    /* in case of degenerate, try to find a domain _containing_ CQ */
    if ((!CQ2 || emptyQ(CQ2)) && CQ->NbBid==0) {
      int r;
      
#ifdef EDEBUG2
      fprintf(stderr,"Trying to call Polyhedron_Preprocess2 : CQ = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,CQ);
#endif
      
      /* Check if CQ is bounded */
      for(r=0;r<CQ->NbRays;r++) {
	if(value_zero_p(CQ->Ray[r][0]) ||
	   value_zero_p(CQ->Ray[r][CQ->Dimension+1]))
	  break;
      }
      if(r==CQ->NbRays) {
	
	/* ok, CQ is bounded */
	/* now find if CQ is contained in a hypercube of size m1 */
	CQ2 = Polyhedron_Preprocess2(CQ,m1,lcm,MAXRAYS);
      }
    }

    if (!CQ2) {
      Polyhedron *tmp;
#ifdef EDEBUG2
      fprintf(stderr,"Homogenize.\n");
#endif
      hom = 1;
      tmp = homogenize(CQ, MAXRAYS);
      CQ2 = Polyhedron_Preprocess(tmp,m1,MAXRAYS);
      Polyhedron_Free(tmp);
      if (!Ph)
	Ph = homogenize(P, MAXRAYS);
      for (np=0; np < nb_param+1; np++)
	  if (value_notzero_p(lcm[np]))
	      value_addto(m1[np],m1[np],lcm[np]);
    }
    
    if (!CQ2 || emptyQ(CQ2)) {
#ifdef EDEBUG2
      fprintf(stderr,"Degenerate.\n");
#endif
      fprintf(stdout,"Degenerate Domain. Can not continue.\n");
      value_init(res->EP.d);
      value_init(res->EP.x.n);
      value_set_si(res->EP.d,1);
      value_set_si(res->EP.x.n,-1);
    }
    else {
      
#ifdef EDEBUG2
      fprintf(stderr,"CQ2 = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,CQ2);
      if( ! PolyhedronIncludes(CQ, CQ2) )
	fprintf( stderr,"CQ does not include CQ2 !\n");
      else
	fprintf( stderr,"CQ includes CQ2.\n");
      if( ! PolyhedronIncludes(res->ValidityDomain, CQ2) )
	fprintf( stderr,"CQ2 is *not* included in validity domain !\n");
      else
	fprintf( stderr,"CQ2 is included in validity domain.\n");
#endif
      
      /* L is used in counting the number of points in the base cases */
      L = Polyhedron_Scan(hom ? Ph : P,CQ2,MAXRAYS);
      U = Universe_Polyhedron(0);
      
      /* LQ is used to scan the parameter space */
      LQ = Polyhedron_Scan(CQ2,U,MAXRAYS); /* bounds on parameters */
      Domain_Free(U);
      if(CT)	/* we did compute another Q->Domain */
	Domain_Free(CQ);
      
      /* Else, CQ was Q->Domain (used in res) */
      Domain_Free(CQ2);
      
#ifdef EDEBUG2
      fprintf(stderr,"L = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,L);
      fprintf(stderr,"LQ = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,LQ);
#endif
#ifdef EDEBUG3
      fprintf(stdout,"\nSystem of Equations:\n");
#endif
      
      value_init(res->EP.d);
      value_set_si(res->EP.d,0);
      
      /* Create a context vector size dim+2 */
      context = (Value *) malloc ((hdim+1+hom)*sizeof(Value));  
      for(i=0;i<=(hdim+hom);i++)
	value_init(context[i]);
      Vector_Set(context,0,(hdim+1+hom));
      
      /* Set context[hdim] = 1  (the constant) */
      value_set_si(context[hdim+hom],1);
      
      CATCH(overflow_error) {
	  fprintf(stderr,"Enumerate: arithmetic overflow error.\n");
	  fprintf(stderr,"You should rebuild PolyLib using GNU-MP "
		  "or increasing the size of integers.\n");
	  overflow_warning_flag = 0;
	  assert(overflow_warning_flag);
	
      }
      TRY {
	  res->EP.x.p = P_Enum(L,LQ,context,1,nb_param+hom,
			       dim+hom,lcm,param_name);
	  UNCATCH(overflow_error);	
      }
      if (hom)
	  dehomogenize_evalue(&res->EP, nb_param+1);
      
      for(i=0;i<=(hdim+hom);i++)
	value_clear(context[i]);
      free(context);
      Domain_Free(L);
      Domain_Free(LQ);
      
#ifdef EDEBUG5
	if( param_name )
   {
      fprintf(stdout,"\nEhrhart Polynomial (before simplification):\n");
      print_evalue(stdout,&res->EP,param_name);
	}
#endif

      /* Try to simplify the result */
      reduce_evalue(&res->EP);
      
      /* Put back the original parameters into result */
      /* (equalities have been eliminated)            */
      if(CT) 
          addeliminatedparams_evalue(&res->EP,CT);
      
      if (param_name) {
	fprintf(stdout,"\nEhrhart Polynomial:\n");
	print_evalue(stdout,&res->EP, param_name);
	fprintf(stdout,"\n");
	/* sometimes the final \n lacks (when a single constant is printed) */
      }
      
    }
  }

  if (Ph)
    Polyhedron_Free(Ph);
  /* Clear all the 'Value' variables */
  for (np=0; np < nb_param+1; np++)
  {
    value_clear(lcm[np]); value_clear(m1[np]);
  }
  value_clear(hdv);
  free(lcm);
  free(m1);
  /* We can't simply call Param_Polyhedron_Free because we've reused the domains */
  Param_Vertices_Free(PP->V);
  while (PP->D) {
    Q = PP->D;
    PP->D = PP->D->next;
    free(Q->F);
    free(Q);
  }
  free(PP);

out:
  if (CEq)
    Polyhedron_Free(CEq);
  if (CT)
    Matrix_Free(CT);
  if( P != Pi )
    Polyhedron_Free( P );

  return res;
} /* Polyhedron_Enumerate */ 

void Enumeration_Free(Enumeration *en)
{
  Enumeration *ee;

  while( en )
  {
          free_evalue_refs( &(en->EP) );
          Domain_Free( en->ValidityDomain );
          ee = en ->next;
          free( en );
          en = ee;
  }
}

/* adds by B. Meister for Ehrhart Polynomial approximation */

/**
 * Divides the evalue e by the integer n<br/>
 * recursive function<br/>
 * Warning :  modifies e
 * @param e an evalue (to be divided by n)
 * @param n
 */

void evalue_div(evalue * e, Value n) {
  int i;
  Value gc;
  value_init(gc);
  if (value_zero_p(e->d)) {
    for (i=0; i< e->x.p->size; i++) {
      evalue_div(&(e->x.p->arr[i]), n);
    }
  }
  else {
    value_multiply(e->d, e->d, n);
    /* simplify the new rational if needed */
    Gcd(e->x.n, e->d, &gc);
    if (value_notone_p(gc)&&(value_notzero_p(gc))) {
      value_division(e->d, e->d, gc);
      value_division(e->x.n, e->x.n, gc);
    }
  }
  value_clear(gc);
} /* evalue_div */


/** 

Ehrhart_Quick_Apx_Full_Dim(P, C, MAXRAYS, param_names)

Procedure to estimate the nubmer of points in a parameterized polytope.
Returns a list of validity domains + evalues EP
B.M.
The most rough and quick approximation by variables expansion  
Deals with the full-dimensional case.                          
@param Pi : Polyhedron to enumerate (approximatively)
@param C : Context Domain
@param MAXRAYS : size of workspace
@param param_name : names for the parameters (char strings)

*/
Enumeration *Ehrhart_Quick_Apx_Full_Dim(Polyhedron *Pi,Polyhedron *C,
					unsigned MAXRAYS, char ** param_name)
{
  Polyhedron *L, *CQ, *CQ2, *LQ, *U, *CEq, *rVD, *P;
  Matrix *CT;
  Param_Polyhedron *PP;
  Param_Domain   *Q;
  int i,j,hdim, dim, nb_param, np;
  Value *lcm, *m1, hdv;
  Value *context;
  Enumeration *en, *res;
  unsigned int nb_vars;
  Matrix * denoms;
  Value expansion_det;
  Value global_var_lcm;
  Matrix * expansion;
  Polyhedron * Expanded;

  /* used to scan the vertices */
  Param_Vertices * V_tmp;

  res = NULL;
  P = Pi;


  value_init(expansion_det);
  value_init(global_var_lcm);

#ifdef EDEBUG2
  fprintf(stderr,"C = \n");
  Polyhedron_Print(stderr,P_VALUE_FMT,C);
  fprintf(stderr,"P = \n");
  Polyhedron_Print(stderr,P_VALUE_FMT,P);
#endif

  hdim		= P->Dimension + 1;
  dim		= P->Dimension;
  nb_param	= C->Dimension;

  /* Don't call Polyhedron2Param_Domain if there are no parameters */
  if(nb_param == 0) {
    
    return(Enumerate_NoParameters(P,C,NULL,NULL,MAXRAYS, param_name));  
  }

#if EDEBUG2
  printf("Enumerating polyhedron : \n");
  Polyhedron_Print(stdout, P_VALUE_FMT, P);
#endif

  PP = Polyhedron2Param_SimplifiedDomain(&P,C,MAXRAYS,&CEq,&CT);
  if(!PP) {
    if( param_name )
      fprintf(stdout, "\nEhrhart Polynomial:\nNULL\n");

    return(NULL);
  }
  
  /* CT : transformation matrix to eliminate useless ("false") parameters */
  if(CT) {
    nb_param -= CT->NbColumns-CT->NbRows;
    dim  -= CT->NbColumns-CT->NbRows;
    hdim -= CT->NbColumns-CT->NbRows;
    
    /* Don't call Polyhedron2Param_Domain if there are no parameters */
    if(nb_param == 0)
    {
	    res = Enumerate_NoParameters(P,C,CT,CEq,MAXRAYS,param_name);
	    if( P != Pi )
		    Polyhedron_Free( P );
	    return( res );
    }  
  }

  /* get memory for Values */
  lcm = (Value *)malloc( nb_param * sizeof(Value));
  m1  = (Value *)malloc( nb_param * sizeof(Value));
  /* Initialize all the 'Value' variables */
  for( np=0 ; np<nb_param ; np++ )
  {
    value_init(lcm[np]); value_init(m1[np]);
  }
  value_init(hdv);

  /* Scan the vertices and make an orthogonal expansion of the variable
     space */
  /* a- prepare the array of common denominators */
  if (!PP->nbV) return 0;
  else {
    nb_vars = P->Dimension-nb_param;
    denoms = Matrix_Alloc(1, nb_vars);
    for (i=0; i< nb_vars; i++) value_set_si(denoms->p[0][i], 0);
  }
  
  /* b- scan the vertices and compute the variables' global lcms */
  for (V_tmp = PP->V; V_tmp; V_tmp=V_tmp->next)
    for (i=0; i< nb_vars; i++) 
      Lcm3(denoms->p[0][i],V_tmp->Vertex->p[i][nb_param+1], &(denoms->p[0][i]));
  printf("denoms = \n");
  Matrix_Print(stderr, P_VALUE_FMT, denoms);
  value_set_si(expansion_det, 1);
  value_set_si(global_var_lcm, 1);
  for (i=0; i< nb_vars;i++) {
    value_multiply(expansion_det, expansion_det, denoms->p[0][i]);
    Lcm3(global_var_lcm, denoms->p[0][i], &global_var_lcm);
  }
  printf("expansion_det:\n");
  value_print(stderr, P_VALUE_FMT, expansion_det);
  printf("\n");
  
  
  /* the expansion can be actually writen as global_var_lcm.L^{-1} */
  /* this is equivalent to multiply the rows of P by denoms_det */
  for (i=0; i< nb_vars; i++) value_division(denoms->p[0][i], 
					    global_var_lcm, denoms->p[0][i]);
  
  /* OPT : we could use a vector instead of a diagonal matrix here (c- and d-).*/
  /* c- make the quick expansion matrix */
  printf("nb vars = %d, nb param = %d", nb_vars, nb_param);
  expansion = Matrix_Alloc(nb_vars+nb_param+1, nb_vars+nb_param+1);
  for (i=0; i< nb_vars; i++) {
    for (j=0; j< nb_vars+nb_param+1; j++) {
      if (i==j) value_assign(expansion->p[i][j], denoms->p[0][i]);
      else value_set_si(expansion->p[i][j], 0);
    }
  }
  for (i=nb_vars; i< nb_vars+nb_param+1; i++) {
    for (j=0; j< nb_vars+nb_param+1; j++) {
      if (i==j) value_assign(expansion->p[i][j], global_var_lcm);
      else value_set_si(expansion->p[i][j], 0);
    }
  }
  value_clear(global_var_lcm);
  printf("expansion = \n");
  Matrix_Print(stderr, P_VALUE_FMT, expansion);
  
  /* d- apply the variable expansion to the polyhedron */
#if EDEBUG2 
  Polyhedron_Print(stderr, P_VALUE_FMT, P);
#endif

  Expanded = Polyhedron_Preimage(P, expansion, MAXRAYS);
  
#if EDEBUG2
  Polyhedron_Print(stderr, P_VALUE_FMT, Expanded);
#endif

  Polyhedron_Free(P);
  P = Expanded;
  /* formerly : Scan the vertices and compute lcm 
     Scan_Vertices_Quick_Apx(PP,Q,CT,lcm,nb_param); */
  /* now : lcm = 1 (by construction) */
  /* OPT : A lot among what happens after this point can be simplified by
     knowing that lcm[i] = 1 for now, we just conservatively fool the rest of
     the function with lcm = 1 */
  for (i=0; i< nb_param; i++) value_set_si(lcm[i], 1);
  
  for(Q=PP->D;Q;Q=Q->next) {
    if(CT) {
      Polyhedron *Dt;
      CQ = Q->Domain;      
      Dt = Polyhedron_Preimage(Q->Domain,CT,MAXRAYS);
      rVD = DomainIntersection(Dt,CEq,MAXRAYS);
      
      /* if rVD is empty or too small in geometric dimension */
      if(!rVD || emptyQ(rVD) ||
	 (rVD->Dimension-rVD->NbEq < Dt->Dimension-Dt->NbEq-CEq->NbEq)) {
	if(rVD)
	  Polyhedron_Free(rVD);
	Polyhedron_Free(Dt);
	continue;		/* empty validity domain */
      }
      Polyhedron_Free(Dt);
    }
    else
      rVD = CQ = Q->Domain;
    en = (Enumeration *)malloc(sizeof(Enumeration));
    en->next = res;
    res = en;
    res->ValidityDomain = rVD;
	
    if( param_name )
      {    
	fprintf(stdout,"---------------------------------------\n");
	fprintf(stdout,"Domain:\n");
    
#ifdef EPRINT_ALL_VALIDITY_CONSTRAINTS
	Print_Domain(stdout,res->ValidityDomain,param_name);
#else    
	{
	  Polyhedron *VD;
	  VD = DomainSimplify(res->ValidityDomain,C,MAXRAYS);
	  Print_Domain(stdout,VD,param_name);
	  Domain_Free(VD);
	}
#endif /* EPRINT_ALL_VALIDITY_CONSTRAINTS */
      }
    
    overflow_warning_flag = 1;

#ifdef EDEBUG2
    fprintf(stderr,"Denominator = ");
    for( np=0;np<nb_param;np++)
      value_print(stderr,P_VALUE_FMT,lcm[np]);
    fprintf(stderr," and hdim == %d \n",hdim);
#endif
    
#ifdef EDEBUG2
    fprintf(stderr,"CQ = \n");
    Polyhedron_Print(stderr,P_VALUE_FMT,CQ);
#endif

    /* Before scanning, add constraints to ensure at least hdim*lcm */
    /* points in every dimension */
    value_set_si(hdv,hdim-nb_param);

    for( np=0;np<nb_param;np++)
    {
		if( value_notzero_p(lcm[np]) )
			value_multiply(m1[np],hdv,lcm[np]);
		else
			value_set_si(m1[np],1);
    }

#ifdef EDEBUG2 
    fprintf(stderr,"m1 == ");
    for( np=0;np<nb_param;np++)
      value_print(stderr,P_VALUE_FMT,m1[np]);
    fprintf(stderr,"\n");
#endif 

    CATCH(overflow_error) {
      fprintf(stderr,"Enumerate: arithmetic overflow error.\n");
      CQ2 = NULL;
    }
    TRY {      
      CQ2 = Polyhedron_Preprocess(CQ,m1,MAXRAYS);

#ifdef EDEBUG2
      fprintf(stderr,"After preprocess, CQ2 = ");
      Polyhedron_Print(stderr,P_VALUE_FMT,CQ2);
#endif
      
      UNCATCH(overflow_error);
    }
    
    /* Vin100, Feb 2001 */
    /* in case of degenerate, try to find a domain _containing_ CQ */
    if ((!CQ2 || emptyQ(CQ2)) && CQ->NbBid==0) {
      int r;
      
#ifdef EDEBUG2
      fprintf(stderr,"Trying to call Polyhedron_Preprocess2 : CQ = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,CQ);
#endif
      
      /* Check if CQ is bounded */
      for(r=0;r<CQ->NbRays;r++) {
	if(value_zero_p(CQ->Ray[r][0]) ||
	   value_zero_p(CQ->Ray[r][CQ->Dimension+1]))
	  break;
      }
      if(r==CQ->NbRays) {
	
	/* ok, CQ is bounded */
	/* now find if CQ is contained in a hypercube of size m1 */
	CQ2 = Polyhedron_Preprocess2(CQ,m1,lcm,MAXRAYS);
      }
    }
    if (!CQ2 || emptyQ(CQ2)) {
#ifdef EDEBUG2
      fprintf(stderr,"Degenerate.\n");
#endif
      fprintf(stdout,"Degenerate Domain. Can not continue.\n");
      value_init(res->EP.d);
      value_init(res->EP.x.n);
      value_set_si(res->EP.d,1);
      value_set_si(res->EP.x.n,-1);
    }
    else {
      
#ifdef EDEBUG2
      fprintf(stderr,"CQ2 = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,CQ2);
      if( ! PolyhedronIncludes(CQ, CQ2) )
	fprintf( stderr,"CQ does not include CQ2 !\n");
      else
	fprintf( stderr,"CQ includes CQ2.\n");
      if( ! PolyhedronIncludes(res->ValidityDomain, CQ2) )
	fprintf( stderr,"CQ2 is *not* included in validity domain !\n");
      else
	fprintf( stderr,"CQ2 is included in validity domain.\n");
#endif
      
      /* L is used in counting the number of points in the base cases */
      L = Polyhedron_Scan(P,CQ,MAXRAYS);
      U = Universe_Polyhedron(0);
      
      /* LQ is used to scan the parameter space */
      LQ = Polyhedron_Scan(CQ2,U,MAXRAYS); /* bounds on parameters */
      Domain_Free(U);
      if(CT)	/* we did compute another Q->Domain */
	Domain_Free(CQ);
      
      /* Else, CQ was Q->Domain (used in res) */
      Domain_Free(CQ2);
      
#ifdef EDEBUG2
      fprintf(stderr,"L = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,L);
      fprintf(stderr,"LQ = \n");
      Polyhedron_Print(stderr,P_VALUE_FMT,LQ);
#endif
#ifdef EDEBUG3
      fprintf(stdout,"\nSystem of Equations:\n");
#endif
      
      value_init(res->EP.d);
      value_set_si(res->EP.d,0);
      
      /* Create a context vector size dim+2 */
      context = (Value *) malloc ((hdim+1)*sizeof(Value));  
      for(i=0;i<=(hdim);i++)
	value_init(context[i]);
      Vector_Set(context,0,(hdim+1));
      
      /* Set context[hdim] = 1  (the constant) */
      value_set_si(context[hdim],1);
      
      CATCH(overflow_error) {
	fprintf(stderr,"Enumerate: arithmetic overflow error.\n");
	fprintf(stderr,"You should rebuild PolyLib using GNU-MP "
		"or increasing the size of integers.\n");
	overflow_warning_flag = 0;
	assert(overflow_warning_flag);
	
      }
      TRY {
	res->EP.x.p = P_Enum(L,LQ,context,1,nb_param,dim,lcm, param_name);
	UNCATCH(overflow_error);	
      }
      
      for(i=0;i<=(hdim);i++)
	value_clear(context[i]);
      free(context);
      Domain_Free(L);
      Domain_Free(LQ);
      
#ifdef EDEBUG5
      if( param_name )
	{
	  fprintf(stdout,"\nEhrhart Polynomial (before simplification):\n");
	  print_evalue(stdout,&res->EP,param_name);
	}

      /* BM: divide EP by denom_det, the expansion factor */
      fprintf(stdout,"\nEhrhart Polynomial (before division):\n");
      print_evalue(stdout,&(res->EP),param_name);
#endif

      evalue_div(&(res->EP), expansion_det);


      /* Try to simplify the result */
      reduce_evalue(&res->EP);
      
      /* Put back the original parameters into result */
      /* (equalities have been eliminated)            */
      if(CT) 
          addeliminatedparams_evalue(&res->EP,CT);
      
      if( param_name )
	{
	  fprintf(stdout,"\nEhrhart Polynomial:\n");
	  print_evalue(stdout,&res->EP, param_name);
	  fprintf(stdout,"\n");
	  /* sometimes the final \n lacks (when a single constant is printed)
	  */
	}
      
    }
  }
  value_clear(expansion_det);

  if( P != Pi )
 	Polyhedron_Free( P );
  /* Clear all the 'Value' variables */
  for( np=0; np<nb_param ; np++ )
  {
    value_clear(lcm[np]); value_clear(m1[np]);
  }
  value_clear(hdv);

  return res;
} /* Ehrhart_Quick_Apx_Full_Dim */ 


/** 
 * Computes the approximation of the Ehrhart polynomial of a polyhedron 
 * (implicit form -> matrix), treating the non-full-dimensional case.
 * @param M a polyhedron under implicit form
 * @param C  M's context under implicit form
 * @param Validity_Lattice a pointer to the parameter's validity lattice
 * @param MAXRAYS the needed "working space" for other polylib functions used here
 * @param param_name the names of the parameters, 
 */
Enumeration *Ehrhart_Quick_Apx(Matrix * M, Matrix * C, 
			       Matrix **Validity_Lattice, unsigned maxRays) { 
			       /* char ** param_name) {*/

  /* 0- compute a full-dimensional polyhedron with the same number of points,
     and its parameter's validity lattice */
  Matrix * M_full;
  Polyhedron * P, * PC;
  Enumeration *en;
  
  M_full = full_dimensionize(M, C->NbColumns-2, Validity_Lattice);
  /* 1- apply the same tranformation to the context that what has been applied
     to the parameters space of the polyhedron. */
  mpolyhedron_compress_last_vars(C, *Validity_Lattice);
  show_matrix(M_full);
  P = Constraints2Polyhedron(M_full, maxRays);
  PC = Constraints2Polyhedron(C, maxRays);
  Matrix_Free(M_full);
  /* compute the Ehrhart polynomial of the "equivalent" polyhedron */
  en = Ehrhart_Quick_Apx_Full_Dim(P, PC, maxRays, NULL);

  /* clean up */
  Polyhedron_Free(P);
  Polyhedron_Free(PC);
  return en;
} /* Ehrhart_Quick_Apx */



/** 
 * Computes the approximation of the Ehrhart polynomial of a polyhedron (implicit
 * form -> matrix), treating the non-full-dimensional case.  If there are some
 * equalities involving only parameters, these are used to eliminate useless 
 * parameters.
 * @param M a polyhedron under implicit form
 * @param C  M's context under implicit form
 * @param validityLattice a pointer to the parameter's validity lattice
 * (returned)
 * @param parmsEqualities Equalities involving only the parameters
 * @param maxRays the needed "working space" for other polylib functions used here
 * @param elimParams  array of the indices of the eliminated parameters (returned)
 */
Enumeration *Constraints_EhrhartQuickApx(Matrix const * M, Matrix const * C, 
					 Matrix ** validityLattice, 
					 Matrix ** parmEqualities,
					 unsigned int ** elimParms, 
					 unsigned maxRays) {
  Enumeration *EP;
  Matrix * Mp = Matrix_Copy(M);
  Matrix * Cp = Matrix_Copy(C);
  /* remove useless equalities involving only parameters, using these 
     equalities to remove parameters. */
  (*parmEqualities) = Constraints_Remove_parm_eqs(&Mp, &Cp, 0, elimParms);
  if (Mp->NbRows>=0) {/* if there is no contradiction */
    EP = Ehrhart_Quick_Apx(Mp, Cp, validityLattice, maxRays);
    return EP;
  }
  else {
    /* if there are contradictions, return a zero Ehrhart polynomial */
    return NULL;
  }
  
}


/** 
 * Returns the array of parameter names after some of them have been eliminated.
 * @param parmNames the initial names of the parameters
 * @param elimParms a list of parameters that have been eliminated from the
 * original parameters. The first element of this array is the number of
 * elements.
 * <p> Note: does not copy the parameters names themselves. </p>
 * @param nbParms the initial number of parameters
*/
char ** parmsWithoutElim(char const ** parmNames, 
			       unsigned int const * elimParms, 
			       unsigned int nbParms) {
  int i=0, j=0,k;
  int newParmNb = nbParms - elimParms[0];
  char ** newParmNames = (char **)malloc (newParmNb * sizeof(char *));
  for (k=1; k<= elimParms[0]; k++) {
    while (i!=elimParms[k]) {
      newParmNames[i-k+1] = (char *) parmNames[i];
      i++;
    }
  }
  return newParmNames;
}


/** 
 * returns a constant Ehrhart polynomial whose value is zero for any value of
 * the parameters.
 * @param nbParms the number of parameters, i.e., the number of arguments to
 * the Ehrhart polynomial
 */
Enumeration * Enumeration_zero(unsigned int nbParms, unsigned int maxRays) {
  Matrix * Mz = Matrix_Alloc(1, nbParms+3);
  Polyhedron * emptyP;
  Polyhedron * universe;
  Enumeration * zero;
  /* 1- build an empty polyhedron with the right dimension */
  /* here we choose to take 2i = -1 */
  value_set_si(Mz->p[0][1], 2);
  value_set_si(Mz->p[0][nbParms+2], 1);
  emptyP = Constraints2Polyhedron(Mz, maxRays);
  Matrix_Free(Mz);
  universe = Universe_Polyhedron(nbParms);
  zero = Polyhedron_Enumerate(emptyP, universe, maxRays, NULL);
  Polyhedron_Free(emptyP);
  Polyhedron_Free(universe);
  return zero;
} /* Enumeration_zero() */

