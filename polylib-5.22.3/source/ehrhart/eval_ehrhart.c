/************************************************/
/*  eval_ehrhart.c                              */
/* functions to evaluate an Ehrhart polynomial. */
/* written by Emmanuel Jeannot (c) 1997.        */
/*  Emmanuel.Jeannot@ens-lyon.fr                */
/*  http://www.ens-lyon.fr/~ejeannot            */
/*                                              */
/* modified 1998, 2000, Vincent Loechner        */
/* (ArithmetiqueLib, Param_Names)               */
/************************************************/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

#include <polylib/polylib.h>

/* #define EVAL_EHRHART_DEBUG  */

/********************************************************/
/* function in domain                                   */
/*    check if the parameters in list_args              */
/*    verifies the constraints of Domain P              */
/********************************************************/
int in_domain(Polyhedron *P, Value *list_args) {
  
  int col,row;
  Value v; /* value of the constraint of a row when
               parameters are instanciated*/

  if( !P )
          return( 0 );
  POL_ENSURE_FACETS(P);
  POL_ENSURE_VERTICES(P);

  value_init(v); 
  
  /* P->Constraint constraint matrice of polyhedron P */  
  for(row=0;row<P->NbConstraints;row++) {
    value_assign(v,P->Constraint[row][P->Dimension+1]); /*constant part*/
    for(col=1;col<P->Dimension+1;col++) {
      value_addmul(v, P->Constraint[row][col], list_args[col-1]); 
    }  
    if (value_notzero_p(P->Constraint[row][0])) {
        
      /*if v is not >=0 then this constraint is not respected */
      if (value_neg_p(v)) {
        value_clear(v);
        return( in_domain(P->next, list_args) );
      }        
    }
    else {
      
      /*if v is not = 0 then this constraint is not respected */
      if (value_notzero_p(v)) {
        value_clear(v);
        return( in_domain(P->next, list_args) );
      }
    }
  }
  
  /* if not return before this point => all the constraints are respected */
  value_clear(v);
  return 1;
} /* in_domain */

/****************************************************/
/* function compute enode                           */
/*     compute the value of enode p with parameters */
/*     list "list_args                              */
/*     compute the polynomial or the periodic       */
/****************************************************/

static double compute_enode(enode *p, Value *list_args) {
  
  int i;
  Value m, param;
  double res=0.0;
    
  if (!p)
    return(0.);

  value_init(m);
  value_init(param);

  if (p->type == polynomial) {
    if (p->size > 1)
                 value_assign(param,list_args[p->pos-1]);
    
    /* Compute the polynomial using Horner's rule */
    for (i=p->size-1;i>0;i--) {
      res +=compute_evalue(&p->arr[i],list_args);
      res *=VALUE_TO_DOUBLE(param);
    }
    res +=compute_evalue(&p->arr[0],list_args);
  }
  else if (p->type == periodic) {
    value_assign(m,list_args[p->pos-1]);
    
    /* Choose the right element of the periodic */
    value_set_si(param,p->size);
    value_pmodulus(m,m,param);
    res = compute_evalue(&p->arr[VALUE_TO_INT(m)],list_args);
  }
  value_clear(m);
  value_clear(param);
  return res;
} /* compute_enode */

/*************************************************/
/* return the value of Ehrhart Polynomial        */
/* It returns a double, because since it is      */
/* a recursive function, some intermediate value */
/* might not be integral                         */
/*************************************************/

double compute_evalue(evalue *e,Value *list_args) {
  
  double res;
  
  if (value_notzero_p(e->d)) {
    if (value_notone_p(e->d)) 
      res = VALUE_TO_DOUBLE(e->x.n) / VALUE_TO_DOUBLE(e->d);
    else 
      res = VALUE_TO_DOUBLE(e->x.n);
  }
  else 
    res = compute_enode(e->x.p,list_args);
  return res;
} /* compute_evalue */


/****************************************************/
/* function compute_poly :                          */
/* Check for the good validity domain               */
/* return the number of point in the Polyhedron     */
/* in allocated memory                              */
/* Using the Ehrhart pseudo-polynomial              */
/****************************************************/
Value *compute_poly(Enumeration *en,Value *list_args) {

  Value *tmp;
  /*        double d; int i; */

  tmp = (Value *) malloc (sizeof(Value));
  assert(tmp != NULL);
  value_init(*tmp);
  value_set_si(*tmp,0);

  if(!en)
    return(tmp);        /* no ehrhart polynomial */
  if(en->ValidityDomain) {
    if(!en->ValidityDomain->Dimension) { /* no parameters */
      value_set_double(*tmp,compute_evalue(&en->EP,list_args)+.25);
      return(tmp);
    }
  }  
  else 
    return(tmp);  /* no Validity Domain */    
  while(en) {
    if(in_domain(en->ValidityDomain,list_args)) {
      
#ifdef EVAL_EHRHART_DEBUG
      Print_Domain(stdout,en->ValidityDomain,NULL);
      print_evalue(stdout,&en->EP,NULL);
#endif
      
      /*                        d = compute_evalue(&en->EP,list_args);
                                i = d;
                                printf("(double)%lf = %d\n", d, i ); */
      value_set_double(*tmp,compute_evalue(&en->EP,list_args)+.25);
      return(tmp);
    }
    else
      en=en->next;
  }
  value_set_si(*tmp,0);
  return(tmp); /* no compatible domain with the arguments */
} /* compute_poly */ 




