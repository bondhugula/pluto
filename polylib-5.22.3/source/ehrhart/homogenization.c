/** homogenization.c 
    copyright 2004-2005 Bavo Nootaert
**/

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <polylib/polylib.h>
#include <polylib/homogenization.h>

static evalue *dehomogenize_periodic(enode *en);
static evalue *dehomogenize_polynomial(enode *en);

Polyhedron *homogenize(Polyhedron *P, unsigned MAXRAYS)
{
    Matrix M, *M2;
    /* Pretend P is a Matrix for a second */
    M.NbRows = P->NbConstraints;
    M.NbColumns = P->Dimension+2;
    M.p_Init = P->p_Init;
    M.p = P->Constraint;
    M2 = AddANullColumn(&M);
    P = Constraints2Polyhedron(M2, MAXRAYS);
    Matrix_Free(M2);
    return P;
}

/** dehomogenize an evalue. The last parameter (nb_param) is replaced by 1.
    This function is mutually recursive with dehomogenize_enode.
**/
void dehomogenize_evalue(evalue *ep, int nb_param){
  evalue *w;

  /** cannot dehomogenize rationals **/
  if (value_zero_p(ep->d)){

    /** we need to replace the last parameter **/
    if (ep->x.p->pos == nb_param){
      if (ep->x.p->type == periodic && ep->x.p->size > 1){
	w = dehomogenize_periodic(ep->x.p); 
      }
      else{
	w = dehomogenize_polynomial(ep->x.p);
      }
      free_evalue_refs(ep);
      memcpy(ep, w, sizeof(evalue));
      free(w);
    }
    else{
      /** Not the last parameter. Recurse **/
      dehomogenize_enode(ep->x.p, nb_param);
    }

  }
}

/** dehomogenize all evalues in an enode. 
    This function is mutually recursive with dehomogenize_evalue.
**/
void dehomogenize_enode(enode *p, int nb_param){
  evalue *temp;
  int i;
  for (i = 0; i < p->size; i++){
    dehomogenize_evalue(&p->arr[i], nb_param);
  }
}


/** return the 1st element of an enode representing a periodic **/
static evalue *dehomogenize_periodic(enode *en){
  evalue *w;
  assert(en->type == periodic);
  assert(en->size > 1);
  assert(value_notzero_p(en->arr[1].d));
  w = (evalue*)malloc(sizeof(evalue));
  value_init(w->d); value_init(w->x.n);
  value_assign(w->d, en->arr[1].d); value_assign(w->x.n, en->arr[1].x.n);
  return w;
}

/** dehomogenize a polynomial. Assume the enode contains a polynomial in 
    one variable, the homogenous parameter. 
    Returns an new evalue, representing a rational.
 **/
static evalue *dehomogenize_polynomial(enode *en){
  evalue *enn;
  evalue *ev;
  int i;
  double som;
  Value num, den, gcd, f1, f2;
  assert(en->type == polynomial);
  som = 0;
  value_init(num); value_init(den); value_init(gcd);
  value_init(f1); value_init(f2);
  value_set_si(den, 1);

  /** enumerate over all coefficients (which are either periodic or rational,
      but not polynomial) **/
  for (i = 0; i < en->size; i++){
    if (value_zero_p(en->arr[i].d)){
      if (en->arr[i].x.p->size > 1)
	ev = &en->arr[i].x.p->arr[1];
      else
	ev = &en->arr[i].x.p->arr[0];
    }
    else{
      ev = &en->arr[i];
    }
    /** add ev (fraction) to num/den **/
    value_multiply(f1, den, ev->x.n);
    value_multiply(f2, num, ev->d);
    value_addto(num, f1, f2);
    value_multiply(den, den, ev->d);
  }
  
  /** simplify num/den **/
  Gcd(num, den, &gcd);
  value_division(num, num, gcd);
  value_division(den, den, gcd);

  /** create new evalue representing num/den**/
  enn = (evalue*)malloc(sizeof(evalue));
  value_init(enn->d); value_init(enn->x.n);
  value_assign(enn->d, den);
  value_assign(enn->x.n, num);

  /** cleanup **/
  value_clear(gcd);
  value_clear(f1); value_clear(f2); 
  value_clear(num); value_clear(den);

  return enn;
}

/** dehomogenize a polyhedron. Assume the polyhedron p is homogenous.
    Returns a new polyhedron.
**/
Polyhedron *dehomogenize_polyhedron(Polyhedron *p, int maxRays){
  Matrix *constr, *constrh;
  Polyhedron *ph;
  int i;
  constr = Polyhedron2Constraints(p);
  constrh = Matrix_Alloc(constr->NbRows, constr->NbColumns - 1);
  for (i = 0; i < constr->NbRows; i++){
    Vector_Copy(constr->p[i], constrh->p[i], constr->NbColumns - 1);
  }
  ph = Constraints2Polyhedron(constrh, maxRays);
  Matrix_Free(constr); Matrix_Free(constrh);
  return ph;
}

/** dehomogenize an enumeration. Replaces each validity domain and 
    Ehrhart polynomial in the Enumeration en with the dehomogenized form.
 **/
void dehomogenize_enumeration(Enumeration* en, int nb_params, int maxRays){
  Enumeration *en2;
  Polyhedron *vd;
  for (en2 = en; en2; en2 = en2->next) {
    vd = dehomogenize_polyhedron(en2->ValidityDomain, maxRays);
    Polyhedron_Free(en2->ValidityDomain);
    en2->ValidityDomain = vd;
    dehomogenize_evalue(&en2->EP, nb_params);
  }
}
