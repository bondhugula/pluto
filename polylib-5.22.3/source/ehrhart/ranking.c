/**
 * Tools to compute the ranking function of an iteration J: the number of
 * integer points in P that are lexicographically inferior to J 
 * B. Meister
 * 6/2005
 * LSIIT-ICPS, UMR 7005 CNRS Université Louis Pasteur
 * HiPEAC Network
 */

#include <polylib/polylib.h>
#include <polylib/ranking.h>

/**
 * Returns a list of polytopes needed to compute
 * the number of points in P that are lexicographically
 * smaller than a given point in D.
 * Only the first dim dimensions are taken into account
 * for computing the lexsmaller relation.
 * The remaining variables are assumed to be extra
 * existential/control variables.
 * When P == D, this is the conventional ranking function.
 * P and D are assumed to have the same parameter domain C.
 *
 * The first polyhedron in the list returned is the
 * updated context: a combination of D and C or an extended C.
 *
 * The order of the variables in the remaining polyhedra is
 * - first dim variables of P
 * - existential variables of P
 * - existential variables of D
 * - first dim variables of D
 * - the parameters
 */
Polyhedron *LexSmaller(Polyhedron *P, Polyhedron *D, unsigned dim,
			Polyhedron *C, unsigned MAXRAYS)
{
  unsigned i, j, k, r;
  unsigned nb_parms = C->Dimension;
  unsigned nb_vars = dim;
  unsigned P_extra = P->Dimension - nb_vars - nb_parms;
  unsigned D_extra = D->Dimension - nb_vars - nb_parms;
  unsigned nb_new_parms;
  unsigned ncons;
  Matrix * cur_element, * C_times_J, * Klon;
  Polyhedron * P1, *C1;
  Polyhedron * lexico_lesser_union = NULL;

  POL_ENSURE_INEQUALITIES(C);
  POL_ENSURE_INEQUALITIES(D);
  POL_ENSURE_INEQUALITIES(P);

  assert(P->Dimension >= C->Dimension + dim);
  assert(D->Dimension >= C->Dimension + dim);
  nb_new_parms = nb_vars;

  /* the number of variables must be positive */
  if (nb_vars<=0) {
    printf("\nRanking > No variables, returning NULL.\n"); 
    return NULL;
  }
  /*
   * if D has extra variables, then we can't squeeze the contraints
   * of D in the new context, so we simply add them to each element.
   */
  if (D_extra)
    cur_element = Matrix_Alloc(P->NbConstraints+D->NbConstraints+nb_new_parms, 
			       P->Dimension+D_extra+nb_new_parms+2);
  else
    cur_element = Matrix_Alloc(P->NbConstraints+nb_new_parms, 
			       P->Dimension+D_extra+nb_new_parms+2);


  /* 0- Put P in the first rows of cur_element */
  for (i=0; i < P->NbConstraints; i++) {
    Vector_Copy(P->Constraint[i], cur_element->p[i], nb_vars+P_extra+1);
    Vector_Copy(P->Constraint[i]+1+nb_vars+P_extra, 
		cur_element->p[i]+1+nb_vars+P_extra+D_extra+nb_new_parms, 
		nb_parms+1);
  }
  ncons = P->NbConstraints;
  if (D_extra) {
    for (i=0; i < D->NbConstraints; i++) {
      r = P->NbConstraints + i;
      Vector_Copy(D->Constraint[i], cur_element->p[r], 1);
      Vector_Copy(D->Constraint[i]+1, 
		  cur_element->p[r]+1+nb_vars+P_extra+D_extra, nb_new_parms);
      Vector_Copy(D->Constraint[i]+1+nb_new_parms, 
		  cur_element->p[r]+1+nb_vars+P_extra, D_extra);
      Vector_Copy(D->Constraint[i]+1+nb_new_parms+D_extra, 
		  cur_element->p[r]+1+nb_vars+P_extra+D_extra+nb_new_parms, 
		  nb_parms+1);
    }
    ncons += D->NbConstraints;
  }

  /* 1- compute the Ehrhart polynomial of each disjoint polyhedron defining the
     lexicographic order */
  for (k=0, r = ncons; k < nb_vars; k++, r++) {

    /* a- build the corresponding matrix
     *  the nb of rows of cur_element is fake, so that we do not have to
     *  re-allocate it. */
    cur_element->NbRows = r+1;

    /* convert the previous (strict) inequality into an equality */
    if (k>=1) {
      value_set_si(cur_element->p[r-1][0], 0);
      value_set_si(cur_element->p[r-1][cur_element->NbColumns-1], 0);
    }
    /* build the k-th inequality from P */
    value_set_si(cur_element->p[r][0], 1);
    value_set_si(cur_element->p[r][k+1], -1);
    value_set_si(cur_element->p[r][nb_vars+P_extra+D_extra+k+1], 1);
    /* we want a strict inequality */
    value_set_si(cur_element->p[r][cur_element->NbColumns-1], -1);
#ifdef ERDEBUG
    show_matrix(cur_element);
#endif

    /* b- add it to the current union
       as Constraints2Polyhedron modifies its input, we must clone cur_element */
    Klon = Matrix_Copy(cur_element);
    P1 = Constraints2Polyhedron(Klon, MAXRAYS);
    Matrix_Free(Klon);
    P1->next = lexico_lesser_union;
    lexico_lesser_union = P1;
  }
  
  /* 2- as we introduce n parameters, we must introduce them into the context
   * as well. 
   * The added constraints are P.M.(J N 1 )^T >=0 */
  if (D_extra)
    C_times_J = Matrix_Alloc(C->NbConstraints, nb_new_parms+nb_parms+2);
  else
    C_times_J = Matrix_Alloc(C->NbConstraints + D->NbConstraints, D->Dimension+2);
  /* copy the initial context while adding the new parameters */
  for (i = 0; i < C->NbConstraints; i++) {
    value_assign(C_times_J->p[i][0], C->Constraint[i][0]);
    Vector_Copy(C->Constraint[i]+1, C_times_J->p[i]+1+nb_new_parms, nb_parms+1);
  }

  /* copy constraints from evaluation domain */
  if (!D_extra)
    for (i = 0; i < D->NbConstraints; i++)
      Vector_Copy(D->Constraint[i], C_times_J->p[C->NbConstraints+i], 
		  D->Dimension+2);

#ifdef ERDEBUG
  show_matrix(C_times_J);
#endif
  C1 = Constraints2Polyhedron(C_times_J, POL_NO_DUAL);

  /* 4- clean up */
  Matrix_Free(cur_element);
  Matrix_Free(C_times_J);

  C1->next = P1;

  return C1;
} /* LexSmaller */


/**
 * Returns the number of points in P that are lexicographically
 * smaller than a given point in D.
 * Only the first dim dimensions are taken into account
 * for computing the lexsmaller relation.
 * The remaining variables are assumed to be extra
 * existential/control variables.
 * When P == D, this is the conventional ranking function.
 * P and D are assumed to have the same parameter domain C.
 * The variables in the Enumeration correspond to the first dim variables
 * in D followed by the parameters of D (the variables of C).
 */
Enumeration *Polyhedron_LexSmallerEnumerate(Polyhedron *P, Polyhedron *D, 
					    unsigned dim,
					    Polyhedron *C, unsigned MAXRAYS)
{
  Enumeration * ranking;
  Polyhedron *RC, *RD;

  RC = LexSmaller(P, D, dim, C, MAXRAYS);
  RD = RC->next;
  RC->next = NULL;

  /* Compute the ranking, which is the sum of the Ehrhart polynomials of the n
     disjoint polyhedra we just put in P1. */
  /* OPT : our polyhdera are (already) disjoint, so Domain_Enumerate does
     probably too much work uselessly */
  ranking = Domain_Enumerate(RD, RC, MAXRAYS, NULL);

  Domain_Free(RD);
  Polyhedron_Free(RC);

  return ranking;
}


/*
 * Returns a function that assigns a unique number to each point in the
 * polytope P ranging from zero to (number of points in P)-1.
 * The order of the numbers corresponds to the lexicographical order.
 *
 * C is the parameter context of the polytope
 */
Enumeration *Polyhedron_Ranking(Polyhedron *P, Polyhedron *C, unsigned MAXRAYS)
{
    return Polyhedron_LexSmallerEnumerate(P, P, P->Dimension-C->Dimension, 
					  C, MAXRAYS);
}
