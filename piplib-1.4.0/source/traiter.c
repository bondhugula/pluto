/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 traiter.c                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
 *                                                                            *
 * This library is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at    *
 * your option) any later version.                                            *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "pip.h"

#define max(x,y) ((x) > (y)? (x) : (y))

extern long int cross_product, limit;
extern int verbose;
extern FILE *dump;
extern int profondeur;
extern int compa_count;

#if !defined(LINEAR_VALUE_IS_MP)
int llog(Entier x)
{int n = 0;
/* x must be positive, you dummy */
 if(x<0) x=-x;
 while(x) x >>= 1, n++;
 return(n);
}
#endif

int chercher(Tableau *p, int masque, int n)
{int i;
 for(i = 0; i<n; i++)
     if(p->row[i].flags & masque) break;
 return(i);
}

/* il est convenu que traiter ne doit modifier ni le tableau, ni le contexte;
   le tableau peut grandir en cas de coupure (+1 en hauteur et +1 en largeur
   si nparm != 0) et en cas de partage (+1 en hauteur)(seulement si nparm != 0).
   le contexte peut grandir en cas de coupure (+2 en hauteur et +1 en largeur)
   (seulement si nparm !=0) et en cas de partage (+1 en hauteur)(nparm !=0).
   On estime le nombre de coupures a llog(D) et le nombre de partages a
   ni.
*/

Tableau *expanser(Tableau *tp, int virt, int reel, int ncol, 
                               int off, int dh, int dw)
{
 int i, j, ff;
 char *q; Entier *pq;
 Entier *pp, *qq;
 Tableau *rp;
 if(tp == NULL) return(NULL);
 rp = tab_alloc(reel+dh, ncol+dw, virt);

 #if defined(LINEAR_VALUE_IS_MP)
 mpz_set(rp->determinant, tp->determinant);
 #else
 rp->l_determinant = tp->l_determinant;
 for(i=0; i<tp->l_determinant; i++)
     rp->determinant[i] = tp->determinant[i];
 #endif
 pq = (Entier *) & (rp->row[virt+reel+dh]);
 for(i = off; i<virt + reel; i++)
     {ff = Flag(rp, i) = Flag(tp, i-off);
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_set(Denom(rp, i), Denom(tp, i-off));
      #else
      Denom(rp, i) = Denom(tp, i-off);
      #endif
      if(ff & Unit) rp->row[i].objet.unit = tp->row[i-off].objet.unit;
      else {
	  rp->row[i].objet.val = pq;
	  pq +=(ncol + dw);
	  pp = tp->row[i-off].objet.val;
	  qq = rp->row[i].objet.val;
	  for(j = 0; j<ncol; j++)
             #if defined(LINEAR_VALUE_IS_MP)
	     mpz_set(*qq++, *pp++);
             #else
	     *qq++ = *pp++;
             #endif
	  }
      }
 return(rp);
}

/* Check for "obvious" signs of the parametric constant terms
 * of the inequalities.  As soon as a negative sign is found
 * we return from this function and handle this constraint
 * in the calling function.  The signs of the other constraints
 * are then mostly irrelevant.
 * If any of the negative signs is due to the "big parameter",
 * then we want to use this constraint first.
 * We therefore check for signs determined by the coefficient
 * of the big parameter first.
 */
int exam_coef(Tableau *tp, int nvar, int ncol, int bigparm)
{int i, j ;
 int ff, fff;
 #if defined(LINEAR_VALUE_IS_MP)
 int x;
 #else
 Entier x;
 #endif
 Entier *p;
 
 if (bigparm >= 0)
    for (i = 0; i<tp->height; i++) {
	if (Flag(tp, i) != Unknown)
	    continue;
	x = entier_sgn(Index(tp,i, bigparm));
	if (x < 0) {
	    Flag(tp, i) = Minus;
	    return i;
	} else if (x > 0)     
	    Flag(tp, i) = Plus;
    }

 for(i = 0; i<tp->height; i++)
     {ff = Flag(tp,i);
      if(ff == 0) break;
      if(ff == Unknown) {
	   ff = Zero;
	   p = &(tp->row[i].objet.val[nvar+1]);
	   for(j = nvar+1; j<ncol; j++) {
                #if defined(LINEAR_VALUE_IS_MP)
	        x = mpz_sgn(*p); p++ ;
	        #else
	        x = *p++;
                #endif
		if(x<0) fff = Minus;
		else if (x>0) fff = Plus;
		else fff = Zero;
		if(fff != Zero && fff != ff)
		    if(ff == Zero) ff = fff;
		    else {ff = Unknown;
			  break;
			 }
	       }
/* bug de'tecte' par [paf], 16/2/93 !
   Si tous les coefficients des parame`tres sont ne'gatifs
   et si le terme constant est nul, le signe est inconnu!!
   On traite donc spe'cialement le terme constant. */
           #if defined(LINEAR_VALUE_IS_MP)
	   x = mpz_sgn(Index(tp, i, nvar));
	   #else
	   x = Index(tp, i, nvar);
           #endif
	   if(x<0) fff = Minus;
	   else if(x>0) fff = Plus;
	   else fff = Zero;
/* ici on a le signe du terme constant */
	   switch(ff){
/* le signe est inconnu si les coefficients sont positifs et
   le terme constant ne'gatif */
	   case Plus: if(fff == Minus) ff = Unknown; break;
/* si les coefficients sont tous nuls, le signe est celui
   du terme constant */
	   case Zero: ff = fff; break;
/* le signe est inconnu si les coefficients sont ne'gatifs,
   sauf si le terme constant est egalement negatif. */
	   case Minus: if(fff != Minus) ff = Unknown; break;
/* enfin, il n'y a rien a` dire si le signe des coefficients est inconnu */
	   }
	   Flag(tp, i) = ff;
	   if(ff == Minus) return(i);
	  }
      }
 return(i);
}

void compa_test(Tableau *tp, Tableau *context,
		int ni, int nvar, int nparm, int nc)
{
 int i, j;
 int ff;
 int cPlus, cMinus, isCritic;
 int verbold;
 Tableau *tPlus, *tMinus;
 int p;
 struct high_water_mark q;

 if(nparm == 0) return;
 if(nparm >= MAXPARM) {
     fprintf(stderr, "Too much parameters : %d\n", nparm);
     exit(1);
     }
 q = tab_hwm();

 for(i = 0; i<ni + nvar; i++)
     {ff = Flag(tp,i);
      if(ff & (Critic | Unknown))
	  {isCritic = Pip_True;
	   for(j = 0; j<nvar; j++)
                 #if defined(LINEAR_VALUE_IS_MP)
		 if(mpz_sgn(Index(tp, i, j)) > 0)
                 #else
	         if(Index(tp, i, j) > 0)
                 #endif
		 {isCritic = Pip_False;
		  break;
		 }
           compa_count++;
	   tPlus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tPlus, nparm+nc) = Unknown;
	   for (j = 0; j < nparm; j++)
	       entier_assign(Index(tPlus, nparm+nc, j), Index(tp, i, j+nvar+1));
	   entier_assign(Index(tPlus, nparm+nc, nparm), Index(tp, i, nvar));
	   if (!isCritic)
	       entier_decrement(Index(tPlus, nparm+nc, nparm),
				    Index(tPlus, nparm+nc, nparm));
	   entier_assign(Denom(tPlus, nparm+nc), UN);
	   
	   p = sol_hwm();
	   traiter(tPlus, NULL, nparm, 0, nc+1, 0, -1, TRAITER_INT);
	   cPlus = is_not_Nil(p);
	   if(verbose>0){
	     fprintf(dump, "\nThe positive case has been found ");
	     fprintf(dump, cPlus? "possible\n": "impossible\n");
	     fflush(dump);
	   }

	   sol_reset(p);
	   tMinus = expanser(context, nparm, nc, nparm+1, nparm, 1, 0);
	   Flag(tMinus, nparm+nc) = Unknown;
	   for (j = 0; j < nparm; j++)
	       entier_oppose(Index(tMinus, nparm+nc, j), Index(tp, i, j+nvar+1));
	   entier_oppose(Index(tMinus, nparm+nc, nparm), Index(tp, i, nvar));
	   entier_decrement(Index(tMinus, nparm+nc, nparm),
				Index(tMinus, nparm+nc, nparm));
	   entier_assign(Denom(tMinus, nparm+nc), UN);
	   traiter(tMinus, NULL, nparm, 0, nc+1, 0, -1, TRAITER_INT);
	   cMinus = is_not_Nil(p);
	   if(verbose>0){
	     fprintf(dump, "\nThe negative case has been found ");
	     fprintf(dump, cMinus? "possible\n": "impossible\n");
	     fflush(dump);
	   }

	   sol_reset(p);
	   if (cPlus && cMinus) {
	       Flag(tp,i) = isCritic ? Critic : Unknown;
	     }
	   else if (cMinus)
	      {Flag(tp,i) = Minus;
	       break;
	      }
	   else {
	     Flag(tp,i) = cPlus ? Plus : Zero;
	   }
	  }
     }
 tab_reset(q);
 
 return;
}

Entier *valeur(Tableau *tp, int i, int j)
{
 if(Flag(tp, i) & Unit)
     return(tp->row[i].objet.unit == j ? &Denom(tp,i) : &ZERO);
 else return(&Index(tp, i, j));
}

void solution(Tableau *tp, int nvar, int nparm)
{int i, j;
 int ncol = nvar + nparm + 1;

 sol_list(nvar);
 for(i = 0; i<nvar; i++)
     {sol_forme(nparm+1);
      for(j = nvar+1; j<ncol; j++)
	 sol_val(*valeur(tp, i, j), Denom(tp,i));
      sol_val(*valeur(tp, i, nvar), Denom(tp,i));
     }
}

static void solution_dual(Tableau *tp, int nvar, int nparm, int *pos)
{
    int i;

    sol_list(tp->height - nvar);
    for (i = 0; i < tp->height - nvar; ++i) {
	sol_forme(1);
	if (Flag(tp, pos[i]) & Unit)
	    sol_val(*valeur(tp, 0, tp->row[pos[i]].objet.unit), Denom(tp, 0));
	else
	    sol_val(ZERO, UN);
    }
}

int choisir_piv(Tableau *tp, int pivi, int nvar, int nligne)
{
 int j, k;
 Entier pivot, foo, x, y;
 int sgn_x, pivj = -1;

 #if defined(LINEAR_VALUE_IS_MP)
 mpz_init(pivot); mpz_init(foo); mpz_init(x); mpz_init(y);
 #endif
 
 for(j = 0; j<nvar; j++) {
    #if defined(LINEAR_VALUE_IS_MP)
    mpz_set(foo, Index(tp, pivi, j));
    if(mpz_sgn(foo) <= 0) continue;
    if(pivj < 0)
	{pivj = j;
         mpz_set(pivot, foo);
	 continue;
	}
    for(k = 0; k<nligne; k++)
        {mpz_mul(x, pivot, *valeur(tp, k, j)); 
         mpz_mul(y, *valeur(tp, k, pivj), foo);
         mpz_sub(x, x, y);
         cross_product++;
         sgn_x = mpz_sgn(x);
         if(sgn_x) break;
	}
    if(sgn_x < 0)
        {pivj = j;
         mpz_set(pivot, foo);
        }
    #else
    if((foo = Index(tp, pivi, j)) <= 0) continue;
    if(pivj < 0)
	{pivj = j;
	 pivot = foo;
	 continue;
	}
    for(k = 0; k<nligne; k++)
	{x = pivot * (*valeur(tp, k, j)) - (*valeur(tp, k, pivj)) * foo;
	 cross_product++;
	 if(x) break;
	}
    if(x < 0)
	{pivj = j;
	 pivot = foo;
	}
    #endif
 }
 
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_clear(pivot); mpz_clear(foo); mpz_clear(x); mpz_clear(y);
 #endif

 return(pivj);
}


int pivoter(Tableau *tp, int pivi, int nvar, int nparm, int ni)

{int pivj;
 int ncol = nvar + nparm + 1;
 int nligne = nvar + ni;
 int i, j, k;
 Entier x, y, d, gcd, u, dpiv;
 int ff, fff;
 Entier pivot, foo, z;
 Entier ppivot, dppiv;
 Entier new[MAXCOL], *p, *q;
 Entier lpiv;
 int sgn_x;
 #if !defined(LINEAR_VALUE_IS_MP)
 char format_format[32];

 sprintf(format_format, "\nPivot %s/%s\n", FORMAT, FORMAT);
 #endif

 if(ncol >= MAXCOL) {
   fprintf(stdout, "Too much variables\n");
   exit(1);
 }
 if(0 > pivi || pivi >= nligne || Flag(tp, pivi) == Unit) {
   fprintf(stdout, "Syserr : pivoter : wrong pivot row\n");
   exit(1);
 }

 pivj = choisir_piv(tp, pivi, nvar, nligne);
 if(pivj < 0) return(-1);
 if(pivj >= nvar) {
   fprintf(stdout, "Syserr : pivoter : wrong pivot\n");
   exit(1);
 }

 #if defined(LINEAR_VALUE_IS_MP)
 mpz_init(x); mpz_init(y); mpz_init(d); 
 mpz_init(gcd); mpz_init(u); mpz_init(dpiv);
 mpz_init(lpiv); mpz_init(pivot); mpz_init(foo);
 mpz_init(z); mpz_init(ppivot); mpz_init(dppiv);

 for(i=0; i<ncol; i++)
   mpz_init(new[i]);

 mpz_set(pivot, Index(tp, pivi, pivj));
 mpz_set(dpiv, Denom(tp, pivi));
 mpz_gcd(d, pivot, dpiv);
 mpz_divexact(ppivot, pivot, d);
 mpz_divexact(dppiv, dpiv, d);
 #else
 pivot = Index(tp, pivi, pivj);
 dpiv = Denom(tp, pivi);
 d = pgcd(pivot, dpiv);
 ppivot = pivot/d;
 dppiv = dpiv/d;
 #endif
 
 if(verbose>1){
   #if defined(LINEAR_VALUE_IS_MP)
   fprintf(dump, "Pivot ");
   mpz_out_str(dump, 10, ppivot);
   putc('/', dump);
   mpz_out_str(dump, 10, dppiv);
   putc('\n', dump);
   #else
   fprintf(dump, format_format, ppivot, dppiv);
   #endif
   fprintf(dump, "%d x %d\n", pivi, pivj);
 }

 #if defined(LINEAR_VALUE_IS_MP)
 mpz_fdiv_qr(x, y, tp->determinant, dppiv); 
 #else
 for(i=0; i< tp->l_determinant; i++){
     d=pgcd(tp->determinant[i], dppiv);
     tp->determinant[i] /= d;
     dppiv /= d;
     }
 #endif

 #if defined(LINEAR_VALUE_IS_MP)
 if(mpz_sgn(y) != 0){ 
 #else
 if(dppiv != 1) {
 #endif
   fprintf(stderr, "Integer overflow\n");
   if(verbose>0) fflush(dump);
   exit(1);
 }
 
 #if defined(LINEAR_VALUE_IS_MP)
 mpz_mul(tp->determinant, x, ppivot);
 #else
 for(i=0; i<tp->l_determinant; i++)
     if(llog(tp->determinant[i]) + llog(ppivot) < 8*sizeof(Entier)){
	 tp->determinant[i] *= ppivot;
	 break;
	 }
 if(i >= tp->l_determinant){
     tp->l_determinant++;
     if(tp->l_determinant >= MAX_DETERMINANT){
	 fprintf(stderr, "Integer overflow : %d\n", tp->l_determinant);
	 exit(1);
	 }
     tp->determinant[i] = ppivot;
     }
 #endif

 if(verbose>1){
   fprintf(dump, "determinant ");
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_out_str(dump, 10, tp->determinant);
   #else
   for(i=0; i<tp->l_determinant; i++)
	fprintf(dump, FORMAT, tp->determinant[i]);
   #endif
   fprintf(dump, "\n");
 }

 
 for(j = 0; j<ncol; j++)
   #if defined(LINEAR_VALUE_IS_MP)
   if(j==pivj)
     mpz_set(new[j], dpiv);
   else 
     mpz_neg(new[j], Index(tp, pivi, j));
   #else
   new[j] = (j == pivj ? dpiv : -Index(tp, pivi, j));
   #endif

 for(k = 0; k<nligne; k++){
   if(Flag(tp,k) & Unit)continue;
   if(k == pivi)continue;
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_set(foo, Index(tp, k, pivj));
   mpz_gcd(d, pivot, foo);
   mpz_divexact(lpiv, pivot, d);
   mpz_divexact(foo, foo, d);
   mpz_set(d, Denom(tp,k));
   mpz_mul(gcd, lpiv, d);
   mpz_set(Denom(tp, k), gcd);
   #else
   foo = Index(tp, k, pivj);
   d = pgcd(pivot, foo);
   lpiv = pivot/d;
   foo /= d;
   d = Denom(tp,k);
   gcd = lpiv * d;
   Denom(tp, k) = gcd;
   #endif
   p = tp->row[k].objet.val;
   q = tp->row[pivi].objet.val;
   for(j = 0; j<ncol; j++){
     if(j == pivj)
     #if defined(LINEAR_VALUE_IS_MP)
       mpz_mul(z, dpiv, foo);
     #else
       z = dpiv * foo;
     #endif
     else {
     #if defined(LINEAR_VALUE_IS_MP)
       mpz_mul(z, *p, lpiv);
       mpz_mul(y, *q, foo);
       mpz_sub(z, z, y);
     #else
       z = (*p) * lpiv - (*q) * foo;
     #endif
     }
     q++;
     cross_product++;
     #if defined(LINEAR_VALUE_IS_MP)
     mpz_set(*p, z);
     p++;
     if(mpz_cmp_ui(gcd, 1) != 0)
       mpz_gcd(gcd, gcd, z);
     #else
     *p++ = z;
     if(gcd != 1)
       gcd = pgcd(gcd, z);
     #endif
   }
   #if defined(LINEAR_VALUE_IS_MP)
   if(mpz_cmp_ui(gcd, 1) != 0){
     p = tp->row[k].objet.val;
     for(j = 0; j<ncol; j++){
       mpz_divexact(*p, *p, gcd);
       p++;
     }
   }
   mpz_divexact(Denom(tp,k), Denom(tp,k), gcd);
   #else
   if(gcd != 1) {
    p = tp->row[k].objet.val;
    for(j = 0; j<ncol; j++)
      *p++ /= gcd;
      Denom(tp,k) = Denom(tp,k)/gcd;
   }
   #endif
 }
 p = tp->row[pivi].objet.val;
 for(k = 0; k<nligne; k++)
   if((Flag(tp, k) & Unit) && tp->row[k].objet.unit == pivj) break;
 Flag(tp, k) = Plus;
 tp->row[k].objet.val = p;
 for(j = 0; j<ncol; j++)
   #if defined(LINEAR_VALUE_IS_MP)
   mpz_set(*p++, new[j]);
   #else
   *p++ = new[j];
   #endif

 #if defined(LINEAR_VALUE_IS_MP)
 mpz_set(Denom(tp, k), pivot);
 Flag(tp, pivi) = Unit | Zero;
 mpz_set(Denom(tp, pivi), UN);
 #else
 Denom(tp, k) = pivot; 
 Flag(tp, pivi) = Unit | Zero;
 Denom(tp, pivi) = UN;
 #endif
 tp->row[pivi].objet.unit = pivj;

 for(k = 0; k<nligne; k++){
   ff = Flag(tp, k);
   if(ff & Unit) continue;
   #if defined(LINEAR_VALUE_IS_MP)
   sgn_x = mpz_sgn(Index(tp, k, pivj));
   #else
   sgn_x = Index(tp, k, pivj);
   #endif
   if(sgn_x < 0) fff = Minus;
   else if(sgn_x == 0) fff = Zero;
   else fff = Plus;
   if(fff != Zero && fff != ff)
     if(ff == Zero) ff = (fff == Minus ? Unknown : fff);
     else ff = Unknown;
   Flag(tp, k) = ff;
 }

 if(verbose>2){
   fprintf(dump, "just pivoted\n");
   tab_display(tp, dump);
 }

 #if defined(LINEAR_VALUE_IS_MP)
 mpz_clear(x); mpz_clear(y); mpz_clear(d); mpz_clear(gcd);
 mpz_clear(u); mpz_clear(dpiv); mpz_clear(lpiv);
 mpz_clear(pivot); mpz_clear(foo); mpz_clear(z);
 mpz_clear(ppivot); mpz_clear(dppiv);

 for(i=0; i<ncol; i++)
   mpz_clear(new[i]);
 #endif

 return(0);
}

/*
 * Sort the rows in increasing order of the largest coefficient
 * and (if TRAITER_DUAL is set) return the new position of the
 * original constraints.
 */
static int *tab_sort_rows(Tableau *tp, int nvar, int nligne, int flags)
{
    int i, j;
    int pivi;
    double s, t, d, smax = 0;
    struct L temp;
    int *pos = NULL, *ineq = NULL;

    if (flags & TRAITER_DUAL) {
	ineq = malloc(tp->height * sizeof(int));
	pos = malloc((tp->height-nvar) * sizeof(int));
	if (!ineq || !pos) {
	    fprintf(stderr, "Memory Overflow.\n") ;
	    exit(1) ;
	}
    }

    for (i = nvar; i < nligne; i++) {
	if (Flag(tp,i) & Unit)
	    continue;
	s = 0;
	d = ENTIER_TO_DOUBLE(Denom(tp, i));
	for (j = 0; j < nvar; j++) {
	    t = ENTIER_TO_DOUBLE(Index(tp,i,j))/d;
	    s = max(s, abs(t));
	}
	tp->row[i].size = s;
	smax = max(s, smax);
	if (flags & TRAITER_DUAL)
	    ineq[i] = i-nvar;
    }

    for (i = nvar; i < nligne; i++) {
	if (Flag(tp,i) & Unit)
	    continue;
	s = smax;
	pivi = i;
	for (j = i; j < nligne; j++) {
	    if (Flag(tp,j) & Unit)
		continue;
	    if (tp->row[j].size < s) {
		s = tp->row[j].size;
		pivi = j;
	    }
	}
	if (pivi != i) {
	    temp = tp->row[pivi];
	    tp->row[pivi] = tp->row[i];
	    tp->row[i] = temp;
	    if (flags & TRAITER_DUAL) {
		j = ineq[i];
		ineq[i] = ineq[pivi];
		ineq[pivi] = j;
	     }
	}
    }

    if (flags & TRAITER_DUAL) {
	for (i = nvar; i < nligne; i++)
	    pos[ineq[i]] = i;
	free(ineq);
    }

    return pos;
}

/* dans cette version, "traiter" modifie ineq; par contre
   le contexte est immediatement recopie' */

void traiter(Tableau *tp, Tableau *ctxt, int nvar, int nparm, int ni, int nc,
	     int bigparm, int flags)
{
 int j;
 int pivi, nligne, ncol;
 struct high_water_mark x;
 Tableau *context;
 int dch, dcw;
 int i;
 int *pos;

 #if !defined(LINEAR_VALUE_IS_MP)
 Entier D = UN;
 #endif

 #if defined(LINEAR_VALUE_IS_MP)
 dcw = mpz_sizeinbase(tp->determinant, 2);
 #else
 dcw = 0;
 for(i=0; i<tp->l_determinant; i++)
   dcw += llog(tp->determinant[i]);
 #endif
 dch = 2 * dcw + 1;
 x = tab_hwm();
 nligne = nvar+ni;

 context = expanser(ctxt, 0, nc, nparm+1, 0, dch, dcw);

 pos = tab_sort_rows(tp, nvar, nligne, flags);

 for(;;) {
   if(verbose>2){
     fprintf(dump, "debut for\n");
     tab_display(tp, dump);
     fflush(dump);
   }
   nligne = nvar+ni; ncol = nvar+nparm+1;
   if(nligne > tp->height || ncol > tp->width) {
     fprintf(stdout, "Syserr : traiter : tableau too small\n");
     exit(1);
   }
   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;	       /* There is a negative row   */
   
   pivi = exam_coef(tp, nvar, ncol, bigparm);

   if(verbose>2){
     fprintf(dump, "coefs examined\n");
     tab_display(tp, dump);
     fflush(dump);
   }

   if(pivi < nligne) goto pirouette;
   /* There is a row whose coefficients are negative */
   compa_test(tp, context, ni, nvar, nparm, nc);
   if(verbose>2){
     fprintf(dump, "compatibility tested\n");
     tab_display(tp, dump);
     fflush(dump);
   }

   pivi = chercher(tp, Minus, nligne);
   if(pivi < nligne) goto pirouette;
   /* The compatibility test has found a negative row */
   pivi = chercher(tp, Critic, nligne);
   if(pivi >= nligne)pivi = chercher(tp, Unknown, nligne);
   /* Here, the problem tree splits        */
   if(pivi < nligne) {
     Tableau * ntp;
     Entier com_dem;
     struct high_water_mark q;
     if(nc >= context->height) {
       #if defined(LINEAR_VALUE_IS_MP)
       dcw = mpz_sizeinbase(context->determinant,2);
       #else
       dcw = 0;
       for(i=0; i<tp->l_determinant; i++)
       dcw += llog(tp->determinant[i]);
       #endif
       dch = 2 * dcw + 1;
       context = expanser(context, 0, nc, nparm+1, 0, dch, dcw);
     }
     if(nparm >= MAXPARM) {
       fprintf(stdout, "Too much parameters : %d\n", nparm);
       exit(2);
     }
     q = tab_hwm();
     if(verbose>1)
       fprintf(stdout,"profondeur %d %lx\n", profondeur, q.top);
     ntp = expanser(tp, nvar, ni, ncol, 0, 0, 0);
     fflush(stdout);
     sol_if();
     sol_forme(nparm+1);
     entier_init_zero(com_dem);
     for (j = 0; j < nparm; j++)
       entier_gcd(com_dem, com_dem, Index(tp, pivi, j + nvar +1));
     if (!(flags & TRAITER_INT))
	 entier_gcd(com_dem, com_dem, Index(tp, pivi, nvar));
     for (j = 0; j < nparm; j++) {
       entier_divexact(Index(context, nc, j), Index(tp, pivi, j + nvar + 1), com_dem);
       sol_val(Index(context, nc, j), UN);
     }
     if (!(flags & TRAITER_INT))
	 entier_divexact(Index(context, nc, nparm), Index(tp, pivi, nvar), com_dem);
     else
	 entier_pdivision(Index(context, nc, nparm), Index(tp, pivi, nvar), com_dem);
     sol_val(Index(context, nc, nparm), UN);
     entier_clear(com_dem);
     Flag(context, nc) = Unknown;
     entier_set_si(Denom(context, nc), 1);
     Flag(ntp, pivi) = Plus;
     profondeur++;
     fflush(stdout);
     if(verbose > 0) fflush(dump);
     #if defined(LINEAR_VALUE_IS_MP)
     traiter(ntp, context, nvar, nparm, ni, nc+1, bigparm, flags);
     profondeur--;
     tab_reset(q);
     if(verbose>1)
       fprintf(stdout, "descente %d %lx\n", profondeur, tab_hwm().top);
     for(j = 0; j<nparm; j++)
       mpz_neg(Index(context, nc, j), Index(context, nc, j));
     mpz_add_ui(Index(context, nc, nparm), Index(context, nc, nparm), 1);
     mpz_neg(Index(context, nc, nparm), Index(context, nc, nparm));
     Flag(tp, pivi) = Minus;
     mpz_set(Denom(context, nc), UN);
     #else
     traiter(ntp, context, nvar, nparm, ni, nc+1, bigparm, flags);
     profondeur--;
     tab_reset(q);
     if(verbose>1)
       fprintf(stderr, "descente %d %lx\n", profondeur, tab_hwm().top);
     for(j = 0; j<nparm; j++)
       Index(context, nc, j) = - Index(context, nc, j);
     Index(context, nc, nparm) = - Index(context, nc, nparm) -1;
     Flag(tp, pivi) = Minus;
     Denom(context, nc) = UN;
     #endif
     nc++;
     goto pirouette;
   }
/* Here, all rows are positive. Do we need an integral solution?      */
   if (!(flags & TRAITER_INT)) {
     solution(tp, nvar, nparm);
     if (flags & TRAITER_DUAL)
	solution_dual(tp, nvar, nparm, pos);
     break;
   }
/* Yes we do! */
   pivi = integrer(&tp, &context, &nvar, &nparm, &ni, &nc, bigparm);
   if(pivi > 0) goto pirouette;
		    /* A cut has been inserted and is always negative */
/* Here, either there is an integral solution, */
   if(pivi == 0) solution(tp, nvar, nparm);
/* or no solution exists */
   else sol_nil();
   break;

/* Here, a negative row has been found. The call to <<pivoter>> executes
      a pivoting step                                                 */

pirouette :
     if (pivoter(tp, pivi, nvar, nparm, ni) < 0) {
       sol_nil();
       break;
     }
 }
/* Danger : a premature return would induce memory leaks   */
 tab_reset(x);
 free(pos);
 return;
}
