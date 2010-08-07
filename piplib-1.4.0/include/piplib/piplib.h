/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.h                                   *
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
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>

/* Premiere version du 18 septembre 2002. */

#if defined(LINEAR_VALUE_IS_LONGLONG) && defined(PIP_WIDTH_MP)
#error
#endif

#if !defined(LINEAR_VALUE_IS_LONGLONG) && !defined(LINEAR_VALUE_IS_INT)
#if !defined(LINEAR_VALUE_IS_MP)
# error Please define LINEAR_VALUE_IS_* or #include polylib32.h or polylib64.h
#endif
#endif

#if defined(LINEAR_VALUE_IS_LONGLONG)

# define Entier   long long
# define FORMAT   "%lld"
# define VAL_UN   1LL
# define VAL_ZERO 0LL

#ifdef VALUE_TO_INT
#undef VALUE_TO_INT
#endif
#define VALUE_TO_INT(val) ((int)(val))
#define ENTIER_TO_DOUBLE(val) ((double)(val))

#elif defined(LINEAR_VALUE_IS_INT) 

# define Entier   long int
# define FORMAT   "%ld"
# define VAL_UN   1L
# define VAL_ZERO 0L

#ifndef VALUE_TO_INT
#define VALUE_TO_INT(val) ((int)(val))
#endif
#define ENTIER_TO_DOUBLE(val) ((double)(val))

#elif defined(LINEAR_VALUE_IS_MP) 

# include <gmp.h>
# define Entier   mpz_t
# define FORMAT   "%d"
# define GMP_INPUT_FORMAT   "%lZd"

#ifndef VALUE_TO_INT
#define VALUE_TO_INT(val) ((int)mpz_get_si(val))
#endif
#define ENTIER_TO_DOUBLE(val) (mpz_get_d(val))

#endif

#if defined(LINEAR_VALUE_IS_MP) 

#define entier_addto(ref,val1,val2)    	(mpz_add((ref),(val1),(val2)))
#define entier_increment(ref,val)	(mpz_add_ui((ref),(val),1))
#define entier_assign(v1,v2)	    	(mpz_set((v1),(v2)))
#define entier_clear(val)       	(mpz_clear((val)))
#define entier_divexact(d,v1,v2)    	(mpz_divexact((d),(v1),(v2)))
#define entier_gcd(g,v1,v2)	    	(mpz_gcd((g),(v1),(v2)))
#define entier_init(val)            	(mpz_init((val)))
#define entier_init_zero(val)		(mpz_init((val)))
#define entier_init_set(v1,v2)	    	(mpz_init_set((v1),(v2)))
#define entier_pdivision(ref,val1,val2)	(mpz_fdiv_q((ref),(val1),(val2)))
#define entier_pmodulus(ref,val1,val2)	(mpz_fdiv_r((ref),(val1),(val2)))
#define entier_oppose(ref,val)       	(mpz_neg((ref),(val)))
#define entier_set_si(val,i)    	(mpz_set_si((val),(i)))    
#define entier_subtract(ref,val1,val2) 	(mpz_sub((ref),(val1),(val2)))
#define entier_decrement(ref,val)	(mpz_sub_ui((ref),(val),1))
#define entier_sgn(val)			(mpz_sgn(val))
#define entier_eq(v1,v2) 	    	(mpz_cmp((v1),(v2)) == 0)
#define entier_ne(v1,v2) 	    	(mpz_cmp((v1),(v2)) != 0)
#define entier_one_p(val)		(mpz_cmp_si(val,1) == 0)
#define entier_llog(val)		(mpz_sizeinbase(val, 2))

#else

#define entier_addto(ref,val1,val2) 	((ref) = (val1)+(val2))
#define entier_increment(ref,val)	((ref) = (val)+1)
#define entier_assign(v1,v2)	    	((v1) = (v2))
#define entier_clear(val)             	do { } while(0)
#define entier_divexact(d,v1,v2)    	((d) = (v1) / (v2))
#define entier_gcd(g,v1,v2)	    	((g) = pgcd((v1),(v2)))
#define entier_init(val)             	do { } while(0)
#define entier_init_zero(v1)	    	((v1) = 0)
#define entier_init_set(v1,v2)	    	((v1) = (v2))
#define entier_pdivision(ref,val1,val2)	((ref) = ((val1)-mod((val1),(val2)))/(val2))
#define entier_pmodulus(ref,val1,val2)	((ref) = mod((val1),(val2)))
#define entier_oppose(ref,val)    	((ref) = -(val))
#define entier_set_si(val,i)        	((val) = (Entier)(i))   
#define entier_subtract(ref,val1,val2) 	((ref) = (val1)-(val2))
#define entier_decrement(ref,val)	((ref) = (val)-1)
#define entier_sgn(val)			(val)
#define entier_eq(v1,v2) 	    	((v1) == (v2))
#define entier_ne(v1,v2) 	    	((v1) != (v2))
#define entier_one_p(val)        	((val) == 1)
#define entier_llog(val)		(llog(val))

#endif

#define entier_zero_p(val)        	(entier_sgn(val) == 0)
#define entier_notzero_p(val)        	(entier_sgn(val) != 0)
#define entier_pos_p(val)        	(entier_sgn(val) > 0)

#ifndef PIPLIB_H
#define PIPLIB_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 


/* Structure PipMatrix :
 * Structure de matrice au format PolyLib. Le premier element d'une ligne
 * indique quand il vaut 1 que la ligne decrit une inequation de la forme
 * p(x)>=0 et quand il vaut 0, que la ligne decrit une egalite de la forme
 * p(x)=0. Le dernier element de chaque ligne correspond au coefficient
 * constant.
 */
struct pipmatrix
{ unsigned NbRows, NbColumns ;
  Entier **p ;
  Entier *p_Init ;
  int p_Init_size;	        /* Only for PolyLib compatibility under MP
                                 * version: PolyLib makes sometimes
				 * overestimates on the size of the matrices,
				 * in order to go faster. Thus
				 * NbRows*NbColumns is not the number of
				 * allocated elements. With MP version, we
				 * have to think to mpz_clear() all the
				 * initialized elements before freing, then
				 * we need to know the number of allocated
				 * elements: p_Init_size.
				 */
} ;
typedef struct pipmatrix PipMatrix ;


/* Structure PipVector :
 * Cette structure contient un Vector de 'nb_elements' la ieme composante de
 * ce vecteur vaut the_vector[i]/the_deno[i].
 */
struct pipvector
{ int nb_elements ;             /* Nombre d'elements du vecteur. */
  Entier * the_vector ;         /* Numerateurs du vecteur. */
  Entier * the_deno ;           /* Denominateurs du vecteur. */
} ;
typedef struct pipvector PipVector ;


/* Structure PipNewparm :
 * Liste chainee de Newparm, les informations d'un newparm etant son rang, un
 * vecteur de coefficients et un denominateur. Le newparm est egal a la division
 * du vecteur par le denominateur.
 */
struct pipnewparm
{ int rank ;                    /* Rang du 'newparm'. */
  PipVector * vector ;          /* Le vector decrivant le newparm. */
  Entier deno ;                 /* Denominateur du 'newparm'. */
  struct pipnewparm * next ;    /* Pointeur vers le newparm suivant. */
} ;
typedef struct pipnewparm PipNewparm ;


/* Structure PipList :
 * Liste chainee de Vector.
 */
struct piplist
{ PipVector * vector ;          /* Le vector contenant la partie de solution. */
  struct piplist * next ;       /* Pointeur vers l'element suivant. */
} ;
typedef struct piplist PipList ;


/* Structure pipquast :
 * Arbre binaire. Conformement a la grammaire de sortie (voir mode d'emploi), un
 * noeud de l'arbre des solutions debute par une liste de 'newparm'. Il continue
 * ensuite soit par une 'list' (alors condition vaut null), soit par un 'if'
 * (alors le champ condition contient la condition).
 */
struct pipquast
{ PipNewparm * newparm ;        /* Les 'newparm'. */
  PipList * list ;              /* La 'list' si pas de 'if'. */
  PipVector * condition ;       /* La condition si 'if'. */
  struct pipquast * next_then ; /* Noeud si condition et si verifiee. */
  struct pipquast * next_else ; /* Noeud si condition et si non verifiee. */
  struct pipquast * father ;    /* Pointeur vers le quast pere. */
} ;      
typedef struct pipquast PipQuast ;


/* Structure pipoptions:
 * This structure contains each option that can be set to change the PIP
 * behaviour.
 */
struct pipoptions
{ int Nq ;                      /* 1 if an integer solution is needed,
                                 * 0 otherwise.
				 */
  int Verbose ;                 /* -1 -> absolute silence,
                                 *  0 -> relative silence,
                                 *  1 -> information on cuts when an integer
				 *       solution is needed,
                                 *  2 -> information sur les pivots et les
				 *       déterminants,
                                 *  3 -> information on arrays,
                                 * Each option include the preceding.
				 */
  int Simplify ;                /* Set to 1 to eliminate some trivial
                                 * solutions, 0 otherwise.
				 */
  int Deepest_cut ;             /* Set to 1 to include deepest cut
                                 * algorithm.
				 */
  int Maximize;                 /* Set to 1 if maximum is needed. */
  int Urs_parms;             	/* -1 -> all parameters may be negative 
				 *  0 -> all parameters are non-negative
				 */
  int Urs_unknowns;             /* -1 -> all unknowns may be negative 
				 *  0 -> all unknowns are non-negative
				 */
  int Compute_dual;
} ;      
typedef struct pipoptions PipOptions ;


/* Prototypes des fonctions d'affichages des structures de la PipLib. */
void pip_matrix_print(FILE *, PipMatrix *) ;
void pip_vector_print(FILE *, PipVector *) ;
void pip_newparm_print(FILE * foo, PipNewparm *, int indent) ;
void pip_list_print(FILE * foo, PipList *, int indent) ;
void pip_quast_print(FILE *, PipQuast *, int) ;


/* Prototypes des fonctions de liberation memoire des structures de la PipLib.*/
void pip_matrix_free(PipMatrix *) ;
void pip_vector_free(PipVector *) ;
void pip_newparm_free(PipNewparm *) ;
void pip_list_free(PipList *) ;
void pip_quast_free(PipQuast *) ;
void pip_options_free(PipOptions *) ;


/* Prototypes des fonctions d'acquisition de matrices de contraintes et
 * options.
 */
PipMatrix * pip_matrix_alloc(unsigned, unsigned) ;
PipMatrix * pip_matrix_read(FILE *) ;
PipOptions * pip_options_init(void) ;
 

/* initialization of pip library */
void pip_init();
void pip_close();


/* Prototype de la fonction de resolution :
 * pip_solve resoud le probleme qu'on lui passe en parametre, suivant les
 * options elles aussi en parametre. Elle renvoie la solution sous forme
 * d'un arbre de PipQuast. Parametres :
 * - probleme :
 * 1 PipMatrix  : systeme des inequations definissant le domaine des inconnues,
 * 2 PipMatrix  : systeme des inequations satisfaites par les parametres,
 * 3 int        : column rank of the bignum, or negative value if there
 *                is no big parameter.
 * 4 PipOptions : options for PIP.
 */ 
PipQuast * pip_solve(PipMatrix *, PipMatrix *, int, PipOptions *) ;

#define SOL_SHIFT		(1 << 0)    /* Shift solution over -bigparam */
#define SOL_NEGATE		(1 << 1)    /* Negate solution */
#define SOL_REMOVE		(1 << 2)    /* Remove big parameter */
#define SOL_MAX			(SOL_SHIFT | SOL_NEGATE)
					    /* Maximum was computed */
#define SOL_DUAL		(1 << 3)    /* Create dual leaf */
PipQuast *sol_quast_edit(int *i, PipQuast *father, int Bg, int Urs_p, int flags);

#if defined(__cplusplus)
  }
#endif 
#endif /* define PIPLIB_H */
