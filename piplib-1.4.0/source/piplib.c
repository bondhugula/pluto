/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.c                                   *
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
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

/* Premiere version du 30 juillet 2001. */

# include <stdlib.h>
# include <stdio.h>
# include <ctype.h>
#include <string.h>
#ifdef WIN32
#include <windows.h>
#endif

#include "pip.h"
#define min(x,y) ((x) < (y)? (x) : (y))

Entier UN;
Entier ZERO;

long int cross_product, limit;
int allocation, comptage;
int verbose = 0;
int profondeur = 0;
int compa_count;
int deepest_cut = 0;

FILE *dump = NULL;

/* Larger line buffer to accomodate Frédo Vivien exemples. A version
handling arbitrary line length should be written ASAP.
*/

#define INLENGTH 2048

char inbuff[INLENGTH];
int inptr = 256;
int proviso = 0;


/******************************************************************************
 *                 Fonctions d'acquisition de données (ex-maind.c)            *
 ******************************************************************************/


int dgetc(FILE *foo)
{
 char *p;
 if(inptr >= proviso)
   {p = fgets(inbuff, INLENGTH, foo);
    if(p == NULL) return EOF;
    proviso = min(INLENGTH, strlen(inbuff));
    inptr = 0;
    if(verbose > 2) fprintf(dump, "-- %s", inbuff);
  }
 return inbuff[inptr++];
}

#ifdef WIN32
static char *create_temp_filename()
{
    static char dump_name[MAX_PATH];

    GetTempFileName(".", "Pip", 0, dump_name);
    return dump_name;
}
#else
static char *create_temp_filename()
{
    static char dump_name_template[] = "PipXXXXXX";
    static char dump_name[sizeof(dump_name_template)];

    strcpy(dump_name, dump_name_template);
    mkstemp(dump_name);
    return dump_name;
}
#endif

FILE *pip_create_dump_file()
{
    char *g;
    FILE *dump;

    g = getenv("DEBUG");
    if (g && *g) {
    	dump = fopen(g, "w");
        if (!dump)
	    fprintf(stderr, "%s unaccessible\n", g);
    } else
        dump = fopen(create_temp_filename(), "w");
    return dump;
}


#if defined(LINEAR_VALUE_IS_MP)
int dscanf(FILE *foo, Entier  val)
#else
int dscanf(FILE *foo, Entier *val)
#endif
{
 char * p;
 int c;

 for(;inptr < proviso; inptr++)
   if(inbuff[inptr] != ' ' && inbuff[inptr] != '\n' && inbuff[inptr] != '\t')
				break;
 while(inptr >= proviso)
   {p = fgets(inbuff, 256, foo);
    if(p == NULL) return EOF;
    proviso = strlen(inbuff);
    if(verbose > 2) {
      fprintf(dump, ".. %s", inbuff);
      fflush(dump);
    }
    for(inptr = 0; inptr < proviso; inptr++)
       if(inbuff[inptr] != ' '
       && inbuff[inptr] != '\n'
       && inbuff[inptr] != '\t') break;
  }
 #if defined(LINEAR_VALUE_IS_MP)
 if(gmp_sscanf(inbuff+inptr, GMP_INPUT_FORMAT, val) != 1)
 #else
 if(sscanf(inbuff+inptr, FORMAT, val) != 1)
 #endif
 return -1;
 
 for(; inptr < proviso; inptr++)
	if((c = inbuff[inptr]) != '-' && !isdigit(c)) break;
 return 0;
}


/******************************************************************************
 *                    Fonctions d'affichage des structures                    *
 ******************************************************************************/


/* Fonction pip_matrix_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipMatrix qu'elle recoit en parametre.
 * Premiere version : Ced. 29 juillet 2001. 
 */
void pip_matrix_print(FILE * foo, PipMatrix * Mat)
{ Entier * p;
  int i, j ;
  unsigned NbRows, NbColumns ;

  fprintf(foo,"%d %d\n", NbRows=Mat->NbRows, NbColumns=Mat->NbColumns) ;
  for (i=0;i<NbRows;i++) 
  { p=*(Mat->p+i) ;
    for (j=0;j<NbColumns;j++)
    #if defined(LINEAR_VALUE_IS_MP)
    { fprintf(foo," ") ;
      mpz_out_str(foo,10,*p++) ;
    }
    #else
    fprintf(foo," "FORMAT, *p++) ;
    #endif
    fprintf(foo, "\n") ;
  }
} 


/* Fonction pip_vector_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipVector qu'elle recoit en parametre.
 * Premiere version : Ced. 20 juillet 2001. 
 */
void pip_vector_print(FILE * foo, PipVector * vector)
{ int i ;
  
  if (vector != NULL)
  { fprintf(foo,"#[") ;
    for (i=0;i<vector->nb_elements;i++)
    { fprintf(foo," ") ;
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_out_str(foo,10,vector->the_vector[i]) ;
      if (mpz_cmp(vector->the_deno[i],UN) != 0)
      #else
      fprintf(foo,FORMAT,vector->the_vector[i]) ;
      if (vector->the_deno[i] != UN)
      #endif
      { fprintf(foo,"/") ;
        #if defined(LINEAR_VALUE_IS_MP)
        mpz_out_str(foo,10,vector->the_deno[i]) ;
        #else
        fprintf(foo,FORMAT,vector->the_deno[i]) ;
        #endif
      }
    }
    fprintf(foo,"]") ;
  }
}
  

/* Fonction pip_newparm_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipNewparm qu'elle recoit en parametre.
 * Le parametre indent est le nombre d'espaces blancs en debut de chaque
 * ligne avant indentation. Une valeur negative de indent signifie qu'on ne
 * desire pas d'indentation.
 * Premiere version : Ced. 18 octobre 2001. 
 */
void pip_newparm_print(FILE * foo, PipNewparm * newparm, int indent)
{ int i ;

  if (newparm != NULL)
  { do
    { for (i=0;i<indent;i++) fprintf(foo," ") ;             /* Indent. */
      fprintf(foo,"(newparm ") ;
      fprintf(foo,"%d",newparm->rank) ;
      fprintf(foo," (div ") ;
      pip_vector_print(foo,newparm->vector) ;
      fprintf(foo," ") ;
      #if defined(LINEAR_VALUE_IS_MP)
      mpz_out_str(foo,10,newparm->deno) ;
      #else
      fprintf(foo,FORMAT,newparm->deno) ;
      #endif
      fprintf(foo,"))\n") ;
    }
    while ((newparm = newparm->next) != NULL) ;
  }  
}


/* Fonction pip_list_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipList qu'elle recoit en parametre.
 * Le parametre indent est le nombre d'espaces blancs en debut de chaque
 * ligne avant indentation. Une valeur negative de indent signifie qu'on ne
 * desire pas d'indentation.
 * Premiere version : Ced. 18 octobre 2001. 
 * 16 novembre 2005 : Ced. Prise en compte du cas list->vector == NULL,
 *                         jusque là impossible.
 */
void pip_list_print(FILE * foo, PipList * list, int indent)
{ int i ;

  if (list == NULL)
  { for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,"()\n") ;
  }
  else
  { for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,"(list\n") ;
    do
    { if (list->vector != NULL)
      { for (i=0;i<indent+1;i++) fprintf(foo," ") ;         /* Indent. */
        pip_vector_print(foo,list->vector) ;
        fprintf(foo,"\n") ;
      }
    }
    while ((list = list->next) != NULL) ;
    for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,")\n") ;
  }
}


/* Fonction pip_quast_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipQuast qu'elle recoit en parametre.
 * Le parametre indent est le nombre d'espaces blancs en debut de chaque
 * ligne avant indentation. Une valeur negative de indent signifie qu'on ne
 * desire pas d'indentation.
 * 20 juillet 2001 : Premiere version, Ced. 
 * 18 octobre 2001 : eclatement. 
 */
void pip_quast_print(FILE * foo, PipQuast * solution, int indent)
{ int i ;
  PipVector * vector ;
  int new_indent = indent >= 0 ? indent+1 : indent;
  
  if (solution != NULL)
  { pip_newparm_print(foo,solution->newparm,indent) ;
    if (solution->condition == NULL) {
	pip_list_print(foo, solution->list, indent);
	/* Possible dual solution */
	if (solution->next_then)
	    pip_quast_print(foo, solution->next_then, new_indent);
    } else
    { for (i=0;i<indent;i++) fprintf(foo," ") ;             /* Indent. */
      fprintf(foo,"(if ") ;
      pip_vector_print(foo,solution->condition) ;
      fprintf(foo,"\n") ;
      pip_quast_print(foo, solution->next_then, new_indent);
      pip_quast_print(foo, solution->next_else, new_indent);
      for (i=0;i<indent;i++) fprintf(foo," ") ;             /* Indent. */
      fprintf(foo,")\n") ;
    }
  }
  else
  { for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,"void\n") ;
  }
}    
  

/* Function pip_options_print:
 * This function prints the content of a PipOptions structure (options)
 * into a file (foo, possibly stdout).
 * March 17th 2003: first version.
 */
void * pip_options_print(FILE * foo, PipOptions * options)
{ fprintf(foo,"Option setting is:\n") ;
  fprintf(foo,"Nq          =%d\n",options->Nq) ;
  fprintf(foo,"Verbose     =%d\n",options->Verbose) ;
  fprintf(foo,"Simplify    =%d\n",options->Simplify) ;
  fprintf(foo,"Deepest_cut =%d\n",options->Deepest_cut) ;
  fprintf(foo,"Maximize    =%d\n",options->Maximize) ;
  fprintf(foo,"Urs_parms   =%d\n",options->Urs_parms);
  fprintf(foo,"Urs_unknowns=%d\n",options->Urs_unknowns);
  fprintf(foo,"\n") ;
}


/******************************************************************************
 *                       Fonctions de liberation memoire                      *
 ******************************************************************************/


/* Fonction pip_matrix_free :
 * Cette fonction libere la memoire reservee a la structure de type PipMatrix
 * que pointe son parametre.
 * Premiere version : Ced. 29 juillet 2001. 
 */
void pip_matrix_free(PipMatrix * matrix)
{ 
  #if defined(LINEAR_VALUE_IS_MP)
  int i, j ;
  Entier * p ;

  p = matrix->p_Init ;
  for (i=0;i<matrix->p_Init_size;i++) 
  mpz_clear(*p++) ;
  #endif

  if (matrix != NULL)
  { free(matrix->p_Init) ;
    free(matrix->p) ;
    free(matrix) ;
  }
}


/* Fonction pip_vector_free :
 * Cette fonction libere la memoire reservee a la structure de type PipVector
 * que pointe son parametre.
 * 20 juillet 2001 : Premiere version, Ced.
 * 18 octobre 2001 : simplification suite a l'eclatement de PipVector.
 * 16 novembre 2005 : Ced. Prise en compte du cas NULL.
 */
void pip_vector_free(PipVector * vector)
{ int i ;
  
  if (vector != NULL)
  { 
    #if defined(LINEAR_VALUE_IS_MP)
    for (i=0;i<vector->nb_elements;i++)
    { mpz_clear(vector->the_vector[i]);
      mpz_clear(vector->the_deno[i]);
    }
    #endif
  
    free(vector->the_vector) ;
    free(vector->the_deno) ;
    free(vector) ;
  }
}


/* Fonction pip_newparm_free :
 * Cette fonction libere la memoire reservee a la structure de type PipNewparm
 * que pointe son parametre. Sont liberes aussi tous les elements de la
 * liste chainee dont il pouvait etre le depart.
 * Premiere version : Ced. 18 octobre 2001. 
 */
void pip_newparm_free(PipNewparm * newparm)
{ PipNewparm * next ;

  while (newparm != NULL)
  { next = newparm->next ;
    #if defined(LINEAR_VALUE_IS_MP)
    mpz_clear(newparm->deno);
    #endif
    pip_vector_free(newparm->vector) ;
    free(newparm) ;
    newparm = next ;
  }
}


/* Fonction pip_list_free :
 * Cette fonction libere la memoire reservee a la structure de type PipList
 * que pointe son parametre. Sont liberes aussi tous les elements de la
 * liste chainee dont il pouvait etre le depart.
 * Premiere version : Ced. 18 octobre 2001. 
 */
void pip_list_free(PipList * list)
{ PipList * next ;

  while (list != NULL)
  { next = list->next ;
    pip_vector_free(list->vector) ;
    free(list) ;
    list = next ;
  }
}


/* Fonction pip_quast_free :
 * Cette fonction libere la memoire reservee a la structure de type 
 * PipSolution que pointe son parametre. Sont liberees aussi toutes les
 * differentes listes chainees qui pouvaient en partir.
 * 20 juillet 2001 : Premiere version, Ced.
 * 18 octobre 2001 : simplification suite a l'eclatement de PipVector.
 */
void pip_quast_free(PipQuast * solution)
{
    if (!solution)
	return;
    pip_newparm_free(solution->newparm) ;
  
    pip_list_free(solution->list) ;

    pip_vector_free(solution->condition);
    pip_quast_free(solution->next_then);
    pip_quast_free(solution->next_else);
    free(solution) ;
}


/* Funtion pip_options_free:
 * This function frees the allocated memory for a PipOptions structure.
 * March 15th 2003: first version.
 */
void pip_options_free(PipOptions * options)
{ free(options) ;
}


/******************************************************************************
 *                     Fonction d'initialisation des options                  *
 ******************************************************************************/


/* Funtion pip_options_init:
 * This function allocates the memory for a PipOptions structure and fill the
 * options with the default values.
 ********
 * Nq est un booleen renseignant si on cherche
 * une solution entiere (vrai=1) ou non (faux=0). Verbose est un booleen
 * permettant de rendre Pip bavard (Verbose a vrai=1), il imprimera
 * alors la plupart de ses traitements dans le fichier dont le nom est
 * dans la variable d'environnement DEBUG, ou si DEBUG
 * n'est pas placee, dans un nouveau fichier de nom genere par mkstemp, si
 * Verbose est a faux=0, Pip restera muet. Simplify est un booleen permettant
 * de demander a Pip de simplifier sa solution (en eliminant les formes de type
 * 'if #[...] () ()') quand il est a vrai=1, ou non quand il est a faux=0. Max
 * n'est pas encore utilise et doit etre mis a 0. 
 ********
 * March 15th 2003: first version.
 */ 
PipOptions * pip_options_init(void)
{ PipOptions * options ;

  options = (PipOptions *)malloc(sizeof(PipOptions)) ;
  /* Default values of the options. */
  options->Nq          = 1 ;  /* Integer solution. */
  options->Verbose     = 0 ;  /* No comments. */
  options->Simplify    = 0 ;  /* Do not simplify solutions. */
  options->Deepest_cut = 0 ;  /* Do not use deepest cut algorithm. */
  options->Maximize    = 0 ;  /* Do not compute maximum. */
  options->Urs_parms   = 0 ;  /* All parameters are non-negative. */
  options->Urs_unknowns= 0 ;  /* All unknows are non-negative. */
  options->Compute_dual= 0;   /* Don't compute dual variables. */
  
  return options ;
}


/******************************************************************************
 *                     Fonctions d'acquisition de matrices                    *
 ******************************************************************************/


/* Fonction pip_matrix_alloc :
 * Fonction (tres) modifiee de Matrix_Alloc de la polylib. Elle alloue l'espace
 * memoire necessaire pour recevoir le contenu d'une matrice de NbRows lignes
 * et de NbColumns colonnes, et initialise les valeurs a 0. Elle retourne un
 * pointeur sur l'espace memoire alloue.
 * Premiere version : Ced. 18 octobre 2001. 
 */
PipMatrix * pip_matrix_alloc(unsigned NbRows, unsigned NbColumns)
{ PipMatrix * matrix ;
  Entier ** p, * q ;
  int i, j ;

  matrix = (PipMatrix *)malloc(sizeof(PipMatrix)) ;
  if (matrix == NULL) 	
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  matrix->NbRows = NbRows ;
  matrix->NbColumns = NbColumns ;
  matrix->p_Init_size = NbRows * NbColumns ;
  if (NbRows == 0) 
  { matrix->p = NULL ;
    matrix->p_Init = NULL ;
  }  
  else 
  { if (NbColumns == 0) 
    { matrix->p = NULL ;
      matrix->p_Init = NULL ;
    }
    else 
    { p = (Entier **)malloc(NbRows*sizeof(Entier *)) ;
      if (p == NULL) 
      { fprintf(stderr, "Memory Overflow.\n") ;
        exit(1) ;
      }
      q = (Entier *)malloc(NbRows * NbColumns * sizeof(Entier)) ;
      if (q == NULL) 
      { fprintf(stderr, "Memory Overflow.\n") ;
        exit(1) ;
      }
      matrix->p = p ;
      matrix->p_Init = q ;
      for (i=0;i<NbRows;i++) 
      { *p++ = q ;
	for (j=0;j<NbColumns;j++)   
        #if defined(LINEAR_VALUE_IS_MP)
	mpz_init_set_si(*(q+j),0) ;
	#else
	*(q+j) = 0 ;
	#endif
	q += NbColumns ;
      }
    }
  }
  return matrix ;
}


/* Fonction pip_matrix_read :
 * Adaptation de Matrix_Read de la polylib. Cette fonction lit les donnees
 * d'une matrice dans un fichier 'foo' et retourne un pointeur vers une
 * structure ou elle a copie les informations de cette matrice. Les donnees
 * doivent se presenter tq :
 * - des lignes de commentaires commencant par # (optionnelles),
 * - le nombre de lignes suivit du nombre de colonnes de la matrice, puis
 *   eventuellement d'un commentaire sur une meme ligne,
 * - des lignes de la matrice, chaque ligne devant etre sur sa propre ligne de
 *   texte et eventuellement suivies d'un commentaire.
 * Premiere version : Ced. 18 octobre 2001. 
 * 24 octobre 2002 : premiere version MP, attention, uniquement capable de
 *                   lire des long long pour l'instant. On utilise pas
 *                   mpz_inp_str car on lit depuis des char * et non des FILE.
 */
PipMatrix * pip_matrix_read(FILE * foo)
{ unsigned NbRows, NbColumns ;
  int i, j, n ;
  #if defined(LINEAR_VALUE_IS_MP)
  long long val ;
  #endif
  char *c, s[1024], str[1024] ;
  PipMatrix * matrix ;
  Entier * p ;

  while (fgets(s,1024,foo) == 0) ;
  while ((*s=='#' || *s=='\n') || (sscanf(s," %u %u",&NbRows,&NbColumns)<2))
  fgets(s, 1024, foo) ;
  
  matrix = pip_matrix_alloc(NbRows,NbColumns) ;

  p = matrix->p_Init ;
  for (i=0;i<matrix->NbRows;i++) 
  { do 
    { c = fgets(s,1024,foo) ;
      while (isspace(*c) && (*c != '\n'))
      c++ ;
    }
    while (c != NULL && (*c == '#' || *c == '\n'));
    
    if (c == NULL) 
    { fprintf(stderr, "Not enough rows.\n") ;
      exit(1) ;
    }
    for (j=0;j<matrix->NbColumns;j++) 
    { if (c == NULL || *c == '#' || *c == '\n')
      { fprintf(stderr, "Not enough columns.\n") ;
        exit(1) ;
      }
      /* NdCed : Dans le n ca met strlen(str). */
      if (sscanf(c,"%s%n",str,&n) == 0) 
      { fprintf(stderr, "Not enough rows.\n") ;
        exit(1) ;
      }
      #if defined(LINEAR_VALUE_IS_MP)
      sscanf(str,"%lld",&val) ;
      mpz_set_si(*p++,val) ;
      #else
      sscanf(str,FORMAT,p++) ;
      #endif
      c += n ;
    }
  }
  return matrix ;
}
 
 
/* initialization of pip */
static int pip_initialized = 0;

void pip_init() {
  /* Avoid initializing (and leaking) several times */
  if (!pip_initialized) {
    #if defined(LINEAR_VALUE_IS_MP)
    mpz_init_set_si(UN, 1);
    mpz_init_set_si(ZERO, 0);
    #else
    UN   = VAL_UN ;
    ZERO = VAL_ZERO ;
    #endif
    sol_init() ;
    tab_init() ;
    pip_initialized = 1;
  }
}

void pip_close() {
  tab_close();
  sol_close();
# if defined(LINEAR_VALUE_IS_MP)
  mpz_clear(UN);
  mpz_clear(ZERO);
# endif
  pip_initialized = 0;
}
 
/*
 * Compute dual variables of equalities from dual variables of
 * the corresponding pair of inequalities.
 *
 * In practice, this means we need to remove one of the two
 * dual variables that corrsponding to the inequalities
 * and negate the value if it corresponds to the negative
 * inequality.
 */
static void pip_quast_equalities_dual(PipQuast *solution, PipMatrix *inequnk)
{
    PipList **list_p, *list;
    int i;

    if (!solution)
	return;
    if (solution->condition) {
	pip_quast_equalities_dual(solution->next_then, inequnk);
	pip_quast_equalities_dual(solution->next_else, inequnk);
    }
    if (!solution->list)
	return;
    if (!solution->next_then)
	return;
    if (!solution->next_then->list)
	return;
    list_p = &solution->next_then->list;
    for (i = 0; i < inequnk->NbRows; ++i) {
	if (entier_zero_p(inequnk->p[i][0])) {
	    if (entier_notzero_p((*list_p)->vector->the_vector[0])) {
		list_p = &(*list_p)->next;
		list = *list_p;
		*list_p = list->next;
		list->next = NULL;
		pip_list_free(list);
	    } else {
		list = *list_p;
		*list_p = list->next;
		list->next = NULL;
		pip_list_free(list);
		entier_oppose((*list_p)->vector->the_vector[0],
			      (*list_p)->vector->the_vector[0]);
		list_p = &(*list_p)->next;
	    }
	} else
	    list_p = &(*list_p)->next;
    }
}


/******************************************************************************
 *                           Fonction de resolution                           *
 ******************************************************************************/


/* Fonction pip_solve :
 * Cette fonction fait un appel a Pip pour la resolution d'un probleme. Le
 * probleme est fourni dans les arguments. Deux matrices de la forme de celles
 * utilisees dans la Polylib forment les systemes d'equations/inequations :
 * un pour les inconnues, l'autre pour les parametres. Bg est le 'bignum'.
 * Le dernier argument contient les options guidant le comportement de PIP.
 * Cette fonction retourne la solution sous la forme d'un arbre de structures
 * PipQuast.
 * 30 juillet 2001 : Premiere version, Ced. 
 * 18 octobre 2001 : suppression de l'argument Np, le nombre de parametres. Il
 *                   est a present deduit de ineqpar. Si ineqpar vaut NULL,
 *                   c'est que Np vaut 0. S'il y a des parametres mais pas de
 *                   contraintes dessus, ineqpar sera une matrice de 0 lignes
 *                   mais du bon nombre de colonnes (Np + 2).
 * 27 février 2003 : Verbose est maintenant gradué.
 *                  -1 -> silence absolu
 *                   0 -> silence relatif
 *                   1 -> information sur les coupures dans le cas ou on
 *                        cherche une solution entière.
 *                   2 -> information sur les pivots et les déterminants
 *                   3 -> information sur les tableaux.
 *                         Chaque option inclut les précédentes. [paf]
 * 15 mars 2003    : passage a une structure d'options.
 */
PipQuast * pip_solve(inequnk, ineqpar, Bg, options)
PipMatrix * inequnk, * ineqpar ;
int Bg ;
PipOptions * options ;
{ Tableau * ineq, * context, * ctxt ;
  int i, Np, Nn, Nl, Nm, p, q, xq, non_vide, Shift = 0, Urs_parms = 0;
  char * g ;
  struct high_water_mark hq ;
  Entier D ;
  PipQuast * solution ;
  int	sol_flags = 0;

  pip_init() ;
   	
  /* initialisations diverses :
   * - la valeur de Verbose et Deepest_cut sont placees dans leurs variables
   *   globales. Dans le cas ou on doit etre en mode verbose, on ouvre le
   *   fichier dans lequel ecrire les tracages. Si la variable d'environnement
   *   DEBUG est placee, on ecrira dans le nom de fichier correspondant, sinon,
   *   dans un nouveau fichier de nom genere par mkstemp,
   * - limit est mis au 0 long int (sa valeur par defaut dans Pip original),
   * - on lance les initialisations pour tab et sol (autres mises en place
   *   de variables globales).
   */
  verbose = options->Verbose ;
  deepest_cut = options->Deepest_cut ;
  if (verbose > 0) {
     dump = pip_create_dump_file();
     if (!dump)
	verbose = 0;
  }
  #if defined(LINEAR_VALUE_IS_MP)
  limit = 0LL ;
  #else
  limit = ZERO ;
  #endif

  /* Si inequnk est NULL, la solution est automatiquement void (NULL). */
  if (inequnk != NULL)
  { /* Np vaut 0 si ineqpar vaut NULL, ineqpar->NbColumns - 2 sinon (on a -1
     * pour la constante et -1 pour le marqueur d'egalite/inegalite = -2).
     */
    Np = (ineqpar == NULL) ? 0 : ineqpar->NbColumns - 2 ;
    /* Calcul du nombre d'inconnues du probleme. Comme les matrices d'entree
     * sont sous la forme de la polylib.
     */
    Nn = inequnk->NbColumns - Np - 2 ;
    /* Calcul du nombre d'inequations du probleme. Le format de matrice de la
     * polylib permet les egalites, on doit donc les compter double quand il
     * y en a.
     */
    Nl = inequnk->NbRows ;
    for (i=0;i<inequnk->NbRows;i++)
    #if defined(LINEAR_VALUE_IS_MP)
    if (mpz_sgn(**(inequnk->p + i)) == 0)
    #else
    if (**(inequnk->p + i) == 0)
    #endif
    Nl ++ ;
      
    /* On prend les 'marques' de debut de traitement. */
    cross_product = 0 ;
    hq = tab_hwm() ;
    xq = p = sol_hwm();
    
    if (options->Maximize) {
      sol_flags |= SOL_MAX;
      Shift = 1;
    } else if (options->Urs_unknowns) {
      sol_flags |= SOL_SHIFT;
      Shift = -1;
    } 

    if (options->Urs_parms) {
      Urs_parms = Np - (Bg >= 0);
      Np += Urs_parms;
    }

    /* Si un maximum est demande, mais sans bignum, on crée le bignum. */
    if (options->Maximize || options->Urs_unknowns) {
      if (Bg < 0) {
	Bg = inequnk->NbColumns - 1 ; /* On choisit sa place. */
	Np ++ ;                       /* On le compte comme parametre. */
	sol_flags |= SOL_REMOVE;      /* On le supprime apres. */
      }
    }
    
    /* On s'assure d'abord que le systeme pour le contexte n'est pas vide
     * avant de commencer le traitement. Si c'est le cas, la solution est
     * void (NULL).
     */
    if (ineqpar != NULL)
    { /* Calcul du nombre d'inequations sur les parametres. Le format de
       * matrice de la polylib permet les egalites, on doit donc les compter
       * double quand il y en a.
       */
      Nm = ineqpar->NbRows ;
      for (i=0;i<ineqpar->NbRows;i++)
      #if defined(LINEAR_VALUE_IS_MP)
      if (mpz_sgn(**(ineqpar->p + i)) == 0)
      #else
      if (**(ineqpar->p + i) == 0)
      #endif
      Nm ++ ;
      
      context = tab_Matrix2Tableau(ineqpar, Nm, Np-Urs_parms, -1,
				   Shift, Bg-Nn-1, Urs_parms);
      if (options->Nq)
	tab_simplify(context, Np);
      
      if (Nm)
      { /* Traduction du format de matrice de la polylib vers celui de
         * traitement de Pip. Puis traitement proprement dit.
         */
	ctxt = expanser(context, Np, Nm, Np+1, Np, 0, 0) ;
        traiter(ctxt, NULL, Np, 0, Nm, 0, -1, TRAITER_INT);
        non_vide = is_not_Nil(p) ;
        sol_reset(p) ;
      }
      else
      non_vide = Pip_True ;
    }
    else
    { Nm = 0 ;
      ineqpar = pip_matrix_alloc(0, 2);
      context = tab_Matrix2Tableau(ineqpar, Nm, Np-Urs_parms, -1,
				   Shift, Bg-Nn-1, Urs_parms);
      pip_matrix_free(ineqpar);
      non_vide = Pip_True ;
    }
        
    if (verbose > 0)
    fprintf(dump, "%d %d %d %d %d %d\n",Nn,Np,Nl,Nm,Bg,options->Nq) ;
    
    /* S'il est possible de trouver une solution, on passe au traitement. */
    if (non_vide) {
      int flags = 0;
      ineq = tab_Matrix2Tableau(inequnk,Nl,Nn,Nn, Shift,Bg, Urs_parms);
      if (options->Nq)
	tab_simplify(ineq, Nn);
  
      compa_count = 0 ;
      if (options->Nq)
	flags |= TRAITER_INT;
      else if (options->Compute_dual) {
	flags |= TRAITER_DUAL;
	sol_flags |= SOL_DUAL;
      }
      traiter(ineq, context, Nn, Np, Nl, Nm, Bg, flags);

      if (options->Simplify)
      sol_simplify(xq) ;
      q = sol_hwm() ;
      /* On traduit la solution du format de solution de Pip vers un arbre
       * de structures de type PipQuast.
       */
      solution = sol_quast_edit(&xq, NULL, Bg-Nn-1, Urs_parms, sol_flags);
      if ((sol_flags & SOL_DUAL) && Nl > inequnk->NbRows)
	  pip_quast_equalities_dual(solution, inequnk);
      sol_reset(p) ;
    }
    else
    return NULL ;
    tab_reset(hq) ;
  }
  else
  return NULL ;
  
  return(solution) ;
}
