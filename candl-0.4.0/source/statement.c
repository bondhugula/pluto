
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                  statement.c                            **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: september 9th 2003               **
    **--- |"-.-"| -------------------------------------------------------**
          |     | 
          |     | 
 ******** |     | *************************************************************
 * CAnDL  '-._,-' the Chunky Analyzer for Dependences in Loops (experimental) *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2003 Cedric Bastoul                                          *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * CAnDL, the Chunky Dependence Analyzer                                      *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/
/* CAUTION: the english used for comments is probably the worst you ever read,
 *          please feel free to correct and improve it !
 */

# include <stdlib.h>
# include <stdio.h>
# include <ctype.h>
# include <string.h>
# include "../include/candl/candl.h"


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * candl_statement_print_structure function:
 * Displays a CandlStatement structure (statement) into a file (file,
 * possibly stdout) in a way that trends to be understandable without falling
 * in a deep depression or, for the lucky ones, getting a headache... It
 * includes an indentation level (level) in order to work with others
 * print_structure functions.
 * - 18/09/2003: first version.
 */
void candl_statement_print_structure(file, statement, level)
FILE * file ;
CandlStatement * statement ;
int level ;
{ int i, j ;

  if (statement != NULL)
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- CandlStatement\n") ;
    
    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Go to the right level and print the label. */
    for(j=0; j<=level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"Label: %d\n",statement->label) ;
  
    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Go to the right level and print the type. */
    for(j=0; j<=level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"Type: ") ;
    switch (statement->type)
    { case CANDL_UNSET        : fprintf(file,"UNSET\n") ;           break ;
      case CANDL_ASSIGNMENT   : fprintf(file,"assignment\n") ;     break ;
      case CANDL_P_REDUCTION  : fprintf(file,"plus-reduction\n") ;  break ;
      case CANDL_M_REDUCTION  : fprintf(file,"minus-reduction\n") ; break ;
      case CANDL_T_REDUCTION  : fprintf(file,"times-reduction\n") ; break ;
      default : fprintf(file,"unknown\n") ;
    }
    
    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Go to the right level and print the outer loops. */
    for(j=0; j<=level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"Depth: %d, outer loops(s):",statement->depth) ;
    for (i=0;i<statement->depth;i++)
    fprintf(file," %d",statement->index[i]) ;
    fprintf(file,"\n") ;
    
    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the iteration domain. */
    candl_matrix_print_structure(file,statement->domain,level+1) ;
    
    /* Print the written data. */
    candl_matrix_print_structure(file,statement->written,level+1) ;
    
    /* Print the read data. */
    candl_matrix_print_structure(file,statement->read,level+1) ;
  }
  else
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- NULL statement\n") ;
  }
   
  /* The last line. */
  for(j=0; j<=level; j++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
}


/** 
 * candl_statement_print function:
 * This function prints the content of a CandlStatement structure (statement)
 * into a file (file, possibly stdout).
 * - 09/09/2003: first version.
 */
void candl_statement_print(FILE * file, CandlStatement * statement)
{ candl_statement_print_structure(file,statement,0) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * candl_statement_free function:
 * This function frees the allocated memory for a CandlStatement structure.
 * - 09/09/2003: first version.
 */
void candl_statement_free(CandlStatement * statement)
{ free(statement->index) ;
  pip_matrix_free(statement->domain) ;
  pip_matrix_free(statement->read) ;
  pip_matrix_free(statement->written) ;
  free(statement) ;
}


/******************************************************************************
 *                              Reading functions                             *
 ******************************************************************************/


/**
 * candl_statement_read function:
 * This function reads statement data from a file (file) and puts them into
 * a CandlStatement structure. This function returns a pointer to this
 * structure.
 * - label is the statement number ;
 * - nb_parameters is the number of parameters.
 ***
 * - 09/09/2003: first version.
 */
CandlStatement * candl_statement_read(FILE * file, int label, int nb_parameters)
{ int i, n, * index ;
  char s[CANDL_MAX_STRING], str[CANDL_MAX_STRING], * c, type ;
  CandlStatement * statement ;
  
  /* Statement data. */
  statement = candl_statement_malloc() ;

  /* We first read the statement type. */
  while (fgets(s,CANDL_MAX_STRING,file) == 0) ;
  while ((*s=='#' || *s=='\n') || (sscanf(s," %c",&type)<1))
  fgets(s,CANDL_MAX_STRING,file) ;
  
  switch (type)
  { case 'A': statement->type = CANDL_ASSIGNMENT ; break ;
    case 'P': statement->type = CANDL_P_REDUCTION ; break ;
    case 'M': statement->type = CANDL_M_REDUCTION ; break ;
    case 'T': statement->type = CANDL_T_REDUCTION ; break ;
    default : fprintf(stderr, "[Candl]ERROR: unknown statement type %c\n",type);
              fprintf(stderr, "              possible types are:\n") ;
              fprintf(stderr, "              - A for assignment       (=),\n") ;
              fprintf(stderr, "              - P for plus-reduction  (+=),\n") ;
              fprintf(stderr, "              - M for minus-reduction (-=),\n") ;
              fprintf(stderr, "              - T for times-reduction (*=).\n") ;
              exit(1) ;
  }

  statement->label = label ;
  statement->domain  = pip_matrix_read(file) ;
  statement->depth = statement->domain->NbColumns - nb_parameters - 2 ;

  index = (int *)malloc(sizeof(int)*statement->depth) ;
  if (index == NULL) 
  { fprintf(stderr, "[Candl]ERROR: memory overflow.\n") ;
    exit(1) ;
  }

  do  /* Skip the comments, spaces and empty lines... */
  { c = fgets(s,CANDL_MAX_STRING,file) ;
    while ((c != NULL) && isspace(*c) && (*c != '\n'))
    c++ ;
  }
  while (c != NULL && (*c == '#' || *c == '\n'));
    
  if (c == NULL) 
  { fprintf(stderr, "[Candl]ERROR: no labels in input file.\n") ;
    exit(1) ;
  }
  for (i=0;i<statement->depth;i++) 
  { /* All iterator labels must be on the same line. */
    while (isspace(*c))
    c++ ;
    if (c == NULL || *c == '#' || *c == '\n')
    { fprintf(stderr, "[Candl]ERROR: not enough labels in input file.\n") ;
      exit(1) ;
    }
    /* n is strlen(str). */
    if (sscanf(c,"%s%n",str,&n) == 0) 
    { fprintf(stderr, "[Candl]ERROR: no labels in input file.\n") ;
      exit(1) ;
    }
    sscanf(str,"%d",&index[i]) ;
    c += n ;
  }
  statement->index   = index ;
  
  statement->written = pip_matrix_read(file) ;
  statement->read    = pip_matrix_read(file) ;

  return statement ;
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * candl_statement_malloc function:
 * This function allocates the memory space for a CandlStatement structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - 09/12/2005: first version.
 */
CandlStatement * candl_statement_malloc()
{ CandlStatement * statement ;

  /* Memory allocation for the CloogProgram structure. */
  statement = (CandlStatement *)malloc(sizeof(CandlStatement)) ;
  if (statement == NULL) 
  { fprintf(stderr, "[Candl]ERROR: memory overflow.\n") ;
    exit(1) ;
  }
  
  /* We set the various fields with default values. */
  statement->label   = CANDL_UNSET ;
  statement->type    = CANDL_UNSET ;
  statement->depth   = CANDL_UNSET ;
  statement->index   = NULL ;
  statement->domain  = NULL ;
  statement->written = NULL ;
  statement->read    = NULL ;
  statement->ref     = NULL ;

  return statement ;
}


/**
 * candl_statement_commute function:
 * This function returns 1 if the two statements given as parameter commute,
 * 0 otherwise. It uses the statement type information to answer the question.
 * - 09/12/2005: first version.
 */
int candl_statement_commute(statement1, statement2)
CandlStatement * statement1, * statement2 ;
{ int type1, type2, nb_rows, i ;
  
  type1 = statement1->type ;
  type2 = statement2->type ;

   /* In the case of self-dependence, a statement commutes with hitself if
    * it is a reduction.
    */
   if ((statement1 == statement2) &&
       ((type1 == CANDL_P_REDUCTION) ||
        (type1 == CANDL_M_REDUCTION) ||
        (type1 == CANDL_T_REDUCTION)))
   return 1 ;
   
   /* Two statement commute when they are a reduction of the same type (or if
    * their left and right members are the same, but it's not exploited here).
    * The type may differ if it is either minus or plus-reduction. Furthermore,
    * they have to write onto the same array (and only one array).
    */
   if (((type1 == CANDL_P_REDUCTION) && (type2 == CANDL_P_REDUCTION)) ||
       ((type1 == CANDL_M_REDUCTION) && (type2 == CANDL_M_REDUCTION)) ||
       ((type1 == CANDL_T_REDUCTION) && (type2 == CANDL_T_REDUCTION)) ||
       ((type1 == CANDL_P_REDUCTION) && (type2 == CANDL_M_REDUCTION)) ||
       ((type1 == CANDL_M_REDUCTION) && (type2 == CANDL_P_REDUCTION)))
   { /* Here we check that there is one, only one and the same array. */
     nb_rows = statement1->written->NbRows ;
     if ((nb_rows == 0) || (nb_rows != statement2->written->NbRows))
     return 0 ;
     
     if (statement1->written->p[0][0] != statement2->written->p[0][0])
     return 0 ;
          
     for (i=1;i<nb_rows;i++)
     if ((statement1->written->p[i][0] != 0) ||
         (statement2->written->p[i][0] != 0))
     return 0 ;
   
     return 1 ;
   }
   
   return 0 ;
}






