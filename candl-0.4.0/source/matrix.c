
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                    matrix.c                             **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: december 9th 2005                **
    **--- |"-.-"| -------------------------------------------------------**
          |     |
          |     |
 ******** |     | *************************************************************
 * CAnDL  '-._,-' the Chunky Analyzer for Dependences in Loops (experimental) *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2005 Cedric Bastoul                                          *
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
 * candl_matrix_print_structure function:
 * Displays a CandlMatrix structure (matrix) into a file (file, possibly stdout)
 * in a way that trends to be understandable without falling in a deep
 * depression or, for the lucky ones, getting a headache... It includes an
 * indentation level (level) in order to work with others print_structure
 * functions.
 * - 09/12/2005: first version (from CLooG 0.14.0).
 */
void candl_matrix_print_structure(FILE * file, CandlMatrix * matrix, int level)
{ int i, j ;

  if (matrix != NULL)
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- CandlMatrix\n") ;

    for(j=0; j<=level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"%d %d\n",matrix->NbRows,matrix->NbColumns) ;

    /* Display the matrix. */
    for (i=0; i<matrix->NbRows; i++)
    { for(j=0; j<=level; j++)
      fprintf(file,"|\t") ;

      fprintf(file,"[ ") ;

      for (j=0; j<matrix->NbColumns; j++)
      { CANDL_print(file,CANDL_FMT,matrix->p[i][j]) ;
        fprintf(file," ") ;
      }

      fprintf(file,"]\n") ;
    }
  }
  else
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- NULL matrix\n") ;
  }

  /* The last line. */
  for(j=0; j<=level; j++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
}


/**
 * candl_matrix_print function:
 * This function prints the content of a CandlMatrix structure (matrix) into a
 * file (file, possibly stdout).
 * - 09/12/2005: first version (from CLooG 0.14.0).
 */
void candl_matrix_print(FILE * file, CandlMatrix * matrix)
{ candl_matrix_print_structure(file,matrix,0) ;
}



/**
 * candl_matrix_print_data function:
 * This function prints the content of a CandlMatrix data (matrix) into a
 * file (file, possibly stdout).
 */
void candl_matrix_print_data(FILE * file, CandlMatrix * matrix)
{
  int i, j;

  fprintf (file, "%d %d\n", matrix->NbRows, matrix->NbColumns);
  for (i = 0; i < matrix->NbRows; ++i)
    {
      for (j = 0; j < matrix->NbColumns; ++j)
	CANDL_print(file,CANDL_FMT,matrix->p[i][j]);
      fprintf (file, "\n");
    }
}


/**
 * candl_matrix_list_print_structure function:
 * Displays a CandlMatrixList structure (list) into a file (file, possibly
 * stdout) in a way that trends to be understandable without falling in a deep
 * depression or, for the lucky ones, getting a headache... It includes an
 * indentation level (level) in order to work with others print_structure
 * functions.
 * - 11/12/2005: first version.
 */
void candl_matrix_list_print_structure(file, list, level)
FILE * file ;
CandlMatrixList *list ;
int level ;
{ int i, j, first=1 ;

  if (list != NULL)
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- CandlMatrixList\n") ;
  }
  else
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- NULL matrix list\n") ;
  }

  while (list != NULL)
  { if (!first)
    { /* Go to the right level. */
      for(j=0; j<level; j++)
      fprintf(file,"|\t") ;
      fprintf(file,"|   CandlMatrixList\n") ;
    }
    else
    first = 0 ;

    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the matrix. */
    candl_matrix_print_structure(file,list->matrix,level+1) ;

    /* Next line. */
    if (list->next != NULL)
    { for(i=0; i<=level; i++)
      fprintf(file,"|\t") ;
      fprintf(file,"V\n") ;
    }
    list = list->next ;
  }

  /* The last line. */
  for(j=0; j<=level; j++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
}


/**
 * candl_matrix_list_print function:
 * This function prints the content of a CandlMatrixList structure (list) into a
 * file (file, possibly stdout).
 * - 11/12/2005: first version.
 */
void candl_matrix_list_print(FILE * file, CandlMatrixList * list)
{ candl_matrix_list_print_structure(file,list,0) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * candl_matrix_free function:
 * This function frees the allocated memory for a CandlMatrix structure.
 * - 09/12/2005: first version.
 */
void candl_matrix_free(CandlMatrix * matrix)
{ pip_matrix_free(matrix) ;
}


/**
 * candl_matrix_list_free function:
 * This function frees the allocated memory for a CandlMatrixList structure.
 * - 11/12/2005: first version.
 */
void candl_matrix_list_free(CandlMatrixList * list)
{ CandlMatrixList * next ;

  while (list != NULL)
  { next = list->next ;
    pip_matrix_free(list->matrix) ;
    free(list) ;
    list = next ;
  }
}


/******************************************************************************
 *                              Reading functions                             *
 ******************************************************************************/


/**
 * candl_matrix_read function:
 * This function reads a matrix into a file (foo, posibly stdin) and returns a
 * pointer to a CandlMatrix containing the read informations.
 * - 09/12/2005: first version.
 */
CandlMatrix * candl_matrix_read(FILE * file)
{ return pip_matrix_read(file) ;
}


/**
 * cloog_domain_list_read function:
 * This function reads a list of matrices into a file (foo, posibly stdin) and
 * returns a pointer to a CandlMatrixList containing the read information.
 * - 11/12/2005: first version (from CLooG 0.14.0's cloog_domain_list_read).
 */
CandlMatrixList * candl_matrix_list_read(FILE * file)
{ int i, nb_matrices ;
  char s[CANDL_MAX_STRING] ;
  CandlMatrixList * list, * now, * next ;

  /* We read first the number of matrices in the list. */
  while (fgets(s,CANDL_MAX_STRING,file) == 0) ;
  while ((*s=='#' || *s=='\n') || (sscanf(s," %d",&nb_matrices)<1))
  fgets(s,CANDL_MAX_STRING,file) ;

  /* Then we read the matrices. */
  list = NULL ;
  if (nb_matrices > 0)
  { list = (CandlMatrixList *)malloc(sizeof(CandlMatrixList)) ;
    list->matrix = candl_matrix_read(file) ;
    list->next = NULL ;
    now = list ;
    for (i=1;i<nb_matrices;i++)
    { next = (CandlMatrixList *)malloc(sizeof(CandlMatrixList)) ;
      next->matrix = candl_matrix_read(file) ;
      next->next = NULL ;
      now->next = next ;
      now = now->next ;
    }
  }

  return(list) ;
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * candl_matrix_malloc function:
 * This function allocates the memory space for a CandlMatrix structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - 09/12/2005: first version.
 */
CandlMatrix * candl_matrix_malloc(int nb_rows, int nb_columns)
{ return pip_matrix_alloc(nb_rows,nb_columns) ;
}


/**
 * candl_matrix_list_malloc function:
 * This function allocates the memory space for a CandlMatrixList structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - 11/12/2005: first version.
 */
CandlMatrixList * candl_matrix_list_malloc()
{ CandlMatrixList * list ;

  /* Memory allocation for the CandlDependence structure. */
  list = (CandlMatrixList *)malloc(sizeof(CandlMatrixList)) ;
  if (list == NULL)
  { fprintf(stderr, "[Candl]ERROR: memory overflow.\n") ;
    exit(1) ;
  }

  /* We set the various fields with default values. */
  list->matrix = NULL ;
  list->next   = NULL ;

  return list ;
}




/**
 * candl_matrix_violation function :
 * this function builds the constraint system corresponding to a violation of a
 * dependence, for a given transformation couple at a given depth.
 * - dependence is the constraint system of a dependence between two
     statements,
 * - t_source is the transformation function for the source statement,
 * - t_target is the transformation function for the target statement,
 * - dimension is the transformation dimension checked for legality,
 * - nb_par is the number of parameters.
 ***
 * - 13/12/2005: first version (extracted from candl_violation).
 */
CandlMatrix * candl_matrix_violation(dependence,t_source,t_target,
                                     dimension,nb_par)
CandlMatrix * dependence, * t_source, * t_target ;
int dimension, nb_par ;
{ int i, j, nb_rows, nb_columns, constraint, s_dims, t_dims ;
  CandlMatrix * system ;
  Entier temp ;

  CANDL_init(temp) ;

  /* The number of dimensions of the source and target domains. */
  s_dims = t_source->NbColumns - nb_par - 2 ;
  t_dims = t_target->NbColumns - nb_par - 2 ;

  /* Size of the constraint system. */
  nb_rows    = dependence->NbRows + dimension + 1 ;
  nb_columns = dependence->NbColumns ;

  /* We allocate memory space for the constraint system. */
  system = candl_matrix_malloc(nb_rows, nb_columns) ;

  /* We fill the constraint system (there is no need to put zeros in the
   * empty zones since candl_matrix_alloc initialized all to 0):
   */

  /* 1. We copy the constraints of the dependence polyhedron. */
  for (i = 0; i < dependence->NbRows; i++)
  for (j = 0; j < dependence->NbColumns; j++)
  CANDL_assign(system->p[i][j],dependence->p[i][j]) ;

  constraint = dependence->NbRows ;

  /* 2. We set the equality constraints (equality tag is already 0). */
  for (i = 0; i < dimension; i++)
  { /* The source dimension part. */
    for (j = 1; j <= s_dims; j++)
    CANDL_assign(system->p[constraint][j],t_source->p[i][j]) ;

    /* The -target dimension part. */
    for (; j <= s_dims + t_dims; j++)
    { CANDL_oppose(temp,t_target->p[i][j - s_dims]) ;
      CANDL_assign(system->p[constraint][j], temp) ;
    }

    /* The source-target parameter/scalar part. */
    for (; j < nb_columns; j++)
    CANDL_subtract(system->p[constraint][j],
                    t_source->p[i][j - t_dims],
                    t_target->p[i][j - s_dims]) ;
    constraint++ ;
  }

  /* 3. We set the target < source constraint. */
  /* This is an inequality. */
  CANDL_set_si(system->p[constraint][0], 1) ;

  /* The source dimension part. */
  for (j = 1; j<= s_dims; j++)
  CANDL_assign(system->p[constraint][j],t_source->p[dimension][j]) ;

  /* The -target dimension part. */
  for (; j<= s_dims + t_dims; j++)
  { CANDL_oppose(temp,t_target->p[dimension][j - s_dims]) ;
    CANDL_assign(system->p[constraint][j],temp) ;
  }

  /* The source-target parameter/scalar part. */
  for (; j < nb_columns; j++)
  CANDL_subtract(system->p[constraint][j],
                  t_source->p[dimension][j - t_dims],
                  t_target->p[dimension][j - s_dims]) ;
  /* We subtract 1 to the scalar to achieve >0 constraint. */
  CANDL_decrement(system->p[constraint][nb_columns - 1],
		  system->p[constraint][nb_columns - 1]) ;

  CANDL_clear(temp) ;
  return system ;
}




/**
 * candl_matrix_check_point function:
 * This function checks if there is an integral point in the set of
 * constraints, provided a given domain (possibly NULL).
 *
 */
int
candl_matrix_check_point (CandlMatrix* domain,
			  CandlMatrix* context)
{
  PipOptions* options;
  PipQuast* solution;
  int ret = 0;
  options = pip_options_init ();
  options->Simplify = 1;
  options->Urs_parms = -1;
  options->Urs_unknowns = -1;
  solution = pip_solve (domain, context, -1, options);

  if ((solution != NULL) &&
      ((solution->list != NULL) || (solution->condition != NULL)))
    ret = 1;
  pip_options_free (options);
  pip_quast_free (solution);

  return ret;
}
