
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                             names.c                               **
    **-------------------------------------------------------------------**
    **                  First version: august 1st 2002                   **
    **-------------------------------------------------------------------**/


/******************************************************************************
 *               CLooG : the Chunky Loop Generator (experimental)             *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2002-2005 Cedric Bastoul                                     *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License as published by the Free Software Foundation; either               *
 * version 2.1 of the License, or (at your option) any later version.         *
 *                                                                            *
 * This library is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
 * Lesser General Public License for more details.                            *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public           *
 * License along with this library; if not, write to the Free Software        *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,                         *
 * Boston, MA  02110-1301  USA                                                *
 *                                                                            *
 * CLooG, the Chunky Loop Generator                                           *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/
/* CAUTION: the english used for comments is probably the worst you ever read,
 *          please feel free to correct and improve it !
 */


# include <stdlib.h>
# include <stdio.h>
# include <ctype.h>
# include "../include/cloog/cloog.h"


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_names_print function:
 * this function is a human-friendly way to display the CloogNames data
 * structure, it shows all the different fields and includes an indentation
 * level (level) in order to work with others print_structure functions.
 * - July 1st 2005: first version based on the old cloog_names_print function,
 *                  it was the first modification in this file since two years !
 */
void cloog_names_print_structure(FILE * file, CloogNames * names, int level)
{ int i ;
  
  /* Go to the right level. */
  for (i=0; i<level; i++)
  fprintf(file,"|\t") ;
  
  if (names != NULL)
  { fprintf(file,"+-- CloogNames\n") ;
    
    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the scalar dimension number. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"Scalar dimension number ---: %d\n",names->nb_scalars) ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the scalar iterators. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    if (names->nb_scalars > 0)
    { fprintf(file,"+-- Scalar iterator strings:") ;
      for (i=0;i<names->nb_scalars;i++)
      fprintf(file," %s",names->scalars[i]) ;
      fprintf(file,"\n") ;
    }
    else
    fprintf(file,"+-- No scalar string\n") ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the scattering dimension number. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"Scattering dimension number: %d\n",names->nb_scattering) ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the scattering iterators. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    if (names->nb_scattering > 0)
    { fprintf(file,"+-- Scattering strings ----:") ;
      for (i=0;i<names->nb_scattering;i++)
      fprintf(file," %s",names->scattering[i]) ;
      fprintf(file,"\n") ;
    }
    else
    fprintf(file,"+-- No scattering string\n") ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;
    
    /* Print the iterator number. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"Iterator number -----------: %d\n",names->nb_iterators) ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the iterators. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    if (names->nb_iterators > 0)
    { fprintf(file,"+-- Iterator strings ------:") ;
      for (i=0;i<names->nb_iterators;i++)
      fprintf(file," %s",names->iterators[i]) ;
      fprintf(file,"\n") ;
    }
    else
    fprintf(file,"+-- No iterators\n") ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;
    
    /* Print the parameter number. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"Parameter number ----------: %d\n",names->nb_parameters) ;

    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the parameters. */
    for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
    if (names->nb_parameters > 0)
    { fprintf(file,"+-- Parameter strings -----:") ;
      for (i=0;i<names->nb_parameters;i++)
      fprintf(file," %s",names->parameters[i]) ;
      fprintf(file,"\n") ;
    }
    else
    fprintf(file,"No parameters\n") ;
    
  }
  else
  fprintf(file,"+-- No CloogNames\n") ;
  fprintf(file, "Number of active references: %d\n", names->references);
}


/**
 * cloog_names_print function:
 * This function prints the content of a CloogNames structure (names) into a
 * file (file, possibly stdout).
 * - July 1st 2005: Now this function is only a frontend to
 *                  cloog_program_print_structure, with a quite better
 *                  human-readable representation.
 */
void cloog_names_print(FILE * file, CloogNames * names)
{ cloog_names_print_structure(file,names,0) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_names_free function:
 * This function decrements the number of active references to 
 * a CloogNames structure and frees the allocated memory for this structure
 * if the count drops to zero.
 */
void cloog_names_free(CloogNames * names)
{ int i ;

  if (--names->references)
    return;

  if (names->scalars != NULL)
  { for (i=0;i<names->nb_scalars;i++)
    free(names->scalars[i]) ;
    free(names->scalars) ;
  }
   
  if (names->scattering != NULL)
  { for (i=0;i<names->nb_scattering;i++)
    free(names->scattering[i]) ;
    free(names->scattering) ;
  }
   
  if (names->iterators != NULL)
  { for (i=0;i<names->nb_iterators;i++)
    free(names->iterators[i]) ;
    free(names->iterators) ;
  }
   
  if (names->parameters != NULL)
  { for (i=0;i<names->nb_parameters;i++)
    free(names->parameters[i]) ;
    free(names->parameters) ;
  }
  free(names) ;
}


/**
 * cloog_names_copy function:
 * As usual in CLooG, "copy" means incrementing the reference count.
 */ 
CloogNames *cloog_names_copy(CloogNames *names)
{
  names->references++;
  return names;
}


/******************************************************************************
 *                              Reading functions                             *
 ******************************************************************************/


/**
 * cloog_names_read_strings function:
 * This function reads names data from a file (file, possibly stdin). It first
 * reads the naming option to know if whether has to automatically generate the
 * names, or to read them. Names are stored into an array of strings, and a
 * pointer to this array is returned.
 * - nb_items is the number of names the function will have to read if the
 *   naming option is set to read.
 * - prefix is the prefix to give to each name in case of automatic generation.
 * - first item is the name of the first suffix in case of automatic generation.
 **
 * - September 9th 2002: first version.
 */
char ** cloog_names_read_strings(file, nb_items, prefix, first_item)
FILE * file ;
int nb_items ;
char * prefix, first_item ;
{ int i, option, n ;
  char s[MAX_STRING], str[MAX_STRING], * c, ** names ;

  /* We first read name option. */
  while (fgets(s,MAX_STRING,file) == 0) ;
  while ((*s=='#' || *s=='\n') || (sscanf(s," %d",&option)<1))
  fgets(s,MAX_STRING,file) ;
  
  /* If there is no item to read, then return NULL. */
  if (nb_items == 0)
  return NULL ;
  
  /* If option is to read them in the file, then we do it and put them into
   * the array.
   */
  if (option)
  { /* Memory allocation. */
    names = (char **)malloc(nb_items*sizeof(char *)) ;
    if (names == NULL) 
      cloog_die("memory overflow.\n");
    for (i=0;i<nb_items;i++)
    { names[i] = (char *)malloc(MAX_NAME*sizeof(char)) ;
      if (names[i] == NULL) 
	cloog_die("memory overflow.\n");
    }
    
    do  /* Skip the comments, spaces and empty lines... */
    { c = fgets(s,MAX_STRING,file) ;
      while ((c != NULL) && isspace(*c) && (*c != '\n'))
      c++ ;
    }
    while (c != NULL && (*c == '#' || *c == '\n'));
    
    if (c == NULL) 
      cloog_die("no names in input file.\n");
    for (i=0;i<nb_items;i++) 
    { /* All names must be on the same line. */
      while (isspace(*c))
      c++ ;
      if (!*c || *c == '#' || *c == '\n')
        cloog_die("not enough names in input file.\n");
      /* n is strlen(str). */
      if (sscanf(c,"%s%n",str,&n) == 0) 
        cloog_die("no names in input file.\n");
      sscanf(str,"%s",names[i]) ;
      c += n ;
    }
  }
  /* Else we create names automatically. */
  else
  names = cloog_names_generate_items(nb_items,prefix,first_item) ;

  return names ;
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * cloog_names_malloc function:
 * This function allocates the memory space for a CloogNames structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - November 21th 2005: first version.
 */
CloogNames * cloog_names_malloc()
{ CloogNames * names ;
  
  /* Memory allocation for the CloogNames structure. */
  names = (CloogNames *)malloc(sizeof(CloogNames)) ;
  if (names == NULL) 
    cloog_die("memory overflow.\n");
  
  /* We set the various fields with default values. */
  names->nb_scalars    = 0 ;
  names->nb_scattering = 0 ;
  names->nb_iterators  = 0 ;
  names->nb_parameters = 0 ;
  names->scalars       = NULL ;
  names->scattering    = NULL ;
  names->iterators     = NULL ;
  names->parameters    = NULL ;
  names->references    = 1;
  
  return names ;
}  


/**
 * cloog_names_alloc function:
 * This function allocates the memory space for a CloogNames structure and
 * sets its fields with those given as input. Then it returns a pointer to the
 * allocated space.
 * - July       7th 2005: first version.
 * - September 11th 2005: addition of both scalar and scattering informations.
 * - November  21th 2005: use of cloog_names_malloc.
 */
CloogNames * cloog_names_alloc()
{ CloogNames * names ;

  /* Memory allocation for the CloogNames structure. */
  names = cloog_names_malloc() ;
  
  names->nb_scalars    = 0;
  names->nb_scattering = 0;
  names->nb_iterators  = 0;
  names->nb_parameters = 0;
  names->scalars       = NULL;
  names->scattering    = NULL;
  names->iterators     = NULL;
  names->parameters    = NULL;
  
  return names ;
}


/**
 * cloog_names_generate_items function:
 * This function returns a pointer to an array of strings with entries set
 * based on the function's parameters.
 * - nb_items will be the number of entries in the string array.
 * - prefix is the name prefix of each item or NULL.
 *   If not NULL, then the remainder of the name will be an integer
 *   in the range [0, nb_items-1].
 * - first_item is the name of the first item (if prefix == NULL),
 *   the nb_items-1 following items will be the nb_items-1
 *   following letters in ASCII code.
 **
 * - September 9th 2002 : first version, extracted from cloog_names_generate.
 */
char ** cloog_names_generate_items(int nb_items, char * prefix, char first_item)
{ int i ;
  char ** names ;
  
  if (nb_items == 0)
  return NULL ;
    
  names = (char **)malloc(nb_items*sizeof(char *)) ;
  if (names == NULL) 
    cloog_die("memory overflow.\n");
  for (i=0;i<nb_items;i++)
  { names[i] = (char *)malloc(MAX_NAME*sizeof(char)) ;
    if (names[i] == NULL) 
      cloog_die("memory overflow.\n");
    if (prefix == NULL)
    sprintf(names[i],"%c",first_item+i) ;
    else
      sprintf(names[i], "%s%d", prefix, 1+i);
  }
  
  return names ;
}


/**
 * cloog_names_generate function:
 * This function returns a pointer to a CloogNames structure with fields set
 * thanks to the function's parameters.
 * - nb_scalars will be the number of scalar dimensions in the structure.
 * - nb_scattering will be the number of scattering dimensions in the structure.
 * - nb_iterators will be the number of iterators in the CloogNames structure.
 * - nb_parameters will be the number of parameters in the CloogNames structure.
 * - first_s is the name of the first scalar iterator, the nb_scalars-1
 *   following iterators will be the nb_scalars-1 following letters in ASCII.
 * - first_t is the name of the first scattering iterator, the nb_scattering-1
 *   following iterators will be the nb_scattering-1 following letters in ASCII.
 * - first_i is the name of the first iterator, the nb_iterators-1 following
 *   iterators will be the nb_iterators-1 following letters in ASCII code.
 * - first_i is the name of the first iterator, the nb_iterators-1 following
 *   iterators will be the nb_iterators-1 following letters in ASCII code.
 * - first_p is the name of the first parameter, the nb_parameters-1 following
 *   parameters will be the nb_parameters-1 following letters in ASCII code.
 **
 * - July       1st 2002 : first version.
 * - September  9th 2002 : use of cloog_names_generate_items.
 * - September 11th 2005 : addition of both scalar and scattering informations.
 */
CloogNames * cloog_names_generate(
     nb_scalars, nb_scattering, nb_iterators, nb_parameters,
     first_s,    first_t,       first_i,      first_p)
int  nb_scalars, nb_scattering, nb_iterators, nb_parameters ;
char first_s,    first_t,       first_i,      first_p ;
{ CloogNames * names ;

  names = (CloogNames *)malloc(sizeof(CloogNames)) ;
  if (names == NULL) 
    cloog_die("memory overflow.\n");
  
  names->nb_scalars    = nb_scalars ;
  names->nb_scattering = nb_scattering ;
  names->nb_parameters = nb_parameters ;
  names->nb_iterators  = nb_iterators ;
  names->scalars       = cloog_names_generate_items(nb_scalars,   NULL,first_s);
  names->scattering    = cloog_names_generate_items(nb_scattering,NULL,first_t);
  names->parameters    = cloog_names_generate_items(nb_parameters,NULL,first_p);
  names->iterators     = cloog_names_generate_items(nb_iterators, NULL,first_i);

  return names ;
}


/* Lastly we update the CLoogNames structure: the iterators corresponding to
 * scalar dimensions have to be removed since these dimensions have been
 * erased and do not need to be print. We copy all the iterator names except
 * the scalar ones in a new string array.
 * - September 12th 2005: first version. 
 */
void cloog_names_scalarize(CloogNames * names, int nb_scattdims, int * scaldims)
{ int  nb_scalars, nb_scattering, i, current_scalar, current_scattering ;
  char ** scalars, ** scattering ;

  if (!nb_scattdims || (scaldims == NULL))
  return ;
  
  nb_scalars = 0 ;
  for (i=0;i<nb_scattdims;i++)
  if (scaldims[i])
  nb_scalars  ++ ;

  if (!nb_scalars)
  return ;
  
  nb_scattering = names->nb_scattering - nb_scalars ;
  scattering = (char **)malloc(nb_scattering * sizeof(char *)) ;
  if (scattering == NULL) 
    cloog_die("memory overflow.\n");
  scalars = (char **)malloc(nb_scalars * sizeof(char *)) ;
  if (scalars == NULL) 
    cloog_die("memory overflow.\n");
  
  current_scalar = 0 ;
  current_scattering  = 0 ;
  for (i=0;i<nb_scattdims;i++)
  { if (!scaldims[i])
    { scattering[current_scattering] = names->scattering[i] ;
      current_scattering ++ ;
    }
    else
    { scalars[current_scalar] = names->scattering[i] ;
      current_scalar ++ ;
    }
  }
  
  free(names->scattering) ;
  names->scattering    = scattering ;
  names->scalars       = scalars ;
  names->nb_scattering = nb_scattering ;
  names->nb_scalars    = nb_scalars ;
}

/**
 * Return the name at a given level (starting at one).
 * May be a scattering dimension or an iterator of the original domain.
 */
const char *cloog_names_name_at_level(CloogNames *names, int level)
{
  if (level <= names->nb_scattering)
    return names->scattering[level - 1];
  else
    return names->iterators[level - names->nb_scattering - 1];
}
