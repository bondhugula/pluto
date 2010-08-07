
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                  violation.c                            **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: december 12th 2005               **
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
# include <string.h>
# include "../include/candl/candl.h"


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * candl_violation_print_structure function:
 * Displays a CandlViolation structure (violation) into a file (file,
 * possibly stdout) in a way that trends to be understandable without falling
 * in a deep depression or, for the lucky ones, getting a headache... It
 * includes an indentation level (level) in order to work with others
 * print_structure functions.
 * - 18/09/2003: first version.
 */
void candl_violation_print_structure(file, violation, level)
FILE * file ;
CandlViolation * violation ;
int level ;
{ int j, first=1 ;
  CandlDependence * next=NULL ;

  if (violation != NULL)
  { /* Go to the right level. */
    for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- CandlViolation\n") ;
  }
  else
  { for(j=0; j<level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- NULL dependence violation\n") ;
  }

  while (violation != NULL)
  { if (!first)
    { /* Go to the right level. */
      for(j=0; j<level; j++)
      fprintf(file,"|\t") ;
      fprintf(file,"|   CandlViolation\n") ;
    }
    else
    first = 0 ;

    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Go to the right level and print the dimension. */
    for(j=0; j<=level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"Dimension: %d\n",violation->dimension) ;

    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the dependence. */
    if (violation->dependence != NULL)
    { next = violation->dependence->next ; /* To not print the whole list... */
      violation->dependence->next = NULL ; /* I know it's not beautiful :-/ ! */
    }
    candl_dependence_print_structure(file,violation->dependence,level+1) ;
    if (violation->dependence != NULL)
    violation->dependence->next = next ;

    /* Print the dependence polyhedron. */
    candl_matrix_print_structure(file,violation->domain,level+1) ;

    violation = violation->next ;

    /* Next line. */
    if (violation != NULL)
    { for (j=0; j<=level; j++)
      fprintf(file,"|\t") ;
      fprintf(file,"V\n") ;
    }
  }

  /* The last line. */
  for(j=0; j<=level; j++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
}


/* candl_violation_print function:
 * This function prints the content of a CandlViolation structure
 * (violation) into a file (file, possibly stdout).
 */
void candl_violation_print(FILE * file, CandlViolation * violation)
{ candl_violation_print_structure(file,violation,0) ;
}


/* candl_violation_pprint function:
 * This function prints the content of a CandlViolation structure (violation)
 * into a file (file, possibly stdout) as a Graphviz input file.
 * See http://www.graphviz.org
 * - 12/12/2005: first version.
 */
void candl_violation_pprint(FILE * file, CandlViolation * violation)
{ int i=0 ;
  CandlDependence * dependence ;

  fprintf(file,"digraph G {\n") ;

  fprintf(file,"# Legality Violation Graph\n") ;
  fprintf(file,"# Generated by Candl "CANDL_RELEASE" "CANDL_VERSION" bits\n");
  if (violation == NULL)
  fprintf(file,"# Congratulations: the transformation is legal !\n");

  while (violation != NULL)
  { dependence = violation->dependence ;

    fprintf(file,"  S%d -> S%d [label=\" ",dependence->source->label,
                                           dependence->target->label) ;
    switch (dependence->type)
    { case CANDL_UNSET : fprintf(file,"UNSET") ; break ;
      case CANDL_RAW   : fprintf(file,"RAW")   ; break ;
      case CANDL_WAR   : fprintf(file,"WAR")   ; break ;
      case CANDL_WAW   : fprintf(file,"WAW")   ; break ;
      case CANDL_RAR   : fprintf(file,"RAR")   ; break ;
      default : fprintf(file,"unknown") ;
    }
    fprintf(file," depth %d, ref %d->%d, viol %d \"];\n",
            dependence->depth,
            dependence->ref_source,
            dependence->ref_target,
            violation->dimension) ;
    violation = violation->next ;
    i++ ;
  }

  if (i>4)
  fprintf(file,"# Number of edges = %i\n}\n",i) ;
  else
  fprintf(file,"}\n") ;
}


/* candl_violation_view function:
 * This function uses dot (see http://www.graphviz.org) and gv (see
 * http://wwwthep.physik.uni-mainz.de/~plass/gv) tools to display the
 * violation graph.
 * - 20/03/2006: first version.
 */
void candl_violation_view(CandlViolation * violation)
{ FILE * temp_output ;

  temp_output = fopen(CANDL_TEMP_OUTPUT,"w") ;
  candl_violation_pprint(temp_output,violation) ;
  fclose(temp_output) ;
  system("(dot -Tps "CANDL_TEMP_OUTPUT" | gv - &) && rm -f "CANDL_TEMP_OUTPUT) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/

/* candl_violation_free function:
 * This function frees the allocated memory for a CandlViolation structure.
 * - 18/09/2003: first version.
 */
void candl_violation_free(CandlViolation * violation)
{ CandlViolation * next ;

  while (violation != NULL)
  { next = violation->next ;
    candl_matrix_free(violation->domain) ;
    free(violation) ;
    violation = next ;
  }
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * candl_violation_malloc function:
 * This function allocates the memory space for a CandlViolation structure
 * and sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - 07/12/2005: first version.
 */
CandlViolation * candl_violation_malloc()
{ CandlViolation * violation ;

  /* Memory allocation for the CandlViolation structure. */
  violation = (CandlViolation *)malloc(sizeof(CandlViolation)) ;
  if (violation == NULL)
  { fprintf(stderr, "[Candl]ERROR: memory overflow.\n") ;
    exit(1) ;
  }

  /* We set the various fields with default values. */
  violation->dependence = NULL ;
  violation->domain     = NULL ;
  violation->next       = NULL ;
  violation->dimension  = CANDL_UNSET ;

  return violation ;
}


/* cloog_violation_add function:
 * This function adds a CloogViolation structure (violation) at a given place
 * (now) of a NULL terminated list of CloogViolation structures. The beginning
 * of this list is (start). This function updates (now) to the end of the loop
 * list (loop), and updates (start) if the added element is the first one -that
 * is when (start) is NULL-.
 * - 12/12/2005: first version (from candl_dependence_add).
 */
void candl_violation_add(start, now, violation)
CandlViolation ** start, ** now, * violation ;
{ if (violation != NULL)
  { if (*start == NULL)
    { *start = violation ;
      *now = *start ;
    }
    else
    { (*now)->next = violation ;
      *now = (*now)->next ;
    }

    while ((*now)->next != NULL)
    *now = (*now)->next ;
  }
}


/**
 * candl_violation function :
 * this function will build the list of violated dependences by a program
 * transformation candidate, according to some user options. It returns
 * the linked list of violations.
 * - program containt the program and transformation candidate informations,
 * - dependence is the dependence graph, if NULL it will be calculated here,
 * - options is the user options data structure.
 **
 * - 12/12/2005: first version.
 */
CandlViolation * candl_violation(program, dependence, options)
CandlProgram * program ;
CandlDependence * dependence ;
CandlOptions * options ;
{ int dimension, max_dimension, violated ;
  CandlMatrix * system, * domain, * t_source, * t_target ;
  CandlStatement * source, * target ;
  CandlViolation * violation=NULL, * now=NULL, * new ;
  PipOptions * pip_options ;
  PipQuast * solution ;

  /* If there is no program or transformation, we consider this legal. */
  if ((program == NULL) || (program->transformation == NULL))
  return NULL ;

  /* If the dependence graph is not already built, do it. */
  if (dependence == NULL)
  dependence = candl_dependence(program,options) ;

  pip_options = pip_options_init() ;
  pip_options->Simplify = 1 ;

  /* We check every edge of the dependence graph. */
  while (dependence != NULL)
  { source = dependence->source ;
    target = dependence->target ;
    domain = dependence->domain ;

    /* We find the source transformation matrix. */
    t_source = program->transformation[source->label] ;

    /* We find the target transformation matrix. */
    t_target = program->transformation[target->label] ;

    /* The maximal dimension we have to check for legality. */
    max_dimension = CANDL_min(t_source->NbRows,t_target->NbRows) ;

    /* We check each dimension for legality. */
    for (dimension = 0; dimension<max_dimension; dimension++)
    { violated = 0 ;
      system = NULL ;

      /* We build the constraint system corresponding to that
       * violation then check if there is an integral point inside,
       * if yes there is actually a dependence violation and we
       * will add this one to the list.
       */
      system = candl_matrix_violation(dependence->domain,t_source,t_target,
                                      dimension,program->context->NbColumns-2) ;

      solution = pip_solve(system,program->context,-1,pip_options) ;

      if ((solution != NULL) &&
          ((solution->list != NULL) || (solution->condition != NULL)))
      violated = 1 ;

      pip_quast_free(solution) ;

      if (violated)
      { new = candl_violation_malloc() ;

	/* We set the various fields with corresponding values. */
        new->dependence = dependence ;
	new->dimension = dimension ;
	new->domain = system ;

	candl_violation_add(&violation,&now,new) ;

	if (!options->fullcheck)
	{ pip_options_free(pip_options) ;
          return violation ;
	}
      }
      else
      candl_matrix_free(system) ;
    }
    dependence = dependence->next ;
  }

  pip_options_free(pip_options) ;
  return violation ;
}

