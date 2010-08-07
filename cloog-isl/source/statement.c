
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                           statement.c                             **
    **-------------------------------------------------------------------**
    **                 First version: november 4th 2001                  **
    **-------------------------------------------------------------------**/


/******************************************************************************
 *               CLooG : the Chunky Loop Generator (experimental)             *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2001-2005 Cedric Bastoul                                     *
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
# include <string.h>
# include "../include/cloog/cloog.h"


/******************************************************************************
 *                             Memory leaks hunting                           *
 ******************************************************************************/


/**
 * These functions and global variables are devoted to memory leaks hunting: we
 * want to know at each moment how many CloogStatement structures had been
 * allocated (cloog_statement_allocated) and how many had been freed
 * (cloog_statement_freed). Each time a CloogStatement structure is allocated,
 * a call to the function cloog_statement_leak_up() must be carried out, and
 * respectively cloog_statement_leak_down() when a CloogStatement structure is
 * freed. The special variable cloog_statement_max gives the maximal number of
 * CloogStatement structures simultaneously alive (i.e. allocated and
 * non-freed) in memory.
 * - July 3rd->11th 2003: first version (memory leaks hunt and correction).
 */


static void cloog_statement_leak_up(CloogState *state)
{
  state->statement_allocated++;
  if ((state->statement_allocated - state->statement_freed) > state->statement_max)
  state->statement_max = state->statement_allocated - state->statement_freed ;
}


static void cloog_statement_leak_down(CloogState *state)
{ 
  state->statement_freed++;
}


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_domain_print_structure :
 * this function is a human-friendly way to display the CloogDomain data
 * structure, it includes an indentation level (level) in order to work with
 * others print_structure functions.
 * - June  16th 2005: first version.
 */
void cloog_statement_print_structure(file, statement, level)
FILE * file ;
CloogStatement * statement ;
int level ;
{ int i ;
      
  if (statement != NULL)
  { /* Go to the right level. */
    for (i=0; i<level; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"+-- CloogStatement %d \n",statement->number) ;
    
    statement = statement->next ;
 
    while (statement != NULL)
    { for (i=0; i<level; i++)
      fprintf(file,"|\t") ;
      fprintf(file,"|          |\n");
      for (i=0; i<level; i++)
      fprintf(file,"|\t") ;
      fprintf(file,"|          V\n");
      
      for (i=0; i<level; i++)
      fprintf(file,"|\t") ;
      fprintf(file,"|   CloogStatement %d \n",statement->number) ;
      statement = statement->next ;
    }
  }
  else
  { for (i=0; i<level; i++)
    fprintf(file,"|\t") ;
    
    fprintf(file,"+-- No CloogStatement\n") ;
  }  
}


/**
 * cloog_statement_print function:
 * This function prints the content of a CloogStatement structure (statement)
 * into a file (file, possibly stdout).
 */
void cloog_statement_print(FILE * file, CloogStatement * statement)
{ cloog_statement_print_structure(file,statement,0) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_statement_free function:
 * This function frees the allocated memory for a CloogStatement structure.
 */
void cloog_statement_free(CloogStatement * statement)
{ CloogStatement * next ;

  while (statement != NULL) {
    cloog_statement_leak_down(statement->state);
    
    next = statement->next ;
    /* free(statement->usr) ; Actually, this is user's job ! */
    free(statement) ;
    statement = next ;
  }
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * cloog_statement_malloc function:
 * This function allocates the memory space for a CloogStatement structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - November 21th 2005: first version.
 */
CloogStatement *cloog_statement_malloc(CloogState *state)
{ CloogStatement * statement ;
  
  /* Memory allocation for the CloogStatement structure. */
  statement = (CloogStatement *)malloc(sizeof(CloogStatement)) ;
  if (statement == NULL) 
    cloog_die("memory overflow.\n");
  cloog_statement_leak_up(state);
  
  /* We set the various fields with default values. */
  statement->state = state;
  statement->number = 0;
  statement->usr  = NULL ; /* To fill it is actually user's job ! */
  statement->next = NULL ;
  
  return statement ;
}  


/**
 * cloog_statement_alloc function:
 * This function allocates the memory space for a CloogStatement structure and
 * sets its fields with those given as input. Then it returns a pointer to the
 * allocated space.
 * - number is the statement number.
 **
 * - September 9th 2002: first version.
 * - March    17th 2003: fix for the usr field in CloogStatement structure.
 * - April    16th 2005: adaptation to new CloogStatement structure (with
 *                       number), cloog_statement_read becomes
 *                       cloog_statement_alloc sincethere is nothing more to
 *                       read on a file.
 * - November 21th 2005: use of cloog_statement_malloc.
 */
CloogStatement *cloog_statement_alloc(CloogState *state, int number)
{ CloogStatement * statement ;
    
  /* Memory allocation and initialization of the structure. */
  statement = cloog_statement_malloc(state);

  statement->number = number ;
  
  return statement ;
}


/**
 * cloog_statement_copy function:
 * This function returns a copy of the CloogStatement structure given as input.
 * - October 28th 2001: first version (in loop.c). 
 * - March   17th 2003: fix for the usr field in CloogStatement structure.
 * - April   16th 2005: adaptation to new CloogStatement struct (with number). 
 */ 
CloogStatement * cloog_statement_copy(CloogStatement * source)
{ CloogStatement * statement, * temp, * now = NULL ;
  
  statement = NULL ;

  while (source != NULL) {
    cloog_statement_leak_up(source->state);

    temp = (CloogStatement *)malloc(sizeof(CloogStatement)) ;
    if (temp == NULL)
      cloog_die("memory overflow.\n");
    
    temp->state  = source->state;
    temp->number = source->number ;
    temp->usr    = source->usr ;
    temp->next   = NULL ;
    
    if (statement == NULL)
    { statement = temp ;
      now = statement ;
    }
    else
    { now->next = temp ;
      now = now->next ;
    }
    source = source->next ;
  }
  return(statement) ;
}


/** 
 * cloog_statement_add function:
 * This function adds a CloogStatement structure (statement) at a given place
 * (now) of a NULL terminated list of CloogStatement structures. The beginning
 * of this list is (start). This function updates (now) to (loop), and
 * updates (start) if the added element is the first one -that is when (start)
 * is NULL-.
 * - March 27th 2004: first version. 
 */ 
void cloog_statement_add(start, now, statement)
CloogStatement ** start, ** now, * statement ;
{ if (*start == NULL)
  { *start = statement ;
    *now = *start ;
  }
  else
  { (*now)->next = statement ;
    *now = (*now)->next ;
  }
}

