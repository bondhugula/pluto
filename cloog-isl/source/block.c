
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                             block.c                               **
    **-------------------------------------------------------------------**
    **                   First version: june 11th 2005                   **
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
# include "../include/cloog/cloog.h"


/******************************************************************************
 *                             Memory leaks hunting                           *
 ******************************************************************************/


/**
 * These functions and global variables are devoted to memory leaks hunting: we
 * want to know at each moment how many CloogBlock structures had been allocated
 * (cloog_block_allocated) and how many had been freed (cloog_block_freed).
 * Each time a CloogBlock structure is allocated, a call to the function
 * cloog_block_leak_up() must be carried out, and respectively
 * cloog_block_leak_down() when a CloogBlock structure is freed. The special
 * variable cloog_block_max gives the maximal number of CloogBlock structures
 * simultaneously alive (i.e. allocated and non-freed) in memory.
 * - June 11th 2005: first version.
 */


static void cloog_block_leak_up(CloogState *state)
{ 
  state->block_allocated++;
  if ((state->block_allocated - state->block_freed) > state->block_max)
    state->block_max = state->block_allocated - state->block_freed;
}


static void cloog_block_leak_down(CloogState *state)
{
  state->block_freed++;
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
void cloog_block_print_structure(FILE * file, CloogBlock * block, int level)
{ int i ;
  
  /* Go to the right level. */
  for (i=0; i<level; i++)
  fprintf(file,"|\t") ;
    
  if (block != NULL)
  { fprintf(file,"+-- CloogBlock\n") ;
   
    /* A blank line. */
    for (i=0; i<level+2; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;    

    /* Print statement list. */
    cloog_statement_print_structure(file,block->statement,level+1) ;
    
    /* A blank line. */
    for (i=0; i<level+2; i++)
    fprintf(file,"|\t") ;    
    fprintf(file,"\n") ;    

    /* Print scattering function. */
    for(i=0; i<level+1; i++)
    fprintf(file,"|\t") ;
  
    fprintf(file,"+-- Null scattering function\n") ;

    /* A blank line. */
    for (i=0; i<level+2; i++)
    fprintf(file,"|\t") ;    
    fprintf(file,"\n") ;    

    /* Print scalar dimensions. */
    for (i=0; i<level+1; i++)
    fprintf(file,"|\t") ;
    
    if (block->nb_scaldims == 0)
    fprintf(file,"No scalar dimensions\n") ;
    else
    { fprintf(file,"Scalar dimensions (%d):",block->nb_scaldims) ;
      for (i = 0; i < block->nb_scaldims; i++) {
	fprintf(file, " ");
	cloog_int_print(file, block->scaldims[i]);
      }
      fprintf(file,"\n") ;
    }
    
    /* A blank line. */
    for (i=0; i<level+2; i++)
    fprintf(file,"|\t") ;    
    fprintf(file,"\n") ;    

    /* Print depth. */
    for (i=0; i<level+1; i++)
    fprintf(file,"|\t") ;
    fprintf(file,"Depth: %d\n",block->depth) ;
    
    /* A blank line. */
    for (i=0; i<level+1; i++)
    fprintf(file,"|\t") ;   
    fprintf(file,"\n") ;    
  }
  else
  fprintf(file,"+-- Null CloogBlock\n") ;
}


/**
 * cloog_block_print function:
 * This function prints the content of a CloogBlock structure (block) into a
 * file (file, possibly stdout).
 * - June 11th 2005: first version.
 * - June  16th 2005: now just a call to cloog_block_print_structure.
 */
void cloog_block_print(FILE * file, CloogBlock * block)
{ cloog_block_print_structure(file,block,0) ;  
}


/**
 * cloog_block_list_print function:
 * This function prints the content of a CloogBlock structure (block) into a
 * file (file, possibly stdout).
 * - June 16th 2005: first version.
 */
void cloog_block_list_print(FILE * file, CloogBlockList * blocklist)
{ int i=0 ;
  
  while (blocklist != NULL)
  { fprintf(file,"+-- CloogBlockList node %d\n",i) ;
    cloog_block_print_structure(file,blocklist->block,1) ;
    blocklist = blocklist->next ;
    i++ ;
  }
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_block_free function:
 * This function frees the allocated memory for a CloogStatement structure.
 * - June 11th 2005: first version.
 * - June 30th 2005: scaldims field management.
 */
void cloog_block_free(CloogBlock * block)
{ int i ;

  if (block != NULL)
  { block->references -- ;
    
    if (block->references == 0) {
      cloog_block_leak_down(block->state);
      if (block->scaldims != NULL)
      { for (i=0;i<block->nb_scaldims;i++)
	  cloog_int_clear(block->scaldims[i]);
      
        free(block->scaldims) ;
      }
      if (block->statement)
	cloog_statement_free(block->statement);
      free(block) ;
    }
  }
}


/**
 * cloog_block_list_free function:
 * This function frees the allocated memory for a CloogBlockList structure.
 * - June 11th 2005: first version.
 */
void cloog_block_list_free(CloogBlockList * blocklist)
{ CloogBlockList * temp ;
  
  while (blocklist != NULL)
  { temp = blocklist->next ;
    cloog_block_free(blocklist->block) ;
    free(blocklist) ;
    blocklist = temp ;
  }
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/

/**
 * cloog_block_malloc function:
 * This function allocates the memory space for a CloogBlock structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - November 21th 2005: first version.
 */
CloogBlock *cloog_block_malloc(CloogState *state)
{ CloogBlock * block ;
  
  /* Memory allocation for the CloogBlock structure. */
  block = (CloogBlock *)malloc(sizeof(CloogBlock)) ;
  if (block == NULL) 
    cloog_die("memory overflow.\n");
  cloog_block_leak_up(state);
  
  /* We set the various fields with default values. */
  block->state = state;
  block->statement = NULL ;
  block->nb_scaldims = 0 ;
  block->scaldims = NULL ;
  block->depth = 0 ;
  block->references = 1 ;
  block->usr = NULL;
  
  return block ;
}  


/**
 * cloog_block_alloc function:
 * This function allocates the memory space for a CloogBlock structure and
 * sets its fields with those given as input. Then it returns a pointer to the
 * allocated space. The two parameters nb_scaldims and scaldims are for internal
 * service, put to respectively 0 and NULL if you don't know what they are
 * useful for !
 * - statement is the statement list of the block,
 * - scattering is the scattering function for the block (NULL if unsure !),
 * - nb_scaldims is the number of scalar dimensions (0 if unsure !),
 * - scaldims is the array with the scalar dimensions values (NULL if unsure !),
 * - depth is the original block depth (the number of outer loops).
 **
 * - June     11th 2005: first version.
 * - June     30th 2005: addition of the nb_scaldims and scaldims parameters.
 * - November 21th 2005: use of cloog_block_malloc.
 */
CloogBlock *cloog_block_alloc(CloogStatement *statement, int nb_scaldims,
				cloog_int_t *scaldims, int depth)
{ CloogBlock * block ;
    
  /* Block allocation. */
  block = cloog_block_malloc(statement->state);

  block->statement = statement ;
  block->nb_scaldims = nb_scaldims ;
  block->scaldims = scaldims ;
  block->depth = depth ;
  block->references = 1 ;
  
  return block ;
}


/**
 * cloog_block_list_malloc function:
 * This function allocates the memory space for a CloogBlockList structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - November 21th 2005: first version.
 */
CloogBlockList * cloog_block_list_malloc()
{ CloogBlockList * blocklist ;
  
  /* Memory allocation for the CloogBlock structure. */
  blocklist = (CloogBlockList *)malloc(sizeof(CloogBlockList)) ;
  if (blocklist == NULL) 
    cloog_die("memory overflow.\n");
  
  /* We set the various fields with default values. */
  blocklist->block = NULL ;
  blocklist->next  = NULL ;
  
  return blocklist ;
}  


/**
 * cloog_block_list_alloc function:
 * This function allocates the memory space for a CloogBlockList structure and
 * sets its fields with those given as input. Then it returns a pointer to the
 * allocated space.
 * - block is the block element of the list node,
 **
 * - June     11th 2005: first version.
 * - November 21th 2005: use of cloog_block_list_malloc.
 */
CloogBlockList * cloog_block_list_alloc(CloogBlock * block)
{ CloogBlockList * blocklist ;
  
  /* Block list node allocation. */
  blocklist = cloog_block_list_malloc() ;

  blocklist->block = block ;
  blocklist->block->references ++ ; /* The block has a new reference to it. */
  blocklist->next = NULL ;
  
  return blocklist ;
}


/**
 * cloog_block_copy function:
 * This function returns a copy of a CloogBlock structure 'block'. To save
 * memory this is not a memory copy but we increment a counter of active
 * references inside the structure, then return a pointer to that structure.
 */ 
CloogBlock * cloog_block_copy(CloogBlock * block)
{ if (block == NULL)
  return NULL ;

  block->references ++ ;
  return block ;
}


/**
 * cloog_block_merge function:
 * this function adds at the end of the statement list of the block 'block',
 * the statement list of the block 'merged'. Then the  CloogBlock structure
 * of 'merged' is freed (obviously not its statement list that is now
 * included in 'block').
 * - June 11th 2005: first version.
 */
void cloog_block_merge(CloogBlock * block, CloogBlock * merged)
{ CloogStatement * statement ;

  if ((block == NULL) || (merged == NULL))
  return ;
  
  if (block->statement != NULL)
  { statement = block->statement ;
    
    while (statement->next != NULL)
    statement = statement->next ;
    
    statement->next = merged->statement ;
  }
  else
  block->statement = merged->statement ;

  merged->statement = NULL;
  cloog_block_free(merged);
}










