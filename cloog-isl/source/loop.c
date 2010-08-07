
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                             loop.c                                **
    **-------------------------------------------------------------------**
    **                  First version: october 26th 2001                 **
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
 * want to know at each moment how many CloogLoop structures had been allocated
 * (cloog_loop_allocated) and how many had been freed (cloog_loop_freed).
 * Each time a CloogLoog structure is allocated, a call to the function
 * cloog_loop_leak_up() must be carried out, and respectively
 * cloog_loop_leak_down() when a CloogLoop structure is freed. The special
 * variable cloog_loop_max gives the maximal number of CloogLoop structures
 * simultaneously alive (i.e. allocated and non-freed) in memory.
 * - July 3rd->11th 2003: first version (memory leaks hunt and correction).
 */
 

static void cloog_loop_leak_up(CloogState *state)
{
  state->loop_allocated++;
  if ((state->loop_allocated - state->loop_freed) > state->loop_max)
    state->loop_max = state->loop_allocated - state->loop_freed;
}


static void cloog_loop_leak_down(CloogState *state)
{
  state->loop_freed++;
}


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_loop_print_structure function:
 * Displays a loop structure in a way that trends to be understandable without
 * falling in a deep depression or, for the lucky ones, getting a headache... 
 * Written by Olivier Chorier, Luc Marchaud, Pierre Martin and Romain Tartiere.
 * - April 24th 2005: Initial version.
 * - May   21rd 2005: - New parameter `F' for destination file (ie stdout),
 *                    - Minor tweaks.
 * - May   26th 2005: Memory leak hunt.
 * - June   2nd 2005: (Ced) Integration and minor fixes.
 * -June  22nd 2005: (Ced) Adaptation for GMP.
 */
void cloog_loop_print_structure(FILE * file, CloogLoop * loop, int level)
{ int i, j, first=1 ;

  if (loop)
  { /* Go to the right level. */
    for (i=0; i<level; i++)
    fprintf(file,"|\t") ;
    
    fprintf(file,"+-- CloogLoop\n") ;
  }
  
  /* For each loop. */
  while (loop)
  { if (!first)
    { /* Go to the right level. */
      for (i=0; i<level; i++)
      fprintf(file,"|\t") ;
    
      fprintf(file,"|   CloogLoop\n") ;
    }
    else
    first = 0 ;
    
    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Print the domain. */
    cloog_domain_print_structure(file, loop->domain, level+1, "CloogDomain");
    
    /* Print the stride. */
    for(j=0; j<=level; j++)
    fprintf(file,"|\t") ;
    fprintf(file, "Stride: ") ;
    cloog_int_print(file, loop->stride);
    fprintf(file, "\n") ;
        
    /* A blank line. */
    for(j=0; j<=level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;
    
    /* Print the block. */
    cloog_block_print_structure(file,loop->block,level+1) ;
    
    /* A blank line. */
    for (i=0; i<=level+1; i++)
    fprintf(file,"|\t") ;  
    fprintf(file,"\n") ;
    
    /* Print inner if any. */
    if (loop->inner)
    cloog_loop_print_structure(file,loop->inner,level+1) ;
    
    /* And let's go for the next one. */
    loop = loop->next ;

    /* One more time something that is here only for a better look. */
    if (!loop)
    { /* Two blank lines if this is the end of the linked list. */
      for (j=0; j<2; j++)
      { for (i=0; i<=level; i++)
        fprintf(file,"|\t") ;
      
        fprintf(file,"\n") ;
      }
    }
    else
    { /* A special blank line if the is a next loop. */
      for (i=0; i<=level; i++)
      fprintf(file,"|\t") ;  
      fprintf(file,"V\n") ;
    }
  }
}


/**
 * cloog_loop_print function:
 * This function prints the content of a CloogLoop structure (start) into a
 * file (file, possibly stdout).
 * - June 2nd 2005: Now this very old function (probably as old as CLooG) is
 *                  only a frontend to cloog_loop_print_structure, with a quite
 *                  better human-readable representation.
 */
void cloog_loop_print(FILE * file, CloogLoop * loop)
{ cloog_loop_print_structure(file,loop,0) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_loop_free function:
 * This function frees the allocated memory for a CloogLoop structure (loop),
 * and frees its inner loops and its next loops.
 * - June 22nd 2005: Adaptation for GMP.
 */
void cloog_loop_free(CloogLoop * loop)
{ CloogLoop * next ;
  
  while (loop != NULL) {
    cloog_loop_leak_down(loop->state);
    
    next = loop->next ;
    cloog_domain_free(loop->domain) ;
    cloog_block_free(loop->block) ;
    if (loop->inner != NULL)
    cloog_loop_free(loop->inner) ;
    
    cloog_int_clear(loop->stride);
    free(loop) ;
    loop = next ;
  }
}


/**
 * cloog_loop_free_parts function:
 * This function frees the allocated memory for some parts of a CloogLoop
 * structure (loop), each other argument is a boolean having to be set to 1 if
 * we want to free the corresponding part, 0 otherwise. This function applies
 * the same freeing policy to its inner ans next loops recursively.
 * - July  3rd 2003: first version.
 * - June 22nd 2005: Adaptation for GMP.
 */
void cloog_loop_free_parts(loop, domain, block, inner, next)
CloogLoop * loop ;
int domain, block, inner, next ;
{ CloogLoop * follow ;

  while (loop != NULL) {
    cloog_loop_leak_down(loop->state);
    follow = loop->next ;
    
    if (domain)
    cloog_domain_free(loop->domain) ;
       
    if (block)
    cloog_block_free(loop->block) ;
     
    if ((inner) && (loop->inner != NULL))
    cloog_loop_free_parts(loop->inner,domain,block,inner,1) ;
    
    cloog_int_clear(loop->stride);
    free(loop) ;
    if (next)
    loop = follow ;
    else
    loop = NULL ;
  }
}


/******************************************************************************
 *                              Reading functions                             *
 ******************************************************************************/


/**
 * cloog_loop_read function:
 * This function reads loop data into a file (foo, possibly stdin) and
 * returns a pointer to a CloogLoop structure containing the read information.
 * This function can be used only for input file reading, when one loop is
 * associated with one statement.
 * - number is the statement block number carried by the loop (-1 if none).
 * - nb_parameters is the number of parameters.
 **
 * - September 9th 2002: first version.
 * - April    16th 2005: adaptation to new CloogStatement struct (with number).
 * - June     11th 2005: adaptation to new CloogBlock structure. 
 * - June     22nd 2005: Adaptation for GMP.
 */
CloogLoop *cloog_loop_read(CloogState *state,
			    FILE * foo, int number, int nb_parameters)
{ int nb_iterators, op1, op2, op3 ;
  char s[MAX_STRING] ;
  CloogLoop * loop ;
  CloogStatement * statement ;

  cloog_loop_leak_up(state);
  
  /* Memory allocation and information reading for the first domain: */
  loop = (CloogLoop *)malloc(sizeof(CloogLoop)) ;
  if (loop == NULL) 
    cloog_die("memory overflow.\n");
  /* domain. */
  loop->state = state;
  loop->domain = cloog_domain_union_read(state, foo, nb_parameters);
  if (loop->domain != NULL)
    nb_iterators = cloog_domain_dimension(loop->domain);
  else
  nb_iterators = 0 ;
  /* stride is initialized to 1. */
  cloog_int_init(loop->stride);
  cloog_int_set_si(loop->stride, 1);
  /* included statement block. */
  statement = cloog_statement_alloc(state, number + 1);
  loop->block = cloog_block_alloc(statement, 0, NULL, nb_iterators);
  loop->usr = NULL;
  /* inner is NULL at beginning. */
  loop->inner = NULL ;
  /* next element. */
  loop->next = NULL ;  
    
  /* To read that stupid "0 0 0" line. */
  while (fgets(s,MAX_STRING,foo) == 0) ;
  while ((*s=='#' || *s=='\n') || (sscanf(s," %d %d %d",&op1,&op2,&op3)<3))
  fgets(s,MAX_STRING,foo) ;

  return loop ;
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * cloog_loop_malloc function:
 * This function allocates the memory space for a CloogLoop structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - November 21th 2005: first version.
 */
CloogLoop *cloog_loop_malloc(CloogState *state)
{ CloogLoop * loop ;
  
  /* Memory allocation for the CloogLoop structure. */
  loop = (CloogLoop *)malloc(sizeof(CloogLoop)) ;
  if (loop == NULL) 
    cloog_die("memory overflow.\n");
  cloog_loop_leak_up(state);
 
  
  /* We set the various fields with default values. */
  loop->state    = state;
  loop->domain = NULL ;
  loop->block  = NULL ;
  loop->usr    = NULL;
  loop->inner  = NULL ;
  loop->next   = NULL ;
  cloog_int_init(loop->stride);
  cloog_int_set_si(loop->stride, 1); 
  
  return loop ;
}  


/**
 * cloog_loop_alloc function:
 * This function allocates the memory space for a CloogLoop structure and
 * sets its fields with those given as input. Then it returns a pointer to the
 * allocated space.
 * - October  27th 2001: first version. 
 * - June     22nd 2005: Adaptation for GMP.
 * - November 21th 2005: use of cloog_loop_malloc.
 */ 
CloogLoop *cloog_loop_alloc(CloogState *state,
	CloogDomain *domain, cloog_int_t stride, CloogBlock *block,
	CloogLoop *inner, CloogLoop *next)
{ CloogLoop * loop ;
    
  loop = cloog_loop_malloc(state);
  
  loop->domain = domain ;
  loop->block  = block ;
  loop->inner  = inner ;
  loop->next   = next ;
  cloog_int_set(loop->stride, stride);
  
  return(loop) ;
}


/**
 * cloog_loop_add function:
 * This function adds a CloogLoop structure (loop) at a given place (now) of a
 * NULL terminated list of CloogLoop structures. The beginning of this list
 * is (start). This function updates (now) to (loop), and updates (start) if the
 * added element is the first one -that is when (start) is NULL-.
 * - October 28th 2001: first version. 
 */ 
void cloog_loop_add(CloogLoop ** start, CloogLoop ** now, CloogLoop * loop)
{ if (*start == NULL)
  { *start = loop ;
    *now = *start ;
  }
  else
  { (*now)->next = loop ;
    *now = (*now)->next ;
  }
}


/**
 * cloog_loop_add function:
 * This function adds a CloogLoop structure (loop) at a given place (now) of a
 * NULL terminated list of CloogLoop structures. The beginning of this list
 * is (start). This function updates (now) to the end of the loop list (loop),
 * and updates (start) if the added element is the first one -that is when
 * (start) is NULL-.
 * - September 9th 2005: first version.
 */ 
void cloog_loop_add_list(CloogLoop ** start, CloogLoop ** now, CloogLoop * loop)
{ if (*start == NULL)
  { *start = loop ;
    *now = *start ;
  }
  else
  { (*now)->next = loop ;
    *now = (*now)->next ;
  }
  
  while ((*now)->next != NULL)
  *now = (*now)->next ;
}


/**
 * cloog_loop_copy function:
 * This function returns a copy of the CloogLoop structure given as input. In
 * fact, there is just new allocations for the CloogLoop structures, but their
 * contents are the same.
 * - October 28th 2001: first version. 
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 */ 
CloogLoop * cloog_loop_copy(CloogLoop * source)
{ CloogLoop * loop ;
  CloogBlock * block ;
  CloogDomain * domain ;

  loop = NULL ;
  if (source != NULL)
  { domain = cloog_domain_copy(source->domain) ;
    block  = cloog_block_copy(source->block) ;
    loop   = cloog_loop_alloc(source->state,domain,source->stride,block,NULL,NULL);
    loop->usr = source->usr;
    loop->inner = cloog_loop_copy(source->inner) ;
    loop->next = cloog_loop_copy(source->next) ;
  }
  return(loop) ;
}


/**
 * cloog_loop_add_disjoint function:
 * This function adds some CloogLoop structures at a given place (now) of a
 * NULL terminated list of CloogLoop structures. The beginning of this list
 * is (start). (loop) can be an union of polyhedra, this function separates the
 * union into a list of *disjoint* polyhedra then adds the list. This function
 * updates (now) to the end of the list and updates (start) if first added
 * element is the first of the principal list -that is when (start) is NULL-.
 * (loop) can be freed by this function, basically when its domain is actually
 * a union of polyhedra, but don't worry, all the useful data are now stored
 * inside the list (start). We do not use PolyLib's Domain_Disjoint function,
 * since the number of union components is often higher (thus code size too).
 * - October   28th 2001: first version. 
 * - November  14th 2001: bug correction (this one was hard to find !).
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      22nd 2005: Adaptation for GMP.
 * - October   27th 2005: (debug) included blocks were not copied for new loops.
 */ 
void cloog_loop_add_disjoint(start, now, loop)
CloogLoop ** start, ** now, * loop ;
{
  CloogLoop * sep, * inner ;
  CloogDomain *domain, *seen, *seen_before, *temp, *rest;
  CloogBlock * block ;
  
  if (cloog_domain_isconvex(loop->domain))
  cloog_loop_add(start,now,loop) ;
  else {
    domain = cloog_domain_simplify_union(loop->domain);
    loop->domain = NULL ;
    
    /* We separate the first element of the rest of the union. */
    domain = cloog_domain_cut_first(domain, &rest);
            
    /* This first element is the first of the list of disjoint polyhedra. */
    sep = cloog_loop_alloc(loop->state, domain, loop->state->one,
			   loop->block, loop->inner, NULL);
    cloog_loop_add(start,now,sep) ;
  
    seen = cloog_domain_copy(domain);
    while (!cloog_domain_isempty(domain = rest)) {
      temp = cloog_domain_cut_first(domain, &rest);
      domain = cloog_domain_difference(temp, seen);
      cloog_domain_free(temp);

      if (cloog_domain_isempty(domain)) {
	cloog_domain_free(domain);
	continue;
      }
      
      /* Each new loop will have its own life, for instance we can free its
       * inner loop and included block. Then each one must have its own copy
       * of both 'inner' and 'block'.
       */
      inner = cloog_loop_copy(loop->inner) ;
      block = cloog_block_copy(loop->block) ;
      
      sep = cloog_loop_alloc(loop->state, cloog_domain_copy(domain),
			     loop->state->one, block, inner, NULL);
      /* domain can be an union too. If so: recursion. */
      if (cloog_domain_isconvex(domain))
	cloog_loop_add(start,now,sep) ;
      else
	cloog_loop_add_disjoint(start,now,sep) ;

      if (cloog_domain_isempty(rest)) {
	cloog_domain_free(domain);
	break;
      }
       
      seen_before = seen;
      seen = cloog_domain_union(seen_before, domain);
      cloog_domain_free(domain);
      cloog_domain_free(seen_before);
    }
    cloog_domain_free(rest);
    cloog_domain_free(seen);
    cloog_loop_free_parts(loop,0,0,0,0) ;  
  }
}


/**
 * cloog_loop_disjoint function:
 * This function returns a list of loops such that each loop with non-convex
 * domain in the input list (loop) is separated into several loops where the
 * domains are the components of the union of *disjoint* polyhedra equivalent
 * to the original non-convex domain. See cloog_loop_add_disjoint comments
 * for more details.
 * - September 16th 2005: first version.
 */
CloogLoop * cloog_loop_disjoint(CloogLoop * loop)
{ CloogLoop *res=NULL, * now=NULL, * next ;
  
  /* Because this is often the case, don't waste time ! */
  if ((loop != NULL) && cloog_domain_isconvex(loop->domain))
  return loop ;

  while (loop != NULL)
  { next = loop->next ;
    loop->next = NULL ;
    cloog_loop_add_disjoint(&res,&now,loop) ;
    loop = next ;
  }
  
  return res ;
}


/**
 * cloog_loop_restrict function:
 * This function returns the (loop) in the context of (context): it makes the
 * intersection between the (loop) domain and the (context), then it returns
 * a pointer to a new loop, with this intersection as domain.
 * - nb_par is the number of parameters.
 **
 * - October 27th 2001: first version. 
 * - June    15th 2005: a memory leak fixed (domain was not freed when empty).
 * - June    22nd 2005: Adaptation for GMP.
 */ 
CloogLoop * cloog_loop_restrict(loop, context, nb_par)
CloogLoop * loop ;
CloogDomain * context ;
int nb_par ;
{ int new_dimension ;
  CloogDomain * domain, * extended_context, * new_domain ;
  CloogLoop * new_loop ;
      
  domain = loop->domain ;
  if (cloog_domain_dimension(domain) > cloog_domain_dimension(context))
  {
    new_dimension = cloog_domain_dimension(domain);
    extended_context = cloog_domain_extend(context, new_dimension);
    new_domain = cloog_domain_intersection(extended_context,loop->domain) ;
    cloog_domain_free(extended_context) ;
  }
  else
  new_domain = cloog_domain_intersection(context,loop->domain) ;
  
  if (cloog_domain_isempty(new_domain))
  { cloog_domain_free(new_domain) ;
    return(NULL) ;
  }
  else {
    new_loop = cloog_loop_alloc(loop->state, new_domain,
				loop->state->one, loop->block, loop->inner, NULL);
    return(new_loop) ;
  }
}


/**
 * cloog_loop_project function:
 * This function returns the projection of (loop) on the (level) first
 * dimensions (outer loops). It makes the projection of the (loop) domain,
 * then it returns a pointer to a new loop, with this projection as domain.
 * - nb_par is the number of parameters.
 **
 * - October   27th 2001: first version. 
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      22nd 2005: Adaptation for GMP.
 */ 
CloogLoop * cloog_loop_project(CloogLoop * loop, int level, int nb_par)
{
  CloogDomain * new_domain ;
  CloogLoop * new_loop, * copy ;

  copy = cloog_loop_alloc(loop->state, loop->domain, loop->stride, loop->block,
                          loop->inner,NULL) ;

  if (cloog_domain_dimension(loop->domain) == level)
  new_domain = cloog_domain_copy(loop->domain) ;  
  else
    new_domain = cloog_domain_project(loop->domain, level);

  new_loop = cloog_loop_alloc(loop->state, new_domain, loop->state->one,
			      NULL, copy, NULL);
  
  return(new_loop) ;
}


/**
 * cloog_loop_concat function:
 * This function returns a pointer to the concatenation of the
 * CloogLoop lists given as input.
 * - October 28th 2001: first version. 
 */
CloogLoop * cloog_loop_concat(CloogLoop * a, CloogLoop * b)
{ CloogLoop * loop, * temp ;

  loop = a  ;
  temp = loop ;
  if (loop != NULL)
  { while (temp->next != NULL)
    temp = temp->next ;
    temp->next = b ;
  }
  else
  loop = b ;
  
  return(loop) ;
} 


/**
 * cloog_loop_combine:
 * Combine consecutive loops with identical domains into
 * a single loop with the concatenation of their inner loops
 * as inner loop.
 */
CloogLoop *cloog_loop_combine(CloogLoop *loop)
{
    CloogLoop *first, *second;

    for (first = loop; first; first = first->next) {
	while (first->next) {
	    if (!cloog_domain_lazy_equal(first->domain, first->next->domain))
		break;
	    second = first->next;
	    first->inner = cloog_loop_concat(first->inner, second->inner);
	    first->next = second->next;
	    cloog_loop_free_parts(second, 1, 0, 0, 0);
	}
    }

    return loop;
}

/**
 * cloog_loop_separate function:
 * This function implements the Quillere algorithm for separation of multiple
 * loops: for a given set of polyhedra (loop), it computes a set of disjoint
 * polyhedra such that the unions of these sets are equal, and returns this set.
 * - October   28th 2001: first version. 
 * - November  14th 2001: elimination of some unused blocks.
 * - August    13th 2002: (debug) in the case of union of polyhedra for one
 *                        loop, redundant constraints are fired.
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      22nd 2005: Adaptation for GMP.
 * - October   16th 2005: Removal of the non-shared constraint elimination when
 *                        there is only one loop in the list (seems to work
 *                        without now, DomainSimplify may have been improved).
 *                        The problem was visible with test/iftest2.cloog.
 */ 
CloogLoop * cloog_loop_separate(CloogLoop * loop)
{ int lazy_equal=0, disjoint = 0;
  CloogLoop * new_loop, * new_inner, * res, * now, * temp, * Q, 
            * inner, * old /*, * previous, * next*/  ;
  CloogDomain * UQ, * old_UQ, * domain ;
  
  if (loop == NULL)
  return NULL ;

  loop = cloog_loop_combine(loop);
  
  if (loop->next == NULL)
  return cloog_loop_disjoint(loop) ;
     
  UQ     = cloog_domain_copy(loop->domain) ;
  domain = cloog_domain_copy(loop->domain) ;
  res    = cloog_loop_alloc(loop->state, domain, loop->state->one,
			    loop->block, loop->inner, NULL);
  	  
  old = loop ;
  while((loop = loop->next) != NULL)
  { temp = NULL ;
    
    /* For all Q, add Q-loop associated with the blocks of Q alone,
     * and Q inter loop associated with the blocks of Q and loop.
     */
    for (Q = res; Q; Q = Q->next) {
        /* Add (Q inter loop). */
        if ((disjoint = cloog_domain_lazy_disjoint(Q->domain,loop->domain)))
	domain = NULL ;
	else
	{ if ((lazy_equal = cloog_domain_lazy_equal(Q->domain,loop->domain)))
	  domain = cloog_domain_copy(Q->domain) ;
          else
	  domain = cloog_domain_intersection(Q->domain,loop->domain) ;
          
	  if (!cloog_domain_isempty(domain))
          { new_inner = cloog_loop_concat(cloog_loop_copy(Q->inner),
                                          cloog_loop_copy(loop->inner)) ;
	    new_loop = cloog_loop_alloc(loop->state, domain, loop->state->one,
					NULL, new_inner, NULL);
            cloog_loop_add_disjoint(&temp,&now,new_loop) ;
          }
          else {
	    disjoint = 1;
	    cloog_domain_free(domain);
	  }
        }
        
	/* Add (Q - loop). */
        if (disjoint)
	domain = cloog_domain_copy(Q->domain) ;
	else
	{ if (lazy_equal)
	  domain = cloog_domain_empty(Q->domain);
	  else
	  domain = cloog_domain_difference(Q->domain,loop->domain) ;
	}
	
	if (!cloog_domain_isempty(domain)) {
          new_loop = cloog_loop_alloc(loop->state, domain, loop->state->one,
				      NULL, Q->inner, NULL);
          cloog_loop_add_disjoint(&temp,&now,new_loop) ;
        }
        else
        { cloog_domain_free(domain) ;
	  /* If Q->inner is no more useful, we can free it. */
          inner = Q->inner ;
          Q->inner = NULL ;
          cloog_loop_free(inner) ;
        }
    }

    /* Add loop-UQ associated with the blocks of loop alone.*/
    if (cloog_domain_lazy_disjoint(loop->domain,UQ))
    domain = cloog_domain_copy(loop->domain) ;
    else
    { if (cloog_domain_lazy_equal(loop->domain,UQ))
      domain = cloog_domain_empty(UQ);
      else
      domain = cloog_domain_difference(loop->domain,UQ) ;
    }
    
    if (!cloog_domain_isempty(domain)) {
      new_loop = cloog_loop_alloc(loop->state, domain, loop->state->one,
				  NULL, loop->inner, NULL);
      cloog_loop_add_disjoint(&temp,&now,new_loop) ;
    }
    else
    { cloog_domain_free(domain) ;
      /* If loop->inner is no more useful, we can free it. */
      cloog_loop_free(loop->inner) ;
    }
    
    loop->inner = NULL ;

    old_UQ = UQ ;
    if (loop->next != NULL)
    UQ = cloog_domain_union(UQ,loop->domain) ;

    cloog_domain_free(old_UQ) ;
    cloog_loop_free_parts(res,1,0,0,1) ;

    res = temp ;
  }  
  cloog_loop_free_parts(old,1,0,0,1) ;

  return(res) ;
}


/**
 * cloog_loop_merge_list
 * Merge two lists of CloogLoops.  The new list contains the
 * elements of the two lists in the same order, but they may
 * be interleaved.
 * In particular, if the elements of a and b are ordered
 * according to the inner loops of the order list, then so are the elements
 * in the new list.
 */
static CloogLoop *cloog_loop_merge_inner_list(CloogLoop *a, CloogLoop *b, 
					      CloogLoop *order)
{
  CloogLoop *loop, **next;
  next = &loop;

  for ( ; order && (a||b); order = order->next) {
    if (a && order->inner->block == a->block) {
      *next = a;
      a = a->next;
      next = &(*next)->next;
      continue;
    }
    if (b && order->inner->block == b->block) {
      *next = b;
      b = b->next;
      next = &(*next)->next;
    }
  }
  return loop;
}

/**
 * cloog_loop_merge function:
 * This function is the 'soft' version of loop_separate if we are looking for
 * a code much simpler (and less efficicient). Here we merge loops if they have
 * common parts in the iteration space (if the intersection of their domains is
 * not empty), and let them isolated otherwise. This function returns the new
 * CloogLoop list.
 * - October 29th 2001: first version. 
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      22nd 2005: Adaptation for GMP.
 */ 
CloogLoop * cloog_loop_merge(CloogLoop * loop, int nb_par, CloogOptions * options)
{
  CloogLoop * res, * merge, * now, * Q, * P, * new_inner, * next, * old ;
  CloogDomain * new_domain, * temp ;

  if (loop == NULL)
  return loop ;

  if (loop->next == NULL)
    return cloog_loop_disjoint(loop);

  /* First loop is added to the target list. */
  res = cloog_loop_alloc(loop->state, loop->domain, loop->state->one,
			 loop->block, loop->inner, NULL);
  old = loop ;
  /* Now the domain is in 'res' and it will be freed. */
  loop->domain = NULL ;
  
  /* And one by one, we see if we have to merge or to add the other loops. */
  while((loop = loop->next) != NULL)
  { merge = NULL ;
    P = cloog_loop_alloc(loop->state, loop->domain, loop->state->one,
			 loop->block, loop->inner, NULL);
    Q = res ;
    /* Now the domain is in 'P' and it will be freed. */
    loop->domain = NULL ;
        
    /* For each loop in the target list, if the intersection with the new loop
     * is empty, we can add the new loop directly, otherwise, we can merge then
     * add the fusion.
     */
    while (Q != NULL)
    { temp = cloog_domain_intersection(Q->domain,P->domain) ;
      next = Q->next ;
      if (cloog_domain_isempty(temp))
      { cloog_domain_free(temp) ;
        cloog_loop_add_disjoint(&merge,&now,Q) ;
      }
      else 
      { cloog_domain_free(temp) ;
        new_inner = cloog_loop_merge_inner_list(Q->inner, P->inner, old);
        temp = cloog_domain_union(P->domain,Q->domain) ;
	if (options->sh)
	  new_domain = cloog_domain_simple_convex(temp);
	else
	  new_domain = cloog_domain_convex(temp);
        cloog_domain_free(temp) ;
        /* Q and P are no more used (but their content yes !).*/
        cloog_loop_free_parts(P,1,0,0,0) ;
        cloog_loop_free_parts(Q,1,0,0,0) ;
        P = cloog_loop_alloc(loop->state, new_domain, loop->state->one,
			     NULL, new_inner, NULL);
      }
      Q = next ;
    }

    /* If there was merging, add it, otherwise add the loop lonely.
     * DEBUG : ici pas besoin de s'assurer que P->next est NULL (possible que
     * non si pas de fusion) car le dernier loop etudie a loop->next = NULL.
     */
    cloog_loop_add_disjoint(&merge,&now,P) ;
    res = merge ;
  }  
  cloog_loop_free_parts(old,0,0,0,1) ;

  return (res);
}


/**
 * cloog_loop_sort function:
 * Adaptation from LoopGen 0.4 by F. Quillere. This function sorts a list of
 * parameterized disjoint polyhedra, in order to not have lexicographic order
 * violation (see Quillere paper).
 * - September 16th 2005: inclusion of cloog_loop_number (October 29th 2001).
 */ 
CloogLoop * cloog_loop_sort(CloogLoop * loop, int level, int nb_par)
{ CloogLoop * res, * now, * temp, ** loop_array ;
  CloogDomain **doms;
  int i, nb_loops=0, * permut ;

  /* There is no need to sort the parameter domains. */
  if (!level)
    return loop;
  
  /* We will need to know how many loops are in the list. */
  temp = loop ;
  while (temp != NULL)
  { nb_loops ++ ;
    temp = temp->next ;
  }

  /* If there is only one loop, it's the end. */
  if (nb_loops == 1)
  return(loop) ;

  /* We have to allocate memory for some useful components:
   * - loop_array: the loop array,
   * - doms: the array of domains to sort,
   * - permut: will give us a possible sort (maybe not the only one).
   */
  loop_array = (CloogLoop **)malloc(nb_loops*sizeof(CloogLoop *)) ;
  doms = (CloogDomain **)malloc(nb_loops*sizeof(CloogDomain *));
  permut = (int *)malloc(nb_loops*sizeof(int)) ;

  /* We fill up the loop and domain arrays. */
  for (i=0;i<nb_loops;i++,loop=loop->next)
  { loop_array[i] = loop ;
    doms[i] = loop_array[i]->domain;
  }
  
  /* cloog_domain_sort will fill up permut. */
  cloog_domain_sort(doms, nb_loops, level, permut);
  
  /* With permut and loop_array we build the sorted list. */
  res = NULL ;
  for (i=0;i<nb_loops;i++)
  { /* To avoid pointer looping... loop_add will rebuild the list. */
    loop_array[permut[i]-1]->next = NULL ;
    cloog_loop_add(&res,&now,loop_array[permut[i]-1]) ;
  }
  
  free(permut) ;
  free(doms);
  free(loop_array) ;

  return res;
}


/**
 * cloog_loop_nest function:
 * This function changes the loop list in such a way that we have no more than
 * one dimension added by level. It returns an equivalent loop list with
 * this property.
 * - October 29th 2001: first version.
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      22nd 2005: Adaptation for GMP.
 * - November  21th 2005: (debug) now OK when cloog_loop_restrict returns NULL.
 */ 
CloogLoop * cloog_loop_nest(loop, context, level, nb_par)
CloogLoop * loop ;
CloogDomain * context ;
int level, nb_par ;
{ int l ;
  CloogLoop * p, * temp, * res, * now, * next ;
  CloogDomain * new_domain ;

  res = NULL ;
  /* Each domain is changed by its intersection with the context. */
  while (loop != NULL)
  { p = cloog_loop_restrict(loop,context,nb_par) ;
    next = loop->next ;
    
    if (p != NULL)
    { cloog_loop_free_parts(loop,1,0,0,0) ;
    
      temp = cloog_loop_alloc(p->state, p->domain, p->state->one,
			      p->block, p->inner, NULL);
    
      /* If the intersection dimension is too big, we make projections smaller
       * and smaller, and each projection includes the preceding projection
       * (thus, in the target list, dimensions are added one by one).
       */
      if (cloog_domain_dimension(p->domain) > level)
      for (l = cloog_domain_dimension(p->domain) - 1; l >= level; l--)
      { new_domain = cloog_domain_project(p->domain, l);
        temp = cloog_loop_alloc(p->state, new_domain, p->state->one,
				NULL, temp, NULL);
      }

       /* p is no more useful (but its content yes !). */ 
      cloog_loop_free_parts(p,0,0,0,0) ;
    
      cloog_loop_add(&res,&now,temp) ;
    }
    else
    cloog_loop_free_parts(loop,1,1,1,0) ;
    
    loop = next ;
  }
  
  return(res) ;
}


/**
 * cloog_loop_stride function:
 * This function will find the stride of a loop for the iterator at the column
 * number 'level' in the constraint matrix. It will update the lower bound of
 * the iterator accordingly. Basically, the function will try to find in the
 * inner loops a common condition on this iterator for the inner loop iterators
 * to be integral. For instance, let us consider a loop with the iterator i,
 * the iteration domain -4<=i<=n, and its two inner loops with the iterator j.
 * The first inner loop has the constraint 3j=i, and the second one has the
 * constraint 6j=i. Then the common constraint on i for j to be integral is
 * i%3=0, the stride for i is 3. Lastly, we have to find the new lower bound
 * for i: the first value satisfying the common constraint: -3. At the end, the
 * iteration domain for i is -3<=i<=n and the stride for i is 3.
 * - loop is the loop including the iteration domain of the considered iterator,
 * - level is the column number of the iterator in the matrix of contraints.
 **
 * - June 29th 2003: first version (work in progress since June 26th 2003).
 * - July 14th 2003: simpler version.
 * - June 22nd 2005: Adaptation for GMP (from S. Verdoolaege's 0.12.1 version).
 */
void cloog_loop_stride(CloogLoop * loop, int level, int nb_par)
{ int first_search ;
  cloog_int_t stride, ref_offset, offset, potential, lower;
  CloogLoop * inner ;

  cloog_int_init(stride);
  cloog_int_init(ref_offset);
  cloog_int_init(offset);
  cloog_int_init(potential);
  cloog_int_init(lower);

  cloog_int_set_si(ref_offset, 0);
  cloog_int_set_si(offset, 0);
  cloog_int_set_si(lower, 0);

  /* Default stride. */
  cloog_int_set_si(stride, 1);
  first_search = 1 ;
  inner = loop->inner ;
    
  if (cloog_domain_integral_lowerbound(loop->domain,level,&lower))
  while (inner != NULL)
  { /* If the minimun stride has not been found yet, find the stride. */
    if ((first_search) || (!cloog_int_is_one(stride)))
    {
      cloog_domain_stride(inner->domain, level, &potential, &offset);
      if (!cloog_int_is_one(potential) && (!first_search))
      { /* Offsets must be the same for common stride. */
	cloog_int_gcd(stride, potential, stride);
	if (!cloog_int_is_zero(stride)) {
	    cloog_int_fdiv_r(offset, offset, stride);
	    cloog_int_fdiv_r(ref_offset, ref_offset, stride);
	}
        if (cloog_int_ne(offset,ref_offset))
	    cloog_int_set_si(stride, 1);
      }
      else {
        cloog_int_set(stride, potential);
        cloog_int_set(ref_offset, offset);
      }
	
      first_search = 0 ;
    }
      
    inner = inner->next ;
  }

  if (cloog_int_is_zero(stride))
    cloog_int_set_si(stride, 1);
    
  /* Update the values if necessary. */
  if (!cloog_int_is_one(stride))
  { /* Update the stride value. */
    cloog_int_set(loop->stride, stride);
    /* The new lower bound l' is such that 
     *      (l' + offset) % s = 0 and l <= l' <= l+(s-1)
     * Let l' = k s - offset, then 
     *	    k s - offset <= l + (s-1) <= k s - offset + (s-1)
     * Or   l' = floor((l+offset+(s-1))/s) * s - offset
     *         = (floor((l+offset-1)/s) + 1) * s - offset
     */
    cloog_int_add(lower, lower, offset);
    cloog_int_sub_ui(lower, lower, 1);
    cloog_int_fdiv_q(lower, lower, stride);
    cloog_int_add_ui(lower, lower, 1);
    cloog_int_mul(lower, lower, stride);
    cloog_int_sub(lower, lower, offset);
    loop->domain = cloog_domain_lowerbound_update(loop->domain, level, lower);
  }
  
  cloog_int_clear(stride);
  cloog_int_clear(ref_offset);
  cloog_int_clear(offset);
  cloog_int_clear(potential);
  cloog_int_clear(lower);
}


/**
 * cloog_loop_stop function:
 * This function implements the 'stop' option : each domain of each loop
 * in the list 'loop' is replaced by 'context'. 'context' should be the
 * domain of the outer loop. By using this method, there are no more dimensions
 * to scan and the simplification step will automaticaly remove the domains
 * since they are the same as the corresponding contexts. The effect of this
 * function is to stop the code generation at the level this function is called,
 * the resulting code do not consider the next dimensions.
 * - January 11th 2005: first version.
 */
CloogLoop * cloog_loop_stop(CloogLoop * loop, CloogDomain * context)
{ if (loop == NULL)
  return NULL ;
  else
  { cloog_domain_free(loop->domain) ;
    loop->domain = cloog_domain_copy(context) ;
    loop->next = cloog_loop_stop(loop->next, context) ;
  }
  
  return loop ;
}


/**
 * cloog_loop_scalar_gt function:
 * This function returns 1 if loop 'l1' is greater than loop 'l2' for the
 * scalar dimension vector that begins at dimension 'scalar', 0 otherwise. What
 * we want to know is whether a loop is scheduled before another one or not.
 * This function solves the problem when the considered dimension for scheduling
 * is a scalar dimension. Since there may be a succession of scalar dimensions,
 * this function will reason about the vector of scalar dimension that begins
 * at dimension 'level+scalar' and finish to the first non-scalar dimension.
 * \param l1 Loop to be compared with l2.
 * \param l2 Loop to be compared with l1.
 * \param level Current non-scalar dimension.
 * \param scaldims Boolean array saying whether a dimension is scalar or not.
 * \param nb_scattdims Size of the scaldims array.
 * \param scalar Current scalar dimension.
 * \return 1 if (l1 > l2), 0 otherwise.
 **
 * - September 9th 2005: first version.
 * - October  15nd 2007: now "greater than" instead of "greater or equal".
 */
int cloog_loop_scalar_gt(l1, l2, level, scaldims, nb_scattdims, scalar)
CloogLoop * l1, * l2 ;
int level, * scaldims, nb_scattdims, scalar ;
{ while ((scalar < l1->inner->block->nb_scaldims) && scaldims[level+scalar-1])
  { if (cloog_int_gt(l1->inner->block->scaldims[scalar],
		     l2->inner->block->scaldims[scalar]))
    scalar ++ ;
    else
    return 0 ;
  }
  return 1 ;
}


/**
 * cloog_loop_scalar_eq function:
 * This function returns 1 if loop 'l1' is equal to loop 'l2' for the scalar
 * dimension vector that begins at dimension 'scalar', 0 otherwise. What we want
 * to know is whether two loops are scheduled for the same time or not.
 * This function solves the problem when the considered dimension for scheduling
 * is a scalar dimension. Since there may be a succession of scalar dimensions,
 * this function will reason about the vector of scalar dimension that begins
 * at dimension 'level+scalar' and finish to the first non-scalar dimension.
 * - l1 and l2 are the loops to compare,
 * - level is the current non-scalar dimension,
 * - scaldims is the boolean array saying whether a dimension is scalar or not,
 * - nb_scattdims is the size of the scaldims array,
 * - scalar is the current scalar dimension.
 **
 * - September 9th 2005 : first version.
 */
int cloog_loop_scalar_eq(l1, l2, level, scaldims, nb_scattdims, scalar)
CloogLoop * l1, * l2 ;
int level, * scaldims, nb_scattdims, scalar ;
{ while ((scalar < l1->inner->block->nb_scaldims) && scaldims[level+scalar-1])
  { if (cloog_int_eq(l1->inner->block->scaldims[scalar],
		     l2->inner->block->scaldims[scalar]))
    scalar ++ ;
    else
    return 0 ;
  }
  return 1 ;
}


/**
 * cloog_loop_scalar_sort function:
 * This function sorts a linked list of loops (loop) with respect to the
 * scalar dimension vector that begins at dimension 'scalar'. Since there may
 * be a succession of scalar dimensions, this function will reason about the
 * vector of scalar dimension that begins at dimension 'level+scalar' and
 * finish to the first non-scalar dimension.
 * \param loop Loop list to sort.
 * \param level Current non-scalar dimension.
 * \param scaldims Boolean array saying whether a dimension is scalar or not.
 * \param nb_scattdims Size of the scaldims array.
 * \param scalar Current scalar dimension.
 * \return A pointer to the sorted list.
 **
 * - July      2nd 2005: first developments.
 * - September 2nd 2005: first version.
 * - October  15nd 2007: complete rewrite to remove bugs, now a bubble sort.
 */
CloogLoop * cloog_loop_scalar_sort(loop, level, scaldims, nb_scattdims, scalar)
CloogLoop * loop ;
int level, * scaldims, nb_scattdims, scalar ;
{ int ok ;
  CloogLoop **current;
  
  do {
    ok = 1;
    for (current = &loop; (*current)->next; current = &(*current)->next) {
      CloogLoop *next = (*current)->next;
      if (cloog_loop_scalar_gt(*current,next,level,scaldims,nb_scattdims,scalar)) {
        ok = 0;
	(*current)->next = next->next;
	next->next = *current;
	*current = next;
      }
    }
  } while (!ok);

  return loop ;
}


/**
 * cloog_loop_generate_backtrack function:
 * adaptation from LoopGen 0.4 by F. Quillere. This function implements the
 * backtrack of the Quillere et al. algorithm (see the Quillere paper).
 * It eliminates unused iterations of the current level for the new one. See the
 * example called linearity-1-1 example with and without this part for an idea.
 * - October 26th 2001: first version in cloog_loop_generate_general.
 * - July    31th 2002: (debug) no more parasite loops (REALLY hard !). 
 * - October 30th 2005: extraction from cloog_loop_generate_general.
 */
CloogLoop * cloog_loop_generate_backtrack(loop, context, level, nb_par, options)
CloogLoop * loop ;
CloogDomain * context ;
int level, nb_par ;
CloogOptions * options ;
{
  CloogDomain * domain ;
  CloogLoop * now, * now2, * next, * next2, * end, * temp, * l, * inner,
            * new_loop ;
  
  temp = loop ;
  loop = NULL ;

  while (temp != NULL)
  { l = NULL ;
    inner = temp->inner ;
      
    while (inner != NULL)
    { next = inner->next ;
      /* This 'if' and its first part is the debug of july 31th 2002. */
      if (inner->block != NULL) {
        end = cloog_loop_alloc(temp->state, inner->domain, temp->state->one,
			       inner->block, NULL, NULL);
        domain = cloog_domain_copy(temp->domain) ;
        new_loop = cloog_loop_alloc(temp->state, domain, temp->state->one,
				    NULL, end, NULL);
      }
      else
      new_loop = cloog_loop_project(inner,level,nb_par) ;

      cloog_loop_free_parts(inner,0,0,0,0) ;
      cloog_loop_add(&l,&now2,new_loop) ;
      inner = next ;
    }
      
    temp->inner = NULL ;
      
    if (l != NULL)
    { l = cloog_loop_separate(l) ;
      l = cloog_loop_sort(l,level,nb_par) ;
      while (l != NULL) {
        cloog_int_set(l->stride, temp->stride);
        cloog_loop_add(&loop,&now,l) ;
        l = l->next ;
      }
    }
    next2 = temp->next ;
    cloog_loop_free_parts(temp,1,0,0,0) ;
    temp = next2 ;
  }

  return loop ;
}


/**
 * cloog_loop_generate_general function:
 * Adaptation from LoopGen 0.4 by F. Quillere. This function implements the
 * Quillere algorithm for polyhedron scanning from step 3 to 5.
 * (see the Quillere paper).
 * - loop is the loop for which we have to generate a scanning code,
 * - context is the context of the current loop (constraints on parameter and/or
 *   on outer loop counters),
 * - level is the current non-scalar dimension,
 * - scalar is the current scalar dimension,
 * - scaldims is the boolean array saying whether a dimension is scalar or not,
 * - nb_scattdims is the size of the scaldims array,
 * - nb_par is the number of parameters,
 * - options are the general code generation options.
 **
 * - October 26th 2001: first version.
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      22nd 2005: Adaptation for GMP.
 * - September  2nd 2005: The function have been cutted out in two pieces:
 *                        cloog_loop_generate and this one, in order to handle
 *                        the scalar dimension case more efficiently with
 *                        cloog_loop_generate_scalar.
 * - November  15th 2005: (debug) the result of the cloog_loop_generate call may
 *                        be a list of polyhedra (especially if stop option is
 *                        used): cloog_loop_add_list instead of cloog_loop_add.
 */ 
CloogLoop * cloog_loop_generate_general(loop, context, level, scalar,
                                        scaldims, nb_scattdims, nb_par, options)
CloogLoop * loop ;
CloogDomain * context ;
int level, scalar, * scaldims, nb_scattdims, nb_par ;
CloogOptions * options ;
{
  CloogLoop * res, * now, * temp, * l, * new_loop, * inner, * now2, * end,
            * next, * into ;
  CloogDomain * domain ;

  /* 3. Separate all projections into disjoint polyhedra. */
  res = (level && ((options->f > level+scalar) || (options->f < 0))) ?
        cloog_loop_merge(loop, nb_par, options) : cloog_loop_separate(loop);
    
  /* 3b. -correction- sort the loops to determine their textual order. */
  res = cloog_loop_sort(res,level,nb_par) ;
  
  /* 4. Recurse for each loop with the current domain as context. */
  temp = res ;
  res = NULL ;
  if (!level || (level+scalar < options->l) || (options->l < 0))
  while(temp != NULL)
  { if (level && options->strides)
    cloog_loop_stride(temp,level,nb_par) ;
    inner = temp->inner ;
    domain = temp->domain ;
    into = NULL ;
    while (inner != NULL)
    { /* 4b. -ced- recurse for each sub-list of non terminal loops. */
      if (cloog_domain_dimension(inner->domain) > level)
      { end = inner ;
        while ((end->next != NULL) &&
               (cloog_domain_dimension(end->next->domain) > level))
        end = end->next ;
        
	next = end->next ;
        end->next = NULL ;

        l = cloog_loop_generate(inner,domain,level+1,scalar,
	                        scaldims,nb_scattdims,nb_par,options) ;
        
	if (l != NULL)
        cloog_loop_add_list(&into,&now,l) ;
        
        inner = next ;
      }
      else
      { cloog_loop_add(&into,&now,inner) ;
        inner = inner->next ;
      }
    }
    next = temp->next ;
    temp->next = NULL ;
    temp->inner = into ;
    cloog_loop_add(&res,&now2,temp) ;
    temp = next ;
  }
  else
  while (temp != NULL)
  { next = temp->next ;
    l = cloog_loop_nest(temp->inner,temp->domain,level+1,nb_par) ;
    new_loop = cloog_loop_alloc(temp->state, temp->domain, temp->state->one,
				NULL, l, NULL);
    temp->inner = NULL ;
    temp->next = NULL ;
    cloog_loop_free_parts(temp,0,0,0,0) ;
    cloog_loop_add(&res,&now,new_loop) ;
    temp = next ;
  }
  
  /* 5. eliminate unused iterations of the current level for the new one. See
   *    the example called linearity-1-1 example with and without this part
   *    for an idea.
   */
  if ((!options->nobacktrack) && level &&
      ((level+scalar < options->l) || (options->l < 0)) &&
      ((options->f <= level+scalar) && !(options->f < 0)))
  res = cloog_loop_generate_backtrack(res,context,level,nb_par,options) ;
  
  /* Pray for my new paper to be accepted somewhere since the following stuff
   * is really amazing :-) !
   * Far long later: The paper has been accepted to PACT 2004 :-))). But there
   * are still some bugs and I have no time to fix them. Thus now you have to 
   * pray for me to get an academic position for that really amazing stuff :-) !
   * Later again: OK, I get my academic position, but still I have not enough
   * time to fix and clean this part... Pray again :-) !!!
   */
  /* res = cloog_loop_unisolate(res,context,level,nb_par) ;*/

  return(res) ;
}


/**
 * cloog_loop_generate_scalar function:
 * This function applies the simplified code generation scheme in the trivial
 * case of scalar dimensions. When dealing with scalar dimensions, there is 
 * no need of costly polyhedral operations for separation or sorting: sorting
 * is a question of comparing scalar vectors and separation amounts to consider
 * only loops with the same scalar vector for the next step of the code
 * generation process. This function achieves the separation/sorting process
 * for the vector of scalar dimension that begins at dimension 'level+scalar'
 * and finish to the first non-scalar dimension.
 * - loop is the loop for which we have to generate a scanning code,
 * - context is the context of the current loop (constraints on parameter and/or
 *   on outer loop counters),
 * - level is the current non-scalar dimension,
 * - scalar is the current scalar dimension,
 * - scaldims is the boolean array saying whether a dimension is scalar or not,
 * - nb_scattdims is the size of the scaldims array,
 * - nb_par is the number of parameters,
 * - options are the general code generation options.
 **
 * - September  2nd 2005: First version.
 */ 
CloogLoop * cloog_loop_generate_scalar(loop, context, level, scalar,
                                       scaldims, nb_scattdims, nb_par, options)
CloogLoop * loop ;
CloogDomain * context ;
int level, scalar, * scaldims, nb_scattdims, nb_par ;
CloogOptions * options ;
{ CloogLoop * res, * now, * temp, * l, * end, * next, * ref ;
  
  /* We sort the loop list with respect to the current scalar vector. */
  res = cloog_loop_scalar_sort(loop,level,scaldims,nb_scattdims,scalar) ;
  
  temp = res ;
  res = NULL ;
  while (temp != NULL)
  { /* Then we will appy the general code generation process to each sub-list
     * of loops with the same scalar vector.
     */
    end = temp ;
    ref = temp ;
    
    while((end->next != NULL) &&
         cloog_loop_scalar_eq(ref,end->next,level,scaldims,nb_scattdims,scalar))
    end = end->next ;

    next = end->next ;
    end->next = NULL ;

    /* For the next dimension, scalar value is updated by adding the scalar
     * vector size, which is stored at scaldims[level+scalar-1].
     */
    l = cloog_loop_generate_general(temp,context,level,
                                    scalar+scaldims[level+scalar-1],
	                            scaldims,nb_scattdims,nb_par,options) ;

    if (l != NULL)
    cloog_loop_add_list(&res,&now,l) ;
      
    temp = next ;
  }
  
  return res ;
}


/* Compare loop with the next loop based on their constant dimensions.
 * The result is < 0, == 0 or > 0 depending on whether the constant
 * dimensions of loop are lexicographically smaller, equal or greater
 * than those of loop->next.
 * If loop is the last in the list, then it is assumed to be smaller
 * than the "next" one.
 */
static int cloog_loop_next_scal_cmp(CloogLoop *loop)
{
    int i;
    int nb_scaldims;

    if (!loop->next)
	return -1;

    nb_scaldims = loop->block->nb_scaldims;
    if (loop->next->block->nb_scaldims < nb_scaldims)
	nb_scaldims = loop->next->block->nb_scaldims;

    for (i = 0; i < nb_scaldims; ++i) {
	int cmp = cloog_int_cmp(loop->block->scaldims[i],
				loop->next->block->scaldims[i]);
	if (cmp)
	    return cmp;
    }
    return loop->block->nb_scaldims - loop->next->block->nb_scaldims;
}


/* Check whether the globally constant dimensions of a and b
 * have the same value for all globally constant dimensions
 * that are situated before any (locally) non-constant dimension.
 */
static int cloog_loop_equal_prefix(CloogLoop *a, CloogLoop *b,
				    int *scaldims, int nb_scattdims)
{
    int i;
    int cst = 0;
    int dim = 0;

    for (i = 0; i < nb_scattdims; ++i) {
	if (!scaldims[i]) {
	    dim++;
	    continue;
	}
	if (!cloog_int_eq(a->block->scaldims[cst], b->block->scaldims[cst]))
	    break;
	cst++;
    }
    for (i = i + 1; i < nb_scattdims; ++i) {
	if (scaldims[i])
	    continue;
	if (!cloog_domain_lazy_isconstant(a->domain, dim))
	    return 0;
	/* No need to check that dim is also constant in b and that the
	 * constant values are equal.  That will happen during the check
	 * whether the two domains are equal.
	 */
	dim++;
    }
    return 1;
}


/* Try to block adjacent loops in the loop list "loop".
 * We only attempt blocking if the constant dimensions of the loops
 * in the least are (not necessarily strictly) increasing.
 * Then we look for a sublist such that the first (begin) has constant
 * dimensions strictly larger than the previous loop in the complete
 * list and such that the loop (end) after the last loop in the sublist
 * has constant dimensions strictly larger than the last loop in the sublist.
 * Furthermore, all loops in the sublist should have the same domain
 * (with globally constant dimensions removed) and the difference
 * (if any) in constant dimensions may only occur after all the
 * (locally) constant dimensions.
 * If we find such a sublist, then the blocks of all but the first
 * are merged into the block of the first.
 *
 * Note that this function can only be called before the global
 * blocklist has been created because it may otherwise modify and destroy
 * elements on that list.
 */
CloogLoop *cloog_loop_block(CloogLoop *loop, int *scaldims, int nb_scattdims)
{
    CloogLoop *begin, *end, *l;
    int begin_after_previous;
    int end_after_previous;

    if (!loop->next)
	return loop;
    for (begin = loop; begin; begin = begin->next) {
	if (!begin->block || !begin->block->scaldims)
	    return loop;
	if (cloog_loop_next_scal_cmp(loop) > 0)
	    return loop;
    }

    begin_after_previous = 1;
    for (begin = loop; begin; begin = begin->next) {
	if (!begin_after_previous) {
	    begin_after_previous = cloog_loop_next_scal_cmp(begin) < 0;
	    continue;
	}

	end_after_previous = cloog_loop_next_scal_cmp(begin) < 0;
	for (end = begin->next; end; end = end->next) {
	    if (!cloog_loop_equal_prefix(begin, end, scaldims, nb_scattdims))
		break;
	    if (!cloog_domain_lazy_equal(begin->domain, end->domain))
		break;
	    end_after_previous = cloog_loop_next_scal_cmp(end) < 0;
	}
	if (end != begin->next && end_after_previous) {
	    for (l = begin->next; l != end; l = begin->next) {
		cloog_block_merge(begin->block, l->block);
		begin->next = l->next;
		cloog_loop_free_parts(l, 1, 0, 1, 0);
	    }
	}

	begin_after_previous = cloog_loop_next_scal_cmp(begin) < 0;
    }

    return loop;
}


/**
 * cloog_loop_generate function:
 * Adaptation from LoopGen 0.4 by F. Quillere. This function implements the
 * Quillere algorithm for polyhedron scanning from step 1 to 2.
 * (see the Quillere paper).
 * - loop is the loop for which we have to generate a scanning code,
 * - context is the context of the current loop (constraints on parameter and/or
 *   on outer loop counters),
 * - level is the current non-scalar dimension,
 * - scalar is the current scalar dimension,
 * - scaldims is the boolean array saying whether a dimension is scalar or not,
 * - nb_scattdims is the size of the scaldims array,
 * - nb_par is the number of parameters,
 * - options are the general code generation options.
 **
 * - October 26th 2001: first version.
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - June      15th 2005: a memory leak fixed (loop was not entirely freed when
 *                        the result of cloog_loop_restrict was NULL).
 * - June      22nd 2005: Adaptation for GMP.
 * - September  2nd 2005: The function have been cutted out in two pieces:
 *                        cloog_loop_generate and this one, in order to handle
 *                        the scalar dimension case more efficiently with
 *                        cloog_loop_generate_scalar.
 * - November  15th 2005: (debug) Condition for stop option no more take care of
 *                        further scalar dimensions.
 */ 
CloogLoop * cloog_loop_generate(loop, context, level, scalar,
                                scaldims, nb_scattdims, nb_par, options)
CloogLoop * loop ;
CloogDomain * context ;
int level, scalar, * scaldims, nb_scattdims, nb_par ;
CloogOptions * options ;
{ CloogLoop * res, * now, * temp, * next, * old ;
  
  /* If the user asked to stop code generation at this level, let's stop. */
  if ((options->stop >= 0) && (level+scalar >= options->stop+1))
  return cloog_loop_stop(loop,context) ;

  res = NULL ;  
  
  /* 1. Replace each polyhedron by its intersection with the context.
   * 2. Compute the projection of each polyhedron onto the outermost
   *    loop variable and the parameters.
   */
  while (loop != NULL)
  { next = loop->next ;
    temp = cloog_loop_restrict(loop,context,nb_par) ;
    
    if (temp != NULL)
    { old = temp ;
      temp = cloog_loop_project(temp,level,nb_par) ;
      cloog_loop_free_parts(old,0,0,0,0) ;
      cloog_loop_add(&res,&now,temp) ;
      cloog_loop_free_parts(loop,1,0,0,0) ;
    }
    else
    { loop->next = NULL ;
      cloog_loop_free(loop) ;
    }
    
    loop = next ;
  }
  if (res == NULL)
  return NULL ;

  /* To save both time and memory, we switch here depending on whether the
   * current dimension is scalar (simplified processing) or not (general
   * processing).
   */
  if (level && (level+scalar <= nb_scattdims) && (scaldims[level+scalar-1]))
  res = cloog_loop_generate_scalar(res,context,level,scalar,
                                   scaldims,nb_scattdims,nb_par,options) ;
  else
  res = cloog_loop_generate_general(res,context,level,scalar,
                                    scaldims,nb_scattdims,nb_par,options) ;

  return res ;
}


/*
 * Internal function for simplifying a single loop in a list of loops.
 * See cloog_loop_simplify.
 */
static CloogLoop *loop_simplify(CloogLoop *loop, CloogDomain *context,
	int level, int nb_par)
{
  int domain_dim;
  CloogBlock * new_block ;
  CloogLoop *simplified, *inner;
  CloogDomain * domain, * simp, * inter, * extended_context ;

  domain = loop->domain ;
  
  domain_dim = cloog_domain_dimension(domain);
  extended_context = cloog_domain_extend(context, domain_dim);
  inter = cloog_domain_intersection(domain,extended_context) ;
  simp = cloog_domain_simplify(inter,extended_context) ;
  cloog_domain_free(extended_context) ;

  /* If the constraint system is never true, go to the next one. */
  if (cloog_domain_never_integral(simp)) {
    cloog_loop_free(loop->inner);
    cloog_domain_free(inter);
    cloog_domain_free(simp);
    return NULL;
  }

  inner = cloog_loop_simplify(loop->inner,inter,level+1,nb_par) ;
  cloog_domain_free(inter) ;
  
  if ((inner == NULL) && (loop->block == NULL)) {
    cloog_domain_free(simp);
    return NULL;
  }

  new_block = cloog_block_copy(loop->block) ;
  
  simplified = cloog_loop_alloc(loop->state,
				simp, loop->stride, new_block, inner, NULL);

  return(simplified) ; 
}


/**
 * cloog_loop_simplify function:
 * This function implements the part 6. of the Quillere algorithm, it
 * recursively simplifies each loop in the context of the preceding loop domain.
 * It returns a pointer to the simplified loop list.
 * The cloog_domain_simplify (DomainSimplify) behaviour is really bad with
 * polyhedra union and some really awful sidesteppings were written, I plan
 * to solve that...
 * - October   31th 2001: first version. 
 * - July 3rd->11th 2003: memory leaks hunt and correction.
 * - April     16th 2005: a memory leak fixed (extended_context was not freed).
 * - June      15th 2005: a memory leak fixed (loop was not conveniently freed
 *                        when the constraint system is never true).
 * - October   27th 2005: - this function called before cloog_loop_fast_simplify
 *                          is now the official cloog_loop_simplify function in
 *                          replacement of a slower and more complex one (after 
 *                          deep changes in the pretty printer).
 *                        - we use cloog_loop_disjoint to fix the problem when
 *                          simplifying gives a union of polyhedra (before, it
 *                          was under the responsibility of the pretty printer).
 */ 
CloogLoop * cloog_loop_simplify(loop, context, level, nb_par)
CloogLoop * loop ;
CloogDomain * context ;
int level, nb_par ;
{
  CloogLoop *now;
  CloogLoop *res = NULL;
  CloogLoop **next = &res;

  for (now = loop; now; now = now->next) {
    *next = loop_simplify(now, context, level, nb_par);

    now->inner = NULL; /* For loop integrity. */
    cloog_domain_free(now->domain);
    now->domain = NULL;

    if (*next)
      next = &(*next)->next;
  }
  cloog_loop_free(loop);
  
  /* Examples like test/iftest2.cloog give unions of polyhedra after
   * simplifying, thus we we have to disjoint them. Another good reason to
   * put the simplifying step in the Quillere backtrack.
   */
  res = cloog_loop_disjoint(res);

  return res;
}


/**
 * cloog_loop_scatter function:
 * This function add the scattering (scheduling) informations in a loop.
 */
void cloog_loop_scatter(CloogLoop * loop, CloogScattering *scatt)
{
  loop->domain = cloog_domain_scatter(loop->domain, scatt);
}

