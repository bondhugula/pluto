
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                           statement.h                             **
    **-------------------------------------------------------------------**
    **                  First version: november 4th 2001                 **
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


#ifndef CLOOG_STATEMENT_H
#define CLOOG_STATEMENT_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 


struct cloogstatement
{
  CloogState *state;             /* State. */
  int number;                    /* The statement unique number. */
  void * usr ;                   /* A pointer for library users convenience. */
  struct cloogstatement * next ; /* Pointer to the next statement with the
                                  * same original domain and the same
				  * scattering function.
				  */
} ;
typedef struct cloogstatement CloogStatement ;


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void cloog_statement_print_structure(FILE *, CloogStatement *, int) ;
void cloog_statement_print(FILE *, CloogStatement *) ;


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/
void cloog_statement_free(CloogStatement *) ;


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
CloogStatement * cloog_statement_malloc(CloogState *state);
CloogStatement * cloog_statement_alloc(CloogState *state, int);
CloogStatement * cloog_statement_copy(CloogStatement *) ;
void cloog_statement_add(CloogStatement**, CloogStatement**, CloogStatement*) ;

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */

