
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                             pprint.h                              **
    **-------------------------------------------------------------------**
    **                 First version: october 26th 2001                  **
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


#ifndef CLOOG_PPRINT_H
#define CLOOG_PPRINT_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 


# define MAX_STRING_VAL 32
# define INDENT_STEP 2

# define EQTYPE_NONE	 0
# define EQTYPE_CONSTANT 1
# define EQTYPE_PUREITEM 2
# define EQTYPE_EXAFFINE 3

# define LANGUAGE_C 0
# define LANGUAGE_FORTRAN 1


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void clast_pprint(FILE *foo, struct clast_stmt *root, int indent, 
	    CloogOptions *options);


#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
