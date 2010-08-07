/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                  type.h                                    *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
 *                                                                            *
 * This library is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at    *
 * your option) any later version.                                            *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this library; if not, write to the Free Software Foundation,    *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA         *
 *                                                                            *
 * Written by Paul Feautrier                                                  *
 *                                                                            *
 ******************************************************************************/

#ifndef TYPE_H
#define TYPE_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

/* Modified by Serge Torres to handle very big problems (since 1.3.4 we can put
 * any value we want: sol_space is allocated dynamically), but it is left by
 * default to 4096 because of time/space reasons for most people.
 * #define SOL_SIZE 67108864
 */
#define SOL_SIZE 4096

extern Entier UN;
extern Entier ZERO;

#define Pip_True 1
#define Pip_False 0


#ifdef TC
#define DEBUG 8
#endif

#define Q if(cross_product>=limit)

#define MAXCOL 1024
#define MAXPARM 50

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
