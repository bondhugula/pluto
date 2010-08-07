/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 funcall.h                                  *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
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
 *****************************************************************************/

#ifndef FUNCALL_H
#define FUNCALL_H
#if defined(__cplusplus)
extern "C" 
  {
#endif 

#define TRAITER_INT		(1 << 0)	/* Compute integer optimum */
#define TRAITER_DUAL		(1 << 1)	/* Compute dual variables */
void traiter(Tableau *tp, Tableau *ctxt, int nvar, int nparm, int ni, int nc,
	     int bigparm, int flags);
int integrer(Tableau **, Tableau **, int *, int *, int *, int *, int);
#if defined(LINEAR_VALUE_IS_MP)
#else
Entier pgcd(Entier, Entier);
Entier mod(Entier,Entier);
int llog(Entier);
#endif

int dgetc(FILE *foo);
FILE *pip_create_dump_file();
int sol_hwm(void);
void sol_simplify(int);
int is_not_Nil(int);
int sol_edit(FILE *, int);
void tab_reset(struct high_water_mark);
void sol_reset(int);
struct high_water_mark tab_hwm(void);
Tableau *tab_get(FILE *, int,int,int);
int tab_simplify(Tableau *tp, int cst);
void sol_init(void);
void sol_close(void);
void tab_init(void);
void tab_close(void);
void sol_if(void);
void sol_forme(int);
void sol_val(Entier, Entier);
void sol_nil(void);
void sol_error(int);
Tableau * tab_alloc(int, int, int);
void sol_list(int);
void tab_display(Tableau *, FILE *);
Tableau * expanser(Tableau *, int, int, int, int, int, int);
void sol_new(int);
void sol_div(void);

#if defined(__cplusplus)
  }
#endif 
#endif /* define _H */
