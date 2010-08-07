
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                  dependence.h                           **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: september 18th 2003              **
    **--- |"-.-"| -------------------------------------------------------**
          |     |
          |     |
 ******** |     | *************************************************************
 * CAnDL  '-._,-' the Chunky Analyzer for Dependences in Loops (experimental) *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2003-2008 Cedric Bastoul                                     *
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


#ifndef CANDL_DEPENDENCE_H
# define CANDL_DEPENDENCE_H

# include <stdio.h>
# include <candl/statement.h>
# include <candl/matrix.h>
# include <candl/program.h>
# include <candl/options.h>

# define CANDL_ARRAY_BUFF_SIZE		2048
# define CANDL_VAR_UNDEF		1
# define CANDL_VAR_IS_DEF		2
# define CANDL_VAR_IS_USED		3
# define CANDL_VAR_IS_DEF_USED		4


# if defined(__cplusplus)
extern "C"
  {
# endif


/**
 * CandlDependence structure:
 * this structure contains all the informations about a data dependence, it is
 * also a node of the linked list of all dependences of the dependence graph.
 */
struct candldependence
{ CandlStatement  * source;      /**< Pointer to source statement. */
  CandlStatement  * target;      /**< Pointer to target statement. */
  int depth;                     /**< Dependence level. */
  int type;                      /**< Dependence type: a dependence from source
                                   *   to target can be:
				   *   - CANDL_UNSET if the dependence type is
				   *     still not set,
				   *   - CANDL_RAW if source writes M and
				   *     target read M (flow-dependence),
				   *   - CANDL_WAR if source reads M and
				   *     target writes M (anti-dependence),
				   *   - CANDL_WAW if source writes M and
				   *     target writes M too (output-dependence)
				   *   - CANDL_RAR if source reads M and
				   *     target reads M too (input-dependence).
				   */
  int ref_source;                /**< Position of source reference. */
  int ref_target;                /**< Position of target reference. */
  CandlMatrix * domain;          /**< Dependence polyhedron. */

  void* usr;			 /**< User field, for library users
				    convenience. */
  struct candldependence * next; /**< Pointer to next dependence */
};
typedef struct candldependence CandlDependence;
typedef struct candldependence candl_dependence_t;
typedef struct candldependence * candl_dependence_p;


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
void candl_dependence_print_structure(FILE *, candl_dependence_p, int);
void candl_dependence_print(FILE *, candl_dependence_p);
void candl_dependence_pprint(FILE *, candl_dependence_p);
void candl_dependence_view(candl_dependence_p);
# ifdef CANDL_SUPPORTS_SCOPLIB
CandlDependence* candl_dependence_read_from_scop(scoplib_scop_p, CandlProgram*);
void candl_dependence_update_scop_with_deps(scoplib_scop_p, CandlDependence*);
void candl_dependence_print_scop(FILE*, FILE*, CandlDependence*);
# endif


/******************************************************************************
 *                         Memory alloc/dealloc function                      *
 ******************************************************************************/
candl_dependence_p      candl_dependence_malloc();
void			candl_dependence_free(candl_dependence_p);


/******************************************************************************
 *                             Processing functions                           *
 ******************************************************************************/
int			candl_dependence_gcd_test(CandlStatement*,
						  CandlStatement*,
						  CandlMatrix*, int);
int			candl_dependence_check(CandlProgram *,
					       candl_dependence_p,
					       CandlOptions *);
candl_dependence_p      candl_dependence(CandlProgram *, CandlOptions *);


/******************************************************************************
 *                          Scalar analysis functions                         *
 ******************************************************************************/
int
candl_dependence_var_is_scalar (candl_program_p, int);

CandlStatement**
candl_dependence_refvar_chain(candl_program_p, CandlStatement*, int, int);

int
candl_dependence_var_is_ref(CandlStatement*, int);

int
candl_dependence_check_domain_is_included(CandlStatement*, CandlStatement*,
					  CandlMatrix*, int);

int
candl_dependence_scalar_is_privatizable_at(candl_program_p, int, int);

int
candl_dependence_is_loop_carried (candl_program_p, CandlDependence*, int);

void
candl_dependence_prune_scalar_waw (candl_program_p, CandlOptions*,
				   CandlDependence**);

void
candl_dependence_prune_with_privatization (candl_program_p, CandlOptions*,
					   CandlDependence**);

int
candl_dependence_scalar_renaming(candl_program_p, CandlOptions*,
				 CandlDependence**);

int
candl_dependence_analyze_scalars(candl_program_p, CandlOptions*);

int
candl_num_dependences(CandlDependence *candl_deps);

void
candl_compute_last_writer (CandlDependence *dep, CandlProgram *prog);

# if defined(__cplusplus)
  }
# endif
#endif /* define CANDL_DEPENDENCE_H */

