
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 matrix.c                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 30/04/2008                     **
    **- [""M# | #  U"U#U  -----------------------------------------------**
         | #  | #  \ .:/
         | #  | #___| #
 ******  | "--'     .-"  ******************************************************
 *     |"-"-"-"-"-#-#-##   Clan : the Chunky Loop Analyzer (experimental)     *
 ****  |     # ## ######  *****************************************************
 *      \       .::::'/                                                       *
 *       \      ::::'/     Copyright (C) 2008 Cedric Bastoul                  *
 *     :8a|    # # ##                                                         *
 *     ::88a      ###      This is free software; you can redistribute it     *
 *    ::::888a  8a ##::.   and/or modify it under the terms of the GNU Lesser *
 *  ::::::::888a88a[]:::   General Public License as published by the Free    *
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 3 of the       *
 *::::::::8::::888:Y8888:: License, or (at your option) any later version.    *
 *::::':::88::::888::Y88a::::::::::::...                                      *
 *::'::..    .   .....   ..   ...  .                                          *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.							      *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with software; if not, write to the Free Software Foundation, Inc.,  *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Clan, the Chunky Loop Analyzer                                             *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <ctype.h>
# include <clan/matrix.h>


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * clan_matrix_tag_array function:
 * this function tags a matrix to explicit it is describing the array index of
 * a given array. This means using SCoP representation that the very first
 * element of the very first row will be the array number instead of zero.
 * It updates directly the matrix provided as parameter.
 * \param matrix The matrix to tag.
 * \param array  The array number.
 **
 * - 02/05/2008: first version.
 */
void
clan_matrix_tag_array(scoplib_matrix_p matrix, int array)
{
  if ((matrix == NULL) || (matrix->NbRows == 0))
  {
    fprintf(stderr,"[Clan] Error: matrix cannot be array-tagged\n");
    exit(1);
  }

  SCOPVAL_set_si(matrix->p[0][0],array);
}


/**
 * clan_matrix_scheduling function:
 * this function builds the scheduling matrix for the clan_statement_t
 * structures thanks to the parser current state of parser_scheduling (rank)
 * and parser_depth (depth). The "rank" vector gives the "position" of the
 * statement for every loop depth (see Feautrier's demonstration of existence
 * of a schedule for any SCoP or CLooG's manual for original scattering
 * function to understand if necessary). This function just "expands" this
 * vector to a (2*n+1)-dimensional schedule for a statement at depth n and
 * returns it.
 * \param rank  The position of the statement at every loop depth.
 * \param depth The depth of the statement.
 **
 * - 01/05/2008: First version.
 */
scoplib_matrix_p
clan_matrix_scheduling(int * rank, int depth)
{
  int i, j, nb_rows, nb_columns;
  scoplib_matrix_p scheduling;

  nb_rows    = 2 * depth + 1;
  nb_columns = CLAN_MAX_DEPTH + CLAN_MAX_PARAMETERS + 2;
  scheduling = scoplib_matrix_malloc(nb_rows,nb_columns);

  j = 0;
  for (i = 0; i < depth; i++)
  {
    SCOPVAL_set_si(scheduling->p[j][nb_columns-1],rank[i]);
    SCOPVAL_set_si(scheduling->p[j+1][i+1],1);
    j += 2;
  }
  SCOPVAL_set_si(scheduling->p[nb_rows-1][nb_columns-1],rank[depth]);

  return scheduling;
}


/**
 * clan_matrix_compact function:
 * This function compacts a matrix such that it uses the right number
 * of columns (during construction we used CLAN_MAX_DEPTH and
 * CLAN_MAX_PARAMETERS to define matrix and vector sizes). It modifies
 * directly the matrix provided as parameter.
 * \param matrix        The matrix to compact.
 * \param nb_iterators  The true number of iterators for this matrix.
 * \param nb_parameters The true number of parameters in the SCoP.
 **
 * - 02/05/2008: first version.
 * - 24/05/2008: nice bug fixed (p_Init_size was not copied, segfaulting later).
 */
void
clan_matrix_compact(scoplib_matrix_p matrix, int nb_iterators, 
		    int nb_parameters)
{
  int i, j, nb_columns;
  scoplib_matrix_p compacted;

  if (matrix == NULL)
    return;

  nb_columns = nb_iterators + nb_parameters + 2;
  compacted = scoplib_matrix_malloc(matrix->NbRows,nb_columns);

  for (i = 0; i < matrix->NbRows; i++)
  {
    /* We copy the equality/inequality tag and the iterator coefficients */
    for (j = 0; j <= nb_iterators; j++)
      SCOPVAL_assign(compacted->p[i][j],matrix->p[i][j]);

    /* Then we copy the parameter coefficients */
    for (j = 0; j < nb_parameters; j++)
      SCOPVAL_assign(compacted->p[i][j + nb_iterators + 1],
		     matrix->p[i][j + CLAN_MAX_DEPTH + 1]);

    /* Lastly the scalar coefficient */
    SCOPVAL_assign(compacted->p[i][nb_columns - 1],
		   matrix->p[i][matrix->NbColumns - 1]);
  }

  scoplib_matrix_free_inside(matrix);

  /* Replace the inside of matrix */
  matrix->NbRows      = compacted->NbRows;
  matrix->NbColumns   = compacted->NbColumns;
  matrix->p           = compacted->p;
  matrix->p_Init      = compacted->p_Init;
  matrix->p_Init_size = compacted->p_Init_size;

  /* Free the compacted "container" */
  free(compacted);
}
