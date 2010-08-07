
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 vector.c                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 01/05/2008                     **
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
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 2.1 of the     *
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
# include <ctype.h>
# include <scoplib/vector.h>


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * scoplib_vector_print_structure function:
 * Displays a scoplib_vector_t structure (*vector) into a file (file, possibly
 * stdout) in a way that trends to be understandable without falling in a deep
 * depression or, for the lucky ones, getting a headache... It includes an
 * indentation level (level) in order to work with others print_structure
 * functions.
 * \param file   File where informations are printed.
 * \param vector The vector whose information have to be printed.
 * \param level  Number of spaces before printing, for each line.
 **
 * - 01/05/2008: first version.
 */
void
scoplib_vector_print_structure(FILE * file, scoplib_vector_p vector, int level)
{
  int j;

  if (vector != NULL)
  {
    /* Go to the right level. */
    for (j = 0; j < level; j++)
      fprintf(file,"|\t");
    fprintf(file,"+-- scoplib_vector_t\n");

    for (j = 0; j <= level; j++)
      fprintf(file,"|\t");
    fprintf(file,"%d\n",vector->Size);

    /* Display the vector. */
    for (j = 0; j <= level; j++)
      fprintf(file,"|\t");

    fprintf(file,"[ ");

    for (j = 0; j < vector->Size; j++)
    {
      SCOPVAL_print(file,SCOPLIB_FMT,vector->p[j]);
      fprintf(file," ");
    }

    fprintf(file,"]\n");
  }
  else
  {
    /* Go to the right level. */
    for (j = 0; j < level; j++)
      fprintf(file,"|\t");
    fprintf(file,"+-- NULL vector\n");
  }

  /* The last line. */
  for (j = 0; j <= level; j++)
    fprintf(file,"|\t");
  fprintf(file,"\n");
}


/**
 * scoplib_vector_print function:
 * This function prints the content of a scoplib_vector_t structure
 * (*vector) into a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param vector The vector whose information have to be printed.
 **
 * - 01/05/2008: first version.
 */
void
scoplib_vector_print(FILE * file, scoplib_vector_p vector)
{
  scoplib_vector_print_structure(file,vector,0);
}


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/


/**
 * scoplib_vector_malloc function:
 * This function allocates the memory space for a scoplib_vector_t structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * \param Size The number of entries of the vector to allocate.
 **
 * - 01/05/2008: first version.
 */
scoplib_vector_p
scoplib_vector_malloc(unsigned Size)
{
  scoplib_vector_p vector;
  scoplib_int_t * p;
  int i ;

  vector = (scoplib_vector_p)malloc(sizeof(scoplib_vector_t));
  if (vector == NULL)
  {
    fprintf(stderr, "[Scoplib] Memory Overflow.\n");
    exit(1);
  }
  vector->Size = Size;
  if (Size == 0)
    vector->p = NULL;
  else
  {
    p = (scoplib_int_t *)malloc(Size * sizeof(scoplib_int_t));
    if (p == NULL)
    {
      fprintf(stderr, "[Scoplib] Memory Overflow.\n");
      exit(1);
    }
    vector->p = p;
    for (i = 0; i < Size; i++)
      SCOPVAL_init_set_si(vector->p[i],0);
  }
  return vector;
}


/**
 * scoplib_vector_free function:
 * This function frees the allocated memory for a scoplib_vector_t structure.
 * \param vector The pointer to the vector we want to free.
 **
 * - 01/05/2008: first version.
 */
void
scoplib_vector_free(scoplib_vector_p vector)
{
  int i;
  scoplib_int_t * p;

  if (vector != NULL)
  {
    p = vector->p;
    for (i = 0; i < vector->Size; i++)
      SCOPVAL_clear(*p++);

    free(vector->p);
    free(vector);
  }
}


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/

/**
 * scoplib_vector_add_scalar function:
 * This function adds a scalar to the vector representation of an affine
 * expression (this means we add the scalar only to the very last entry of the
 * vector). It returns a new vector resulting from this addition.
 * \param vector The basis vector.
 * \param scalar The scalar to add to the vector.
 **
 * - 01/05/2008: first version.
 */
scoplib_vector_p
scoplib_vector_add_scalar(scoplib_vector_p vector, int scalar)
{
  int i;
  scoplib_vector_p result;

  if ((vector == NULL) || (vector->Size < 2))
  {
    fprintf(stderr,"[Scoplib] Error: incompatible vector for addition\n");
    exit(1);
  }

  result = scoplib_vector_malloc(vector->Size);
  for (i = 0; i < vector->Size; i++)
    SCOPVAL_assign(result->p[i],vector->p[i]);
  SCOPVAL_add_int(result->p[vector->Size - 1],
		  vector->p[vector->Size - 1],scalar);

  return result;
}


/**
 * scoplib_vector_add function:
 * This function achieves the addition of two vectors and returns the
 * result as a new vector (the addition means the ith entry of the new vector
 * is equal to the ith entry of vector v1 plus the ith entry of vector v2).
 * \param v1 The first vector for the addition.
 * \param v2 The second vector for the addition (result is v1+v2).
 **
 * - 01/05/2008: first version.
 */
scoplib_vector_p
scoplib_vector_add(scoplib_vector_p v1, scoplib_vector_p v2)
{
  int i;
  scoplib_vector_p v3;

  if ((v1 == NULL) || (v2 == NULL) || (v1->Size != v2->Size))
  {
    fprintf(stderr,"[Scoplib] Error: incompatible vectors for addition\n");
    exit(1);
  }

  v3 = scoplib_vector_malloc(v1->Size);
  for (i = 0; i < v1->Size; i++)
    SCOPVAL_addto(v3->p[i],v1->p[i],v2->p[i]);

  return v3;
}


/**
 * scoplib_vector_sub function:
 * This function achieves the subtraction of two vectors and returns the
 * result as a new vector (the addition means the ith entry of the new vector
 * is equal to the ith entry of vector v1 minus the ith entry of vector v2).
 * \param v1 The first vector for the subtraction.
 * \param v2 The second vector for the subtraction (result is v1-v2).
 **
 * - 01/05/2008: first version.
 */
scoplib_vector_p
scoplib_vector_sub(scoplib_vector_p v1, scoplib_vector_p v2)
{
  int i;
  scoplib_vector_p v3;

  if ((v1 == NULL) || (v2 == NULL) || (v1->Size != v2->Size))
  {
    fprintf(stderr,"[Scoplib] Error: incompatible vectors for subtraction\n");
    exit(1);
  }

  v3 = scoplib_vector_malloc(v1->Size);
  for (i = 0; i < v1->Size; i++)
    SCOPVAL_subtract(v3->p[i],v1->p[i],v2->p[i]);

  return v3;
}


/**
 * scoplib_vector_tag_inequality function:
 * This function tags a vector representation of a contraint as being an
 * inequality >=0. This means in the PolyLib format, to set to 1 the very first
 * entry of the vector. It modifies directly the vector provided as an argument.
 * \param vector The vector to be tagged.
 **
 * - 01/05/2008: first version.
 */
void
scoplib_vector_tag_inequality(scoplib_vector_p vector)
{
  if ((vector == NULL) || (vector->Size < 1))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot be tagged\n");
    exit(1);
  }
  SCOPVAL_set_si(vector->p[0],1);
}


/**
 * scoplib_vector_tag_equality function:
 * This function tags a vector representation of a contraint as being an
 * equality ==0. This means in the PolyLib format, to set to 0 the very first
 * entry of the vector. It modifies directly the vector provided as an argument.
 * \param vector The vector to be tagged.
 **
 * - 01/05/2008: first version.
 */
void
scoplib_vector_tag_equality(scoplib_vector_p vector)
{
  if ((vector == NULL) || (vector->Size < 1))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot be tagged\n");
    exit(1);
  }
  SCOPVAL_set_si(vector->p[0],0);
}
