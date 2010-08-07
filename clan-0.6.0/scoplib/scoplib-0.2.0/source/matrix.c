
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
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
# include <string.h>
# include <ctype.h>
# include <scoplib/matrix.h>


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * scoplib_matrix_print_structure function:
 * Displays a scoplib_matrix_t structure (*matrix) into a file (file, possibly
 * stdout) in a way that trends to be understandable without falling in a deep
 * depression or, for the lucky ones, getting a headache... It includes an
 * indentation level (level) in order to work with others print_structure
 * functions.
 * \param file   File where informations are printed.
 * \param matrix The matrix whose information have to be printed.
 * \param level  Number of spaces before printing, for each line.
 **
 * - 30/04/2008: first version (from CLooG 0.14.0).
 */
void
scoplib_matrix_print_structure(FILE * file, scoplib_matrix_p matrix, int level)
{
  int i, j;

  /* Go to the right level. */
  for (j = 0; j < level; j++)
    fprintf(file,"|\t");

  if (matrix != NULL)
  {
    fprintf(file,"+-- scoplib_matrix_t\n");

    for(j = 0; j <= level; j++)
      fprintf(file,"|\t");
    fprintf(file,"%d %d\n",matrix->NbRows,matrix->NbColumns);

    /* Display the matrix. */
    for (i = 0; i < matrix->NbRows; i++)
    {
      for (j = 0; j <= level; j++)
        fprintf(file,"|\t");

      fprintf(file,"[ ");

      for (j = 0; j < matrix->NbColumns; j++)
      {
        SCOPVAL_print(file,SCOPLIB_FMT,matrix->p[i][j]);
        fprintf(file," ");
      }

      fprintf(file,"]\n");
    }
  }
  else
    fprintf(file,"+-- NULL matrix\n");

  /* The last line. */
  for (j = 0; j <= level; j++)
    fprintf(file,"|\t");
  fprintf(file,"\n");
}


/**
 * scoplib_matrix_print function:
 * This function prints the content of a scoplib_matrix_t structure
 * (*matrix) into a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param matrix The matrix whose information have to be printed.
 **
 * - 30/04/2008: first version (from CLooG 0.14.0).
 */
void
scoplib_matrix_print(FILE * file, scoplib_matrix_p matrix)
{
  scoplib_matrix_print_structure(file,matrix,0);
}



/**
 * scoplib_matrix_expression_element function:
 * This function returns a string containing the printing of a value (possibly
 * an iterator or a parameter with its coefficient or a constant).
 * \param val   The coefficient or constant value.
 * \param first Pointer to a boolean set to 1 if the current value is the first
 *              of an expresion, 0 otherwise (this function may update it).
 * \param cst   A boolean set to 1 if the value is a constant, 0 otherwise.
 * \param name  String containing the name of the iterator or of the parameter.
 **
 * - 03/05/2008: first version (from CLooG 0.14.0, glorious pprint_val).
 */
static
char *
scoplib_matrix_expression_element(scoplib_int_t val, int * first, int cst,
				  char * name)
{
  char * sval, * body, * temp;

  temp = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
  body = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
  sval = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
  body[0] = '\0';
  sval[0] = '\0';

  /* statements for the 'normal' processing. */
  if (SCOPVAL_notzero_p(val) && (!cst))
  {
    if ((*first) || SCOPVAL_neg_p(val))
    {
      if (SCOPVAL_one_p(val))           /* case 1 */
        sprintf(sval,"%s",name);
      else
      {
        if (SCOPVAL_mone_p(val))        /* case -1 */
          sprintf(sval,"-%s",name);
	else                         /* default case */
	{
	  SCOPVAL_sprint(sval,SCOPLIB_FMT_TXT,val);
	  sprintf(temp,"*%s",name);
	  strcat(sval,temp);
        }
      }
      *first = 0;
    }
    else
    {
      if (SCOPVAL_one_p(val))
        sprintf(sval,"+%s",name);
      else
      {
        sprintf(sval,"+");
	SCOPVAL_sprint(temp,SCOPLIB_FMT_TXT,val);
	strcat(sval,temp);
	sprintf(temp,"*%s",name);
	strcat(sval,temp);
      }
    }
  }
  else
  {
    if (cst)
    {
      if ((SCOPVAL_zero_p(val) && (*first)) || SCOPVAL_neg_p(val))
        SCOPVAL_sprint(sval,SCOPLIB_FMT_TXT,val);
      if (SCOPVAL_pos_p(val))
      {
        if (!(*first))
        {
	  SCOPVAL_sprint(sval,"+"SCOPLIB_FMT_TXT,val); /* Block macro ! */
	}
	else
          SCOPVAL_sprint(sval,SCOPLIB_FMT_TXT,val);
      }
    }
  }
  free(temp);
  free(body);

  return(sval);
}


/**
 * scoplib_matrix_expression function:
 * This function returns a string corresponding to an affine expression
 * stored at the "row"^th row of the matrix pointed by "matrix".
 * \param matrix        A set of linear expressions.
 * \param row           The row of the matrix corresponding to the expression.
 * \param nb_iterators  The number of iterators for the considered statement.
 * \param iterators     An array containing iterator names for the statement.
 * \param nb_parameters The number of parameters in the SCoP.
 * \param parameters    An array containing all parameters names.
 **
 * - 03/05/2008: first version (from CLooG 0.14.0, glorious pprint_val).
 */
static
char *
scoplib_matrix_expression(scoplib_matrix_p matrix, int row,
			  int nb_iterators,  char ** iterators,
			  int nb_parameters, char ** parameters)
{
  int i, first = 1;
  char * sline, * sval;

  sline = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char)) ;
  sline[0] = '\0' ;

  /* First the iterator part. */
  for (i = 1; i <= nb_iterators; i++)
  {
    sval = scoplib_matrix_expression_element(matrix->p[row][i],&first,0,
					     iterators[i-1]);
    strcat(sline,sval);
    free(sval);
  }

  /* Next the parameter part. */
  for (i = nb_iterators + 1; i <= nb_iterators + nb_parameters; i++)
  {
    sval = scoplib_matrix_expression_element(matrix->p[row][i],&first,0,
					     parameters[i - nb_iterators - 1]);
    strcat(sline,sval);
    free(sval);
  }

  /* Finally the constant part (yes, I reused it). */
  sval = scoplib_matrix_expression_element(matrix->p[row][i],&first,1,NULL);
  strcat(sline,sval);
  free(sval);

  return sline;
}


/**
 * scoplib_matrix_print_dot_scop function:
 * This function prints the content of a scoplib_matrix_t structure
 * (*matrix) into a file (file, possibly stdout) for the .scop format.
 * \param file          File where informations are printed.
 * \param matrix        The matrix whose information have to be printed.
 * \param type          A bit of semantic about this matrix (domain, access...).
 * \param nb_iterators  The number of iterators for the considered statement.
 * \param iterators     An array containing iterator names for the statement.
 * \param nb_parameters The number of parameters in the SCoP.
 * \param parameters    An array containing all parameters names.
 * \param nb_arrays     The number of arrays accessed in the SCoP.
 * \param arrays        An array containing all accessed array names.
 **
 * - 02/05/2008: first version.
 */
void
scoplib_matrix_print_dot_scop(FILE * file, scoplib_matrix_p matrix, int type,
			      int nb_iterators,  char ** iterators,
			      int nb_parameters, char ** parameters,
			      int nb_arrays,     char ** arrays)
{
  int i, j, k;
  char * expression;

  if (matrix == NULL)
  {
    fprintf(file,"0 %d\n",nb_iterators+nb_parameters+2);
    return;
  }

  fprintf(file,"%d %d\n",matrix->NbRows,matrix->NbColumns);

  for (i = 0; i < matrix->NbRows; i++)
  {
    for (j = 0; j < matrix->NbColumns; j++)
    {
      SCOPVAL_print(file,SCOPLIB_FMT,matrix->p[i][j]);
      fprintf(file," ");
    }

    if (type == SCOPLIB_TYPE_DOMAIN)
    {
      expression = scoplib_matrix_expression(matrix,i,nb_iterators,iterators,
					     nb_parameters,parameters);
      fprintf(file,"   ## %s",expression);
      free(expression);
      if (SCOPVAL_zero_p(matrix->p[i][0]))
        fprintf(file," == 0");
      else
        fprintf(file," >= 0");
    }

    if (type == SCOPLIB_TYPE_SCATTERING)
    {
      expression = scoplib_matrix_expression(matrix,i,nb_iterators,iterators,
					     nb_parameters,parameters);
      fprintf(file,"   ## %s",expression);
      free(expression);
    }

    if (type == SCOPLIB_TYPE_ACCESS)
    {
      if (SCOPVAL_notzero_p(matrix->p[i][0]))
      {
	if (strncmp(arrays[SCOPVAL_get_si(matrix->p[i][0]) - 1],
		    SCOPLIB_FAKE_ARRAY, strlen(SCOPLIB_FAKE_ARRAY)))
	  fprintf(file,"   ## %s",arrays[SCOPVAL_get_si(matrix->p[i][0]) - 1]);
	k = i;
	do
	{
	  expression = scoplib_matrix_expression(matrix,k,nb_iterators,
						 iterators,
						 nb_parameters,parameters);
          fprintf(file,"[%s]",expression);
          free(expression);
	  k++;
	}
	while ((k < matrix->NbRows) && SCOPVAL_zero_p(matrix->p[k][0]));
      }
      else
        fprintf(file,"   ##");
    }

    fprintf(file,"\n");
  }
}


/**
 * scoplib_matrix_list_print_structure function:
 * Displays a scoplib_matrix_list_t structure (a list of matrices) into a
 * file (file, possibly stdout). See scoplib_matrix_print_structure for
 * more details.
 * \param file   File where informations are printed.
 * \param l	 The list of matrices whose information have to be printed.
 * \param level  Number of spaces before printing, for each line.
 */
void
scoplib_matrix_list_print_structure(FILE * file, scoplib_matrix_list_p l,
				    int level)
{
  int j, first = 1;

  /* Go to the right level. */
  for (j = 0; j < level; j++)
    fprintf(file,"|\t");

  if (l != NULL)
    fprintf(file,"+-- scoplib_matrix_list_t\n");
  else
    fprintf(file,"+-- NULL matrix list\n");

  while (l != NULL)
  {
    if (!first)
    {
      /* Go to the right level. */
      for (j = 0; j < level; j++)
        fprintf(file,"|\t");
      fprintf(file,"|   scoplib_matrix_list_t\n");
    }
    else
      first = 0;

    /* A blank line. */
    for (j = 0; j <= level+1; j++)
      fprintf(file,"|\t");
    fprintf(file,"\n");

    /* Print a matrix. */
    scoplib_matrix_print_structure(file,l->elt,level+1);

    l = l->next;

    /* Next line. */
    if (l != NULL)
    {
      for (j = 0; j <= level; j++)
        fprintf(file,"|\t");
      fprintf(file,"V\n");
    }
  }

  /* The last line. */
  for (j = 0; j <= level; j++)
    fprintf(file,"|\t");
  fprintf(file,"\n");
}


/**
 * scoplib_matrix_list_print function:
 * This function prints the content of a scoplib_matrix_list_t into
 * a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param list The matrix whose information have to be printed.
 **
 * - 30/04/2008: first version (from CLooG 0.14.0).
 */
void
scoplib_matrix_list_print(FILE * file, scoplib_matrix_list_p list)
{
  scoplib_matrix_list_print_structure(file, list, 0);
}


/**
 * scoplib_matrix_list_print_dot_scop function:
 * This function prints the content of a scoplib_matrix_list_t structure into
 * a file (file, possibly stdout) for the .scop format.
 * \param file          File where informations are printed.
 * \param list          The matrix list whose information have to be printed.
 * \param type          A bit of semantic about this matrix (domain, access...).
 * \param nb_iterators  The number of iterators for the considered statement.
 * \param iterators     An array containing iterator names for the statement.
 * \param nb_parameters The number of parameters in the SCoP.
 * \param parameters    An array containing all parameters names.
 * \param nb_arrays     The number of arrays accessed in the SCoP.
 * \param arrays        An array containing all accessed array names.
 */
void
scoplib_matrix_list_print_dot_scop(FILE * file, scoplib_matrix_list_p list,
				   int type,
				   int nb_iterators,  char ** iterators,
				   int nb_parameters, char ** parameters,
				   int nb_arrays,     char ** arrays)
{
  int i;
  scoplib_matrix_list_p head = list;

  /* Count the number of elements in the list. */
  for (i = 0; list; list = list->next, i++)
    ;
  /* Print it. */
  fprintf(file,"%d\n", i);
  /* Print each element of the matrix list. */
  while (head)
  {
    scoplib_matrix_print_dot_scop(file, head->elt, type,
				  nb_iterators, iterators,
				  nb_parameters, parameters,
				  nb_arrays, arrays);
    head = head->next;
  }
}


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/


/**
 * scoplib_matrix_read function:
 * Adaptation from the PolyLib. This function reads a matrix into a file (foo,
 * posibly stdin) and returns a pointer this matrix.
 * October 18th 2001: first version.
 * - April 17th 2005: this function moved from domain.c to here.
 * - June  21rd 2005: Adaptation for GMP (based on S. Verdoolaege's version of
 *                    CLooG 0.12.1).
 * - July 9th 2008: Grabbed from CLooG and adapted for Scoplib.
 */
scoplib_matrix_p
scoplib_matrix_read(FILE* foo)
{
  unsigned NbRows, NbColumns;
  int i, j, n;
  char* c, s[SCOPLIB_MAX_STRING], str[SCOPLIB_MAX_STRING];
  scoplib_matrix_p matrix;
  scoplib_int_t* p;

  while (fgets(s, SCOPLIB_MAX_STRING, foo) == 0)
    ;
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, " %u %u", &NbRows, &NbColumns) < 2))
    fgets(s, SCOPLIB_MAX_STRING, foo);

  matrix = scoplib_matrix_malloc(NbRows, NbColumns);

  p = matrix->p_Init;
  for (i = 0; i < matrix->NbRows; i++)
    {
      do
	{
	  c = fgets(s, SCOPLIB_MAX_STRING, foo);
	  while ((c != NULL) && isspace(*c) && (*c != '\n'))
	    c++;
	}
      while (c != NULL && (*c == '#' || *c == '\n'));

      if (c == NULL)
	{
	  fprintf(stderr, "[Scoplib] Error: not enough rows\n");
	  exit(1);
	}
      for (j = 0; j < matrix->NbColumns; j++)
	{
	  if (c == NULL || *c == '#' || *c == '\n')
	    {
	      fprintf(stderr, "[Scoplib] Error: not enough columns\n");
	      exit(1);
	    }
	  if (sscanf(c, "%s%n", str, &n) == 0)
	    {
	      fprintf(stderr, "[Scoplib] Error: not enough rows\n");
	      exit(1);
	    }
#if defined(SCOPLIB_INT_T_IS_MP)
      long long val;
	  sscanf(str, "%lld", &val);
	  mpz_set_si(*p++, val);
#else
	  sscanf(str, SCOPLIB_FMT_TXT, p++);
#endif
	  c += n;
	}
    }

  return matrix;
}


/**
 * scoplib_matrix_list_read function:
 * This function reads a list of matrices into a file (foo,
 * posibly stdin) and returns a pointer this matrix list.
 * \param file   File where informations are stored.
 */
scoplib_matrix_list_p
scoplib_matrix_list_read(FILE* file)
{
  char s[SCOPLIB_MAX_STRING];
  int i;
  scoplib_matrix_list_p list;
  scoplib_matrix_list_p res;
  int nb_mat;

  /* Skip blank/commented lines. */
  while (fgets(s, SCOPLIB_MAX_STRING, file) == 0 || s[0] == '#' ||
	 isspace(s[0]))
    ;
  /* Read the number of matrices to read. */
  sscanf(s, "%d", &nb_mat);

  /* Allocate the header of the list. */
  res = list = scoplib_matrix_list_malloc();
  for (i = 0; i < nb_mat; ++i)
    {
      list->elt = scoplib_matrix_read(file);
      if (i < nb_mat - 1)
	list->next = scoplib_matrix_list_malloc();
      list = list->next;
    }

  return res;
}


/**
 * scoplib_matrix_read_arrays function:
 * This function reads a matrix into a file (foo, posibly stdin) and
 * returns a pointer this matrix. In addition, it reads the arrays as
 * comments at the end of the line.
 */
scoplib_matrix_p
scoplib_matrix_read_arrays(FILE* foo, char*** arrays, int* nb_arr)
{
  unsigned NbRows, NbColumns;
  int i, j, n;
  int count;
  char* c, s[SCOPLIB_MAX_STRING], str[SCOPLIB_MAX_STRING],
    buff[SCOPLIB_MAX_STRING];
  scoplib_matrix_p matrix;
  scoplib_int_t* p;

  while (fgets(s, SCOPLIB_MAX_STRING, foo) == 0)
    ;
  while ((*s=='#' || *s=='\n') ||
	 (sscanf(s, " %u %u", &NbRows, &NbColumns) < 2))
    fgets(s, SCOPLIB_MAX_STRING, foo);

  matrix = scoplib_matrix_malloc(NbRows, NbColumns);

  p = matrix->p_Init;
  for (i = 0; i < matrix->NbRows; i++)
    {
      do
	{
	  c = fgets(s, SCOPLIB_MAX_STRING, foo);
	  while ((c != NULL) && isspace(*c) && (*c != '\n'))
	    c++;
	}
      while (c != NULL && (*c == '#' || *c == '\n'));

      if (c == NULL)
	{
	  fprintf(stderr, "[Scoplib] Error: not enough rows\n");
	  exit(1);
	}

      for (j = 0; j < matrix->NbColumns; j++)
	{
	  if (c == NULL || *c == '#' || *c == '\n')
	    {
	      fprintf(stderr, "[Scoplib] Error: not enough columns\n");
	      exit(1);
	    }
	  if (sscanf(c, "%s%n", str, &n) == 0)
	    {
	      fprintf(stderr, "[Scoplib] Error: not enough rows\n");
	      exit(1);
	    }
#if defined(SCOPLIB_INT_T_IS_MP)
      long long val;
	  sscanf(str, "%lld", &val);
	  mpz_set_si(*p++, val);
#else
	  sscanf(str, SCOPLIB_FMT_TXT, p++);
#endif
	  c += n;
	}
      /* Read the array, passed as a comment at the end of the line. */
      if (c)
	{
	  while (c && (isspace(*c) || *c == '#'))
	    ++c;
	  for (count = 0; c && *c != '[' && *c != '\n'; ++count)
	    buff[count] = *(c++);
	  buff[count] = '\0';
	  if (count && SCOPVAL_get_si(matrix->p[i][0]))
	    {
	      /* Increase the buffer size if we run out of space. */
	      if (SCOPVAL_get_si(matrix->p[i][0]) - 1 > *nb_arr)
		{
		  *nb_arr = SCOPVAL_get_si(matrix->p[i][0]) - 1;
		  *arrays = (char**) realloc(*arrays,
					     sizeof(char*) * (*nb_arr + 1));
		}
	      /* Backup the array name. */
	      (*arrays)[SCOPVAL_get_si(matrix->p[i][0]) - 1] = strdup(buff);
	    }
	}
    }

  return matrix;
}


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/

/**
 * scoplib_matrix_malloc function:
 * This function allocates the memory space for a scoplib_matrix_t structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * \param NbRows    The number of row of the matrix to allocate.
 * \param NbColumns The number of columns of the matrix to allocate.
 **
 * - 30/04/2008: first version (from PipLib 1.4.0).
 */
scoplib_matrix_p
scoplib_matrix_malloc(unsigned NbRows, unsigned NbColumns)
{
  scoplib_matrix_p matrix;
  scoplib_int_t ** p, * q;
  int i, j;

  matrix = (scoplib_matrix_p)malloc(sizeof(scoplib_matrix_t));
  if (matrix == NULL)
  {
    fprintf(stderr, "[Scoplib] Memory Overflow.\n");
    exit(1);
  }
  matrix->NbRows = NbRows;
  matrix->NbColumns = NbColumns;
  matrix->p_Init_size = NbRows * NbColumns;
  if (matrix->p_Init_size == 0)
  {
    matrix->p = NULL;
    matrix->p_Init = NULL;
  }
  else
  {
    p = (scoplib_int_t **)malloc(NbRows*sizeof(scoplib_int_t *));
    if (p == NULL)
    {
      fprintf(stderr, "[Scoplib] Memory Overflow.\n");
      exit(1);
    }
    q = (scoplib_int_t *)malloc(NbRows * NbColumns * sizeof(scoplib_int_t));
    if (q == NULL)
    {
      fprintf(stderr, "[Scoplib] Memory Overflow.\n");
      exit(1);
    }
    matrix->p = p;
    matrix->p_Init = q;
    for (i = 0; i < NbRows; i++)
    {
      *p++ = q;
      for (j = 0; j < NbColumns; j++)
        SCOPVAL_init_set_si(*(q+j),0);
      q += NbColumns;
    }
  }
  return matrix;
}


/**
 * scoplib_matrix_free_inside function:
 * This function frees the allocated memory for the inside of a
 * scoplib_matrix_t structure, i.e. only p and p_Init.
 * \param matrix The pointer to the matrix we want to free.
 **
 * - 02/05/2008: first version.
 */
void
scoplib_matrix_free_inside(scoplib_matrix_p matrix)
{
  int i;
  scoplib_int_t * p;

  if (matrix == NULL)
    return;

  p = matrix->p_Init;
  for (i = 0; i < matrix->p_Init_size; i++)
    SCOPVAL_clear(*p++);

  if (matrix->p_Init != NULL)
    free(matrix->p_Init);

  if (matrix->p != NULL)
    free(matrix->p);
}


/**
 * scoplib_matrix_free function:
 * This function frees the allocated memory for a scoplib_matrix_t structure.
 * \param matrix The pointer to the matrix we want to free.
 **
 * - 30/04/2008: first version.
 * - 02/05/2008: now uses scoplib_matrix_free_inside.
 */
void
scoplib_matrix_free(scoplib_matrix_p matrix)
{
  if (matrix == NULL)
    return;

  scoplib_matrix_free_inside(matrix);
  free(matrix);
}


/**
 * scoplib_matrix_list_malloc function:
 * This function allocates the memory space for a scoplib_matrix_list_t
 * structure and sets its fields with default values. Then it returns
 * a pointer to the allocated space.
 */
scoplib_matrix_list_p
scoplib_matrix_list_malloc()
{
  scoplib_matrix_list_p res =
    (scoplib_matrix_list_p) malloc(sizeof(scoplib_matrix_list_t));

  if (res == NULL)
    {
      fprintf(stderr, "[Scoplib] Memory Overflow.\n");
      exit(1);
    }

  res->elt = NULL;
  res->next = NULL;

  return res;
}



/**
 * scoplib_matrix_list_free function:
 * This function frees the allocated memory for a scoplib_matrix_list_t
 * structure, and all the matrices stored in the list.
 * \param list The pointer to the matrix list we want to free.
 */
void
scoplib_matrix_list_free(scoplib_matrix_list_p list)
{
  scoplib_matrix_list_p tmp;

  if (list == NULL)
    return;

  while (list)
    {
      if (list->elt)
	scoplib_matrix_free(list->elt);
      tmp = list->next;
      free(list);
      list = tmp;
    }
}


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * scoplib_matrix_ncopy function:
 * this functions builds and returns a "hard copy" (not a pointer copy) of a
 * scoplib_matrix_t data structure such that the copy is restricted to the "n"
 * first rows of the matrix.
 * \param matrix The pointer to the matrix we want to copy.
 * \param n      The number of row of the matrix we want to copy.
 **
 * - 02/05/2008: first version.
 */
scoplib_matrix_p
scoplib_matrix_ncopy(scoplib_matrix_p matrix, int n)
{
  int i, j;
  scoplib_matrix_p copy;

  if (matrix == NULL)
    return NULL;

  if (n > matrix->NbRows)
  {
    fprintf(stderr,"[Scoplib] Error: not enough rows in the matrix\n");
    exit(1);
  }

  copy = scoplib_matrix_malloc(n,matrix->NbColumns);

  for (i = 0; i < n; i++)
    for (j = 0; j < matrix->NbColumns; j++)
      SCOPVAL_assign(copy->p[i][j],matrix->p[i][j]);

  return copy;
}


/**
 * scoplib_matrix_copy function:
 * this functions builds and returns a "hard copy" (not a pointer copy) of a
 * scoplib_matrix_t data structure.
 * \param matrix The pointer to the matrix we want to copy.
 **
 * - 30/04/2008: first version (from CLooG 0.14.0).
 * - 02/05/2008: no uses scoplib_matrix_ncopy.
 */
scoplib_matrix_p
scoplib_matrix_copy(scoplib_matrix_p matrix)
{
  if (matrix == NULL)
    return NULL;

  return scoplib_matrix_ncopy(matrix,matrix->NbRows);
}


/**
 * scoplib_matrix_replace_vector function:
 * this function replaces the "row"^th row of a matrix "matrix" with the
 * vector "vector". It directly updates "matrix".
 * \param matrix The matrix we want to change a row.
 * \param vector The vector that will replace a row of the matrix.
 * \param row    The row of the matrix to be replaced.
 **
 * - 02/05/2008: first version.
 */
void
scoplib_matrix_replace_vector(scoplib_matrix_p matrix, scoplib_vector_p vector,
			      int row)
{
  int i;

  if ((matrix == NULL) || (vector == NULL) ||
      (matrix->NbColumns != vector->Size) ||
      (row >= matrix->NbRows) || (row < 0))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot replace a matrix row\n");
    exit(1);
  }

  for (i = 0; i < vector->Size; i++)
    SCOPVAL_assign(matrix->p[row][i],vector->p[i]);
}


/**
 * scoplib_matrix_add_vector function:
 * this function adds the "row"^th row of a matrix "matrix" with the
 * vector "vector". It directly updates "matrix".
 * \param matrix The matrix we want to change a row.
 * \param vector The vector that will replace a row of the matrix.
 * \param row    The row of the matrix to be replaced.
 **
 * - 08/28/2008: first version.
 */
void
scoplib_matrix_add_vector(scoplib_matrix_p matrix, scoplib_vector_p vector,
			  int row)
{
  int i;

  if ((matrix == NULL) || (vector == NULL) ||
      (matrix->NbColumns != vector->Size) ||
      (row >= matrix->NbRows) || (row < 0))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot replace a matrix row\n");
    exit(1);
  }

  if (SCOPVAL_get_si(matrix->p[row][0]) == 0)
    SCOPVAL_assign(matrix->p[row][0],vector->p[0]);
  for (i = 1; i < vector->Size; i++)
    SCOPVAL_addto(matrix->p[row][i],matrix->p[row][i],vector->p[i]);
}


/**
 * scoplib_matrix_sub_vector function:
 * this function substracts the vector "vector" to the "row"^th row of
 * a matrix "matrix. It directly updates "matrix".
 * \param matrix The matrix we want to change a row.
 * \param vector The vector that will replace a row of the matrix.
 * \param row    The row of the matrix to be replaced.
 **
 * - 08/28/2008: first version.
 */
void
scoplib_matrix_sub_vector(scoplib_matrix_p matrix, scoplib_vector_p vector,
			  int row)
{
  int i;

  if ((matrix == NULL) || (vector == NULL) ||
      (matrix->NbColumns != vector->Size) ||
      (row >= matrix->NbRows) || (row < 0))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot replace a matrix row\n");
    exit(1);
  }

  if (SCOPVAL_get_si(matrix->p[row][0]) == 0)
    SCOPVAL_assign(matrix->p[row][0],vector->p[0]);
  for (i = 1; i < vector->Size; i++)
    SCOPVAL_subtract(matrix->p[row][i],matrix->p[row][i],vector->p[i]);
}


/**
 * scoplib_matrix_insert_vector function:
 * this function adds a new row corresponding to the vector "vector" to
 * the matrix "matrix" by inserting it at the "row"^th row. It directly
 * updates "matrix". If "vector" (or "matrix") is NULL, the matrix is left
 * unmodified.
 * \param matrix The matrix we want to extend.
 * \param vector The vector that will be added matrix.
 * \param row The row where to insert the vector.
 **
 * - 02/05/2008: first version.
 */
void
scoplib_matrix_insert_vector(scoplib_matrix_p matrix, scoplib_vector_p vector,
			     int row)
{
  int i, j;
  scoplib_matrix_p new;

  if ((vector == NULL) || (matrix == NULL))
    return;

  if ((matrix->NbColumns != vector->Size) ||
      (row > matrix->NbRows) || (row < 0))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot be inserted\n");
    exit(1);
  }

  /* We use a temporary matrix just to reuse existing functions. Cleaner. */
  new = scoplib_matrix_malloc(matrix->NbRows+1,matrix->NbColumns);

  for (i = 0; i < row; i++)
    for (j = 0; j < matrix->NbColumns; j++)
      SCOPVAL_assign(new->p[i][j],matrix->p[i][j]);

  scoplib_matrix_replace_vector(new,vector,row);

  for (i = row+1; i < matrix->NbRows; i++)
    for (j = 0; j < matrix->NbColumns; j++)
      SCOPVAL_assign(new->p[i][j],matrix->p[i-1][j]);

  scoplib_matrix_free_inside(matrix);

  /* Replace the inside of matrix */
  matrix->NbRows = new->NbRows;
  matrix->NbColumns = new->NbColumns;
  matrix->p = new->p;
  matrix->p_Init = new->p_Init;
  /* Free the new "shell" */
  free(new);
}


/**
 * scoplib_matrix_from_vector function:
 * this function converts a vector "vector" to a matrix with a single row
 * and returns a pointer to that matrix.
 * \param vector The vector to convert to a matrix.
 **
 * - 02/05/2008: first version.
 */
scoplib_matrix_p
scoplib_matrix_from_vector(scoplib_vector_p vector)
{
  scoplib_matrix_p matrix;

  if (vector == NULL)
    return NULL;

  matrix = scoplib_matrix_malloc(1,vector->Size);
  scoplib_matrix_replace_vector(matrix,vector,0);
  return matrix;
}


/**
 * scoplib_matrix_replace_matrix function:
 * this function replaces some rows of a matrix "m1" with the rows of
 * the matrix "m2". It begins at the "row"^th row of "m1". It directly
 * updates "m1".
 * \param m1  The matrix we want to change some row1.
 * \param m2  The matrix containing the new rows.
 * \param row The first row of the matrix m1 to be replaced.
 **
 * - 02/05/2008: first version.
 */
void
scoplib_matrix_replace_matrix(scoplib_matrix_p m1, scoplib_matrix_p m2, int row)
{
  int i, j;

  if ((m1 == NULL) || (m2 == NULL) ||
      (m1->NbColumns != m1->NbColumns) ||
      ((row + m2->NbRows) > m1->NbRows) || (row < 0))
  {
    fprintf(stderr,"[Scoplib] Error: vector cannot replace a matrix row\n");
    exit(1);
  }

  for (i = 0; i < m2->NbRows; i++)
    for (j = 0; j < m2->NbColumns; j++)
      SCOPVAL_assign(m1->p[i+row][j],m2->p[i][j]);
}


/**
 * scoplib_matrix_insert_matrix function:
 * this function adds new rows corresponding to the matrix "m1" to
 * the matrix "m2" by inserting it at the "row"^th row. It directly
 * updates "m1". If "m2" (or "m1") is NULL, the matrix is left
 * unmodified.
 * \param m1  The matrix we want to extend.
 * \param m2  The matrix to be inserted.
 * \param row The row where to insert the matrix
 **
 * - 02/05/2008: first version.
 */
void
scoplib_matrix_insert_matrix(scoplib_matrix_p m1, scoplib_matrix_p m2, int row)
{
  int i, j;
  scoplib_matrix_p new;

  if ((m1 == NULL) || (m2 == NULL))
    return;

  if ((m1->NbColumns != m2->NbColumns) ||
      (row > m1->NbRows) || (row < 0))
  {
    fprintf(stderr,"[Scoplib] Error: matrix cannot be inserted\n");
    exit(1);
  }

  /* We use a temporary matrix just to reuse existing functions. Cleaner. */
  new = scoplib_matrix_malloc(m1->NbRows+m2->NbRows,m1->NbColumns);

  for (i = 0; i < row; i++)
    for (j = 0; j < m1->NbColumns; j++)
      SCOPVAL_assign(new->p[i][j],m1->p[i][j]);

  scoplib_matrix_replace_matrix(new,m2,row);

  for (i = row + m2->NbRows; i < m1->NbRows; i++)
    for (j = 0; j < m1->NbColumns; j++)
      SCOPVAL_assign(new->p[i][j],m1->p[i-m2->NbRows][j]);

  scoplib_matrix_free_inside(m1);

  /* Replace the inside of matrix */
  m1->NbRows = new->NbRows;
  m1->NbColumns = new->NbColumns;
  m1->p = new->p;
  m1->p_Init = new->p_Init;

  /* Free the new "container" */
  free(new);
}


/**
 * scoplib_matrix_concat function:
 * this function builds a new matrix as the concatenation of the rows of
 * two other matrices sent as parameters.
 * \param m1  The first matrix.
 * \param m2  The second matrix.
 **
 * - 02/05/2008: first version.
 */
scoplib_matrix_p
scoplib_matrix_concat(scoplib_matrix_p m1, scoplib_matrix_p m2)
{
  scoplib_matrix_p new;

  if (m1 == NULL)
    return scoplib_matrix_copy(m2);

  if (m2 == NULL)
    return scoplib_matrix_copy(m1);

  if (m1->NbColumns != m2->NbColumns)
  {
    fprintf(stderr,"[Scoplib] Error: matrices cannot be concatenated\n");
    exit(1);
  }

  new = scoplib_matrix_malloc(m1->NbRows+m2->NbRows,m1->NbColumns);
  scoplib_matrix_replace_matrix(new,m1,0);
  scoplib_matrix_replace_matrix(new,m2,m1->NbRows);

  return new;
}


/**
 * scoplib_matrix_equal function:
 * this function returns true if the two matrices are the same, false
 * otherwise.
 * \param m1  The first matrix.
 * \param m2  The second matrix.
 *
 */
int
scoplib_matrix_equal(scoplib_matrix_p m1, scoplib_matrix_p m2)
{
  if (m1 == m2)
    return 1;
  if (m1 == NULL || m2 == NULL)
    return 0;
  if (m1->NbRows != m2->NbRows || m1->NbColumns != m2->NbColumns)
    return 0;
  int i, j;
  for (i = 0; i < m1->NbRows; ++i)
    for (j = 0; j < m1->NbColumns; ++j)
      if (SCOPVAL_ne(m1->p[i][j], m2->p[i][j]))
	return 0;
  return 1;
}
