
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                   program.c                             **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: september 9th 2003               **
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
 * terms of the GNU Lesser General Public License as published by the Free    *
 * Software Foundation; either version 3 of the License, or (at your option)  *
 * any later version.                                                         *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with software; if not, write to the Free Software Foundation, Inc.,  *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * CAnDL, the Chunky Dependence Analyzer                                      *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <candl/candl.h>
#include <candl/program.h>


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/

/**
 * candl_program_print_structure function:
 * Displays a candl_program_t structure (program) into a file (file,
 * possibly stdout) in a way that trends to be understandable without falling
 * in a deep depression or, for the lucky ones, getting a headache... It
 * includes an indentation level (level) in order to work with others
 * print_structure functions.
 * - 09/09/2003: first version.
 */
void candl_program_print_structure(FILE* file, candl_program_p program,
				   int level)
{
  int i, j;

  if (program != NULL)
  {
    /* Go to the right level. */
    for (j = 0; j < level; j++)
      fprintf(file,"|\t");
    fprintf(file,"+-- candl_program_t\n");

    /* A blank line. */
    for (j = 0; j <= level + 1; j++)
      fprintf(file,"|\t");
    fprintf(file,"\n");

    /* Print the context. */
    candl_matrix_print_structure(file, program->context, level+1);

    /* A blank line. */
    for (j = 0; j <= level+1; j++)
      fprintf(file, "|\t");
    fprintf(file, "\n");

    /* Go to the right level and print the statement number. */
    for (j = 0; j <= level; j++)
      fprintf(file, "|\t");
    fprintf(file, "Statement number: %d\n", program->nb_statements);

    /* A blank line. */
    for (j = 0; j <= level+1; j++)
      fprintf(file, "|\t");
    fprintf(file, "\n");

    /* Print the statements. */
    for (i = 0; i < program->nb_statements; ++i)
      candl_statement_print_structure(file, program->statement[i], level+1);

    /* Print the transformation candidate. */
    if (program->transformation != NULL)
      for (i = 0; i < program->nb_statements; i++)
	candl_matrix_print_structure(file, program->transformation[i], level+1);
    else
      {
	/* Go to the right level. */
	for (j = 0; j <= level; j++)
	  fprintf(file, "|\t");
	fprintf(file, "+-- No transformation candidate\n");

	/* A blank line. */
	for (j = 0; j <= level+1; j++)
	  fprintf(file, "|\t");
	fprintf(file, "\n");
      }
  }
  else
    {
      /* Go to the right level. */
      for (j = 0; j < level; j++)
	fprintf(file, "|\t");
      fprintf(file, "+-- NULL candl_program_t\n");
    }

  /* The last line. */
  for (j = 0; j <= level; j++)
    fprintf(file, "|\t");
  fprintf(file, "\n");
}


/**
 * candl_program_print function:
 * This function prints the content of a candl_program_t structure
 * (program) into a file (file, possibly stdout).
 */
void candl_program_print(FILE * file,  candl_program_p program)
{
  candl_program_print_structure(file, program,0);
}



/**
 * candl_program_print function:
 * This function prints a candl_program_t structure (program) into a
 * candl-formatted file (file, possibly stdout).
 */
void candl_program_print_candl_file(FILE * file, candl_program_p program)
{
  int i, j;

  fprintf (file, "# -------------------\n");
  fprintf (file, "# Context\n");
  candl_matrix_print_data(file, program->context);
  fprintf (file, "\n");
  fprintf (file, "# Number of statements\n");
  fprintf (file, "%d\n", program->nb_statements);
  for (i = 0; i < program->nb_statements; ++i)
    {
      fprintf (file, "# -------------------\n");
      fprintf (file, "# Statement %d\n", i + 1);
      fprintf (file, "# Statement type\n");
      /* All types set to Assignment. */
      fprintf (file, "A\n");
      fprintf (file, "\n");
      fprintf (file, "# Iteration domain\n");
      candl_matrix_print_data(file, program->statement[i]->domain);
      fprintf (file, "\n");
      fprintf (file, "# Loop labels\n");
      for (j = 0; j < program->statement[i]->depth; ++j)
	fprintf (file, "%d ", program->statement[i]->index[j]);
      fprintf (file, "\n");
      fprintf (file, "# Written items\n");
      candl_matrix_print_data(file, program->statement[i]->written);
      fprintf (file, "\n");
      fprintf (file, "# Read items\n");
      candl_matrix_print_data(file, program->statement[i]->read);
      fprintf (file, "\n");
    }
  fprintf (file, "# -------------------\n");
  fprintf (file, "# Transformation candidate\n");
  fprintf (file, "0\n");
}


/******************************************************************************
 *                         Memory alloc/dealloc function                      *
 ******************************************************************************/


/**
 * candl_program_malloc function:
 * This function allocates the memory space for a candl_program_t structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - 09/12/2005: first version.
 */
candl_program_p candl_program_malloc()
{
  candl_program_p program;

  /* Memory allocation for the candl_program_t structure. */
  program = (candl_program_p)malloc(sizeof(candl_program_t));
  if (program == NULL)
    CANDL_FAIL("Error: memory overflow");

  /* We set the various fields with default values. */
  program->context        = NULL;
  program->nb_statements  = 0;
  program->statement      = NULL;
  program->transformation = NULL;
  program->scalars_privatizable = NULL;

  return program;
}


/**
 * candl_program_free function:
 * This function frees the allocated memory for a candl_program_t structure, it
 * recursively frees everything inside.
 */
void candl_program_free(candl_program_p program)
{
  int i;

  candl_matrix_free(program->context);

  if (program->statement != NULL)
    {
      for (i = 0; i < program->nb_statements; i++)
	candl_statement_free(program->statement[i]);
      free(program->statement);
    }

  if (program->transformation != NULL)
    {
      for (i = 0; i < program->nb_statements; i++)
	candl_matrix_free(program->transformation[i]);
      free(program->transformation);
    }

  if (program->scalars_privatizable != NULL)
    free(program->scalars_privatizable);

  free(program);
}


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/


/**
 * candl_program_read function:
 * This function reads the informations to put in a candl_program_t
 * structure from a file (file, possibly stdin). It returns a pointer
 * to a candl_program_t structure containing the read informations.
 * September 10th 2003: first version.
 */
candl_program_p candl_program_read(FILE * file)
{
  int i,  nb_statements,  nb_parameters,  nb_functions;
  char s[CANDL_MAX_STRING];
  CandlStatement ** statement;
  CandlMatrix    ** transformation;
  candl_program_p program;

  /* Memory allocation for the candl_program_t structure. */
  program = candl_program_malloc();

  /* First of all,  we read the context data. */
  program->context  = candl_matrix_read(file);
  nb_parameters = program->context->NbColumns - 2;

  /* We read the number of statements. */
  while (fgets(s, CANDL_MAX_STRING, file) == 0)
    ;
  while ((*s=='#'||*s=='\n') || (sscanf(s, " %d", &nb_statements) < 1))
    fgets(s, CANDL_MAX_STRING, file);

  program->nb_statements = nb_statements;

  /* Reading of each statement. */
  if (nb_statements > 0)
    {
      /* Memory allocation for the array of pointers to the statements. */
      statement = (CandlStatement**) malloc(nb_statements *
					    sizeof(CandlStatement*));
      if (statement == NULL)
	CANDL_FAIL("Error: memory overflow");

      for (i = 0; i < nb_statements; i++)
	statement[i] = candl_statement_read(file, i, nb_parameters);

      program->statement = statement;
    }

  /* We read the number of transformation functions. */
  while (fgets(s, CANDL_MAX_STRING, file) == 0)
    ;
  while ((*s=='#' || *s=='\n') || (sscanf(s, " %d", &nb_functions) < 1))
    fgets(s, CANDL_MAX_STRING, file);

  /* Reading of each transformation function. */
  if (nb_functions > 0)
    {
      /* The function number must be the same as statement number. */
      if (nb_functions != nb_statements)
	{
	  fprintf(stderr,
		  "[Candl]ERROR: the numbers of transformations (%d) and "
		  "statements (%d) differ.\n", nb_functions, nb_statements);
	  exit(1);
	}

      /* Memory allocation for the array of pointers to the functions. */
      transformation = (CandlMatrix **)malloc(nb_functions *
					      sizeof(CandlMatrix *));
      if (transformation == NULL)
	CANDL_FAIL("Error: memory overflow");

      for (i = 0; i < nb_functions; i++)
	transformation[i] = candl_matrix_read(file);

      program->transformation = transformation;
    }

  return(program);
}


/**
 * This function reads the .scop formatted file 'file', check for the
 * existence of the <candl> tag in the file, and retrieve the loop
 * index information, if any.
 * This function is built only if candl was configured with ScopLib support.
 *
 */
#ifdef CANDL_SUPPORTS_SCOPLIB
static
int** candl_program_scop_get_opt_indices(scoplib_scop_p scop)
{
  /* Get the <candl></candl> tag content. */
  char* candl_opts = scoplib_scop_tag_content(scop, "<candl>", "</candl>");
  if (! candl_opts)
    return NULL;
  /* Get the <candl><indices></indices></candl> tag content. */
  char* indices = scoplib_scop_tag_content_from_string(candl_opts, "<indices>",
						       "</indices>");
  free (candl_opts);
  if (! indices)
    return NULL;

  /* Tag was found. Scan it. */
  int buffer_size = 128;
  /* Assume maximum loop nest depth of 128. */
  int line[128];
  char buff[32];
  int** res = malloc(buffer_size * sizeof(int*));
  int i, j;
  int count, idx = 0;
  int line_added = 0;
  char* s = indices;

  while (s && *s != '\0')
    {
      for (i = 0; i < 128; ++i)
	{
	  while (*s != '\0' && *s != '\n' && isspace(*s))
	    ++s;
	  if (*s != '\0' && *s != '#' && *s != '\n')
	    {
	      for (count = 0; *s >= '0' && *s <= '9'; ++count)
		buff[count] = *(s++);
	      buff[count] = '\0';
	      line[i] = atoi(buff);
	      line_added = 1;
	    }
	  else
	    break;
	}
      if (line_added)
	{
	  if (idx == buffer_size)
	    res = realloc(res, (buffer_size *= 2) * sizeof(int*));
	  res[idx] = (int*) malloc(i * sizeof(int));
	  for (j = 0; j < i; ++j)
	    res[idx][j] = line[j];
	  ++idx;
	  line_added = 0;
	}
      while (s && *s != '\0' && *s != '\n')
	++s;
      if (s && *s != '\0' && *s == '\n')
	++s;
    }
  res = realloc(res, idx * sizeof(int*));
  free (indices);

  return res;
}
#endif


/**
 * candl_program_read_scop function:
 * This function reads the informations to put in a candl_program_t
 * structure from a file (file, possibly stdin) following the .scop
 * format.  It returns a pointer to a candl_program_t structure
 * containing the read informations.
 * This function is built only if candl was configured with ScopLib support.
 *
 */
#ifdef CANDL_SUPPORTS_SCOPLIB
candl_program_p candl_program_read_scop(FILE * file)
{
  int i;

  /* Read the scop. */
  scoplib_scop_p scop = scoplib_scop_read(file);
  /* Check for the <candl> tag in the options of the .scop file. */
  int** indices = candl_program_scop_get_opt_indices(scop);
  /* Convert the scop. */
  candl_program_p res = candl_program_convert_scop(scop, indices);

  /* Clean temporary data. */
  if (indices)
    {
      for (i = 0; i < res->nb_statements; ++i)
	free(indices[i]);
      free(indices);
    }
  scoplib_scop_free(scop);

  return res;
}
#endif


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * candl_program_convert_scop function:
 * This function extracts the useful information of a scoplib_scop_t
 * structure to a fresh, independent candl_program_t structure.
 * This function is built only if candl was configured with ScopLib support.
 *
 */
#ifdef CANDL_SUPPORTS_SCOPLIB
candl_program_p candl_program_convert_scop(scoplib_scop_p scop, int** indices)
{
  int i, j, k, l;
  candl_program_p res = candl_program_malloc();
  scoplib_statement_p s = scop->statement;

  /* Duplicate the context. */
  res->context = (CandlMatrix*) scoplib_matrix_copy(scop->context);
  if (res->context == NULL)
    res->context = candl_matrix_malloc(0, 2);

  /* Count the number of statements. */
  for (res->nb_statements = 0; s; s = s->next, res->nb_statements++)
    ;

  /* Allocate the statements array. */
  res->statement = (CandlStatement**) malloc(res->nb_statements *
					     sizeof(CandlStatement*));

  /* Initialize structures used in iterator indices computation. */
  int max = 0;
  int max_loop_depth = 128;
  int cur_index[max_loop_depth];
  int last[max_loop_depth];
  for (i = 0; i < max_loop_depth; ++i)
    {
      cur_index[i] = i;
      last[i] = 0;
    }
  /* Create the statements. */
  for (i = 0, s = scop->statement; s; s = s->next, ++i)
    {
      CandlStatement* statement = candl_statement_malloc();
      statement->label = i;
      statement->ref = s;
      if (s->domain->next != NULL)
	CANDL_FAIL("Error: union of domains not supported");

      statement->domain = (CandlMatrix*) scoplib_matrix_copy(s->domain->elt);
      /* For the moment, we do not parse the statement to extract its type. */
      statement->type = CANDL_ASSIGNMENT;
      statement->depth = statement->domain->NbColumns - 2 - scop->nb_parameters;
      statement->written = (CandlMatrix*) scoplib_matrix_copy(s->write);
      if (statement->written == NULL)
	statement->written =
	  candl_matrix_malloc(0, statement->domain->NbColumns);
      statement->read = (CandlMatrix*) scoplib_matrix_copy(s->read);
      if (statement->read == NULL)
	statement->read = candl_matrix_malloc(0, statement->domain->NbColumns);
      statement->index = (int*) malloc(statement->depth * sizeof(int));
      if (indices != NULL)
	/* Iterator indices are provided. */
	for (j = 0; j < statement->depth; ++j)
	  statement->index[j] = indices[i][j];
      else
	{
	  /* Iterator indices must be computed from the scattering matrix. */
	  scoplib_matrix_p m = s->schedule;
	  if (m == NULL)
	    CANDL_FAIL("Error: No scheduling matrix and no loop "
		       "indices specification");

	  /* FIXME: Sort the statements in their execution order. */
	  /* It must be a 2d+1 identity scheduling matrix, and
	     statement must be sorted in their execution order. */
	  /* Check that it is a identity matrix. */
	  int error = 0;
	  if (m->NbRows != 2 * statement->depth + 1)
	    error = 1;
	  else
	    for (l = 0; l < m->NbRows; ++l)
	      {
		for (k = 1; k < m->NbColumns - 1; ++k)
		  switch (CANDL_get_si(m->p[l][k]))
		    {
		    case 0:
		      if (l % 2 && k == (l / 2) + 1) error = 1;
		      break;
		    case 1:
		      if ((l % 2 && k != (l / 2) + 1) || (! l % 2)) error = 1;
		      break;
		    default:
		      error = 1;
		      break;
		    }
		if (l % 2 && CANDL_get_si(m->p[l][k]))
		  error = 1;
	      }
	  if (error)
	    CANDL_FAIL("Error: schedule is not identity 2d+1 shaped.\n"
		       "Consider using the <indices> option tag to declare "
		       " iterator indices");

	  /* Compute the value of the iterator indices. */
	  for (j = 0; j < statement->depth; ++j)
	    {
	      int val = CANDL_get_si(m->p[2 * j][m->NbColumns - 1]);
	      if (last[j] < val)
		{
		  last[j] = val;
		  for (k = j + 1; k < max_loop_depth; ++k)
		    last[k] = 0;
		  for (k = j; k < max_loop_depth; ++k)
		    cur_index[k] = max + (k - j) + 1;
		  break;
		}
	    }
	  for (j = 0; j < statement->depth; ++j)
	    statement->index[j] = cur_index[j];
	  max = max < cur_index[j - 1] ? cur_index[j - 1] : max;
	}
      /* Store the statement. */
      res->statement[i] = statement;
    }

  return res;
}
#endif
