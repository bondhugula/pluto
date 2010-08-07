
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                  dependence.c                           **
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
/**
 * \file dependence.c
 * \author Cedric Bastoul and Louis-Noel Pouchet
 */

# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <candl/candl.h>

#include <assert.h>


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * candl_dependence_print_structure function:
 * Displays a CandlDependence structure (dependence) into a file (file,
 * possibly stdout) in a way that trends to be understandable without falling
 * in a deep depression or, for the lucky ones, getting a headache... It
 * includes an indentation level (level) in order to work with others
 * print_structure functions.
 * - 18/09/2003: first version.
 */
void candl_dependence_print_structure(FILE* file,
				      candl_dependence_p dependence,
				      int level)
{
  int j, first = 1;

  if (dependence != NULL)
    { /* Go to the right level. */
      for(j=0; j<level; j++)
	fprintf(file, "|\t");
      fprintf(file, "+-- CandlDependence\n");
    }
  else
    { for(j=0; j<level; j++)
	fprintf(file, "|\t");
      fprintf(file, "+-- NULL dependence\n");
    }

  while (dependence != NULL)
    { if (!first)
	{ /* Go to the right level. */
	  for(j=0; j<level; j++)
	    fprintf(file, "|\t");
	  fprintf(file, "|   CandlDependence\n");
	}
      else
	first = 0;

      /* A blank line. */
      for(j=0; j<=level+1; j++)
	fprintf(file, "|\t");
      fprintf(file, "\n");

      /* Go to the right level and print the type. */
      for(j=0; j<=level; j++)
	fprintf(file, "|\t");
      fprintf(file, "Type: ");
      switch (dependence->type)
	{
	case CANDL_UNSET : fprintf(file, "UNSET\n");        break;
	case CANDL_RAW   : fprintf(file, "RAW (flow)\n");   break;
	case CANDL_WAR   : fprintf(file, "WAR (anti)\n");   break;
	case CANDL_WAW   : fprintf(file, "WAW (output)\n"); break;
	case CANDL_RAR   : fprintf(file, "RAR (input)\n");  break;
	default : fprintf(file, "unknown\n"); break;
	}

      /* A blank line. */
      for(j=0; j<=level+1; j++)
	fprintf(file, "|\t");
      fprintf(file, "\n");

      /* Go to the right level and print the depth. */
      for(j=0; j<=level; j++)
	fprintf(file, "|\t");
      fprintf(file, "Depth: %d\n", dependence->depth);

      /* A blank line. */
      for(j=0; j<=level+1; j++)
	fprintf(file, "|\t");
      fprintf(file, "\n");

      /* Print the source statement. */
      candl_statement_print_structure(file, dependence->source, level+1);

      /* Print the target statement. */
      candl_statement_print_structure(file, dependence->target, level+1);

      /* Print the dependence polyhedron. */
      candl_matrix_print_structure(file, dependence->domain, level+1);

      dependence = dependence->next;

      /* Next line. */
      if (dependence != NULL)
	{ for (j=0; j<=level; j++)
	    fprintf(file, "|\t");
	  fprintf(file, "V\n");
	}
    }

  /* The last line. */
  for(j=0; j<=level; j++)
    fprintf(file, "|\t");
  fprintf(file, "\n");
}


/* candl_dependence_print function:
 * This function prints the content of a CandlDependence structure (dependence)
 * into a file (file, possibly stdout).
 */
void candl_dependence_print(FILE * file, candl_dependence_p dependence)
{
  candl_dependence_print_structure(file, dependence, 0);
}


/* candl_dependence_pprint function:
 * This function prints the content of a CandlDependence structure (dependence)
 * into a file (file, possibly stdout) as a Graphviz input file.
 * See http://www.graphviz.org
 * - 08/12/2005: first version.
 */
void candl_dependence_pprint(FILE * file, candl_dependence_p dependence)
{
  int i = 0;

  fprintf(file, "digraph G {\n");

  fprintf(file, "# Data Dependence Graph\n");
  fprintf(file, "# Generated by Candl "CANDL_RELEASE" "CANDL_VERSION" bits\n");
  while (dependence != NULL)
    {
      fprintf(file, "  S%d -> S%d [label=\" ", dependence->source->label,
	      dependence->target->label);
      switch (dependence->type)
	{
	case CANDL_UNSET : fprintf(file, "UNSET"); break;
	case CANDL_RAW   : fprintf(file, "RAW")  ; break;
	case CANDL_WAR   : fprintf(file, "WAR")  ; break;
	case CANDL_WAW   : fprintf(file, "WAW")  ; break;
	case CANDL_RAR   : fprintf(file, "RAR")  ; break;
	default : fprintf(file, "unknown"); break;
	}
      fprintf(file, " depth %d, ref %d->%d \"];\n", dependence->depth,
	      dependence->ref_source,
	      dependence->ref_target);
      dependence = dependence->next;
      i++;
    }

  if (i>4)
    fprintf(file, "# Number of edges = %i\n}\n", i);
  else
    fprintf(file, "}\n");
}


/* candl_dependence_view function:
 * This function uses dot (see http://www.graphviz.org) and gv (see
 * http://wwwthep.physik.uni-mainz.de/~plass/gv) tools to display the
 * dependence graph.
 * - 20/03/2006: first version.
 */
void candl_dependence_view(candl_dependence_p dependence)
{
  FILE * temp_output;

  temp_output = fopen(CANDL_TEMP_OUTPUT, "w");
  candl_dependence_pprint(temp_output, dependence);
  fclose(temp_output);
  system("(dot -Tps "CANDL_TEMP_OUTPUT" | gv - &) && rm -f "CANDL_TEMP_OUTPUT);
}


/**
 * Returns a string containing the dependence, formatted to fit the
 * .scop representation.
 *
 */
static
char*
candl_program_deps_to_string(CandlDependence* dependence)
{
  CandlDependence* tmp = dependence;
  int refs = 0, reft = 0;
  int i, j, k;
  int nb_deps;
  int buffer_size = 2048;
  int szbuff;
  char* buffer = (char*) malloc(buffer_size * sizeof(char));
  char buff[1024];
  char* type;
  char* pbuffer;

  for (tmp = dependence, nb_deps = 0; tmp; tmp = tmp->next, ++nb_deps)
    ;
  sprintf(buffer, "# Number of dependences\n%d\n", nb_deps);
  if (nb_deps)
    {
      for (tmp = dependence, nb_deps = 1; tmp; tmp = tmp->next, ++nb_deps)
	{
	  /* Compute the type of the dependence, and the array id
	     accessed. */
	  switch (tmp->type)
	    {
	    case CANDL_UNSET:
	      type = "UNSET";
	      break;
	    case CANDL_RAW:
	      type = "RAW #(flow)";
	      refs = CANDL_get_si(tmp->source->written->p[tmp->ref_source][0]);
	      reft = CANDL_get_si(tmp->target->read->p[tmp->ref_target][0]);
	      break;
	    case CANDL_WAR:
	      type = "WAR #(anti)";
	      refs = CANDL_get_si(tmp->source->read->p[tmp->ref_source][0]);
	      reft = CANDL_get_si(tmp->target->written->p[tmp->ref_target][0]);
	      break;
	    case CANDL_WAW:
	      type = "WAW #(output)";
	      refs = CANDL_get_si(tmp->source->written->p[tmp->ref_source][0]);
	      reft = CANDL_get_si(tmp->target->written->p[tmp->ref_target][0]);
	      break;
	    case CANDL_RAR:
	      type = "RAR #(input)";
	      refs = CANDL_get_si(tmp->source->read->p[tmp->ref_source][0]);
	      reft = CANDL_get_si(tmp->target->read->p[tmp->ref_target][0]);
	      break;
	    default:
	      exit(1);
	      break;
	    }
	  /* Quick consistency check. */
	  if (refs != reft)
	    CANDL_FAIL("Internal error. refs != reft\n");

	  /* Output dependence information. */
	  sprintf(buff, "# Description of dependence %d\n"
		  "# type\n%s\n# From statement id\n%d\n"
		  "# To statement id\n%d\n# Depth \n%d\n# Array id\n%d\n"
		  "# Dependence domain\n%d %d\n", nb_deps, type,
		  tmp->source->label, tmp->target->label, tmp->depth,
		  refs, tmp->domain->NbRows, tmp->domain->NbColumns);
	  strcat(buffer, buff);
	  /* Output dependence domain. */
	  pbuffer = buffer + strlen(buffer);
	  for (i = 0; i < tmp->domain->NbRows; ++i)
	    {
	      for (j = 0; j < tmp->domain->NbColumns; ++j)
		{
		  sprintf(buff, "%lld ", CANDL_get_si(tmp->domain->p[i][j]));
		  szbuff = strlen(buff);
		  if (szbuff == 2)
		    *(pbuffer++) = ' ';
		  for (k = 0; k < szbuff; ++k)
		    *(pbuffer++) = buff[k];
		}
	      *(pbuffer++) = '\n';
	    }
	  *(pbuffer++) = '\0';
	  /* Increase the buffer size if needed. Conservatively assume a
	     dependence is never larger than 2k. */
	  szbuff = strlen(buffer);
	  if (szbuff + 2048 > buffer_size)
	    {
	      buffer = (char*) realloc(buffer, (buffer_size *= 2) *
				       sizeof(char));
	      if (buffer == NULL)
		CANDL_FAIL("Error: memory overflow");
	      buffer[szbuff] = '\0';
	    }
	}
    }

  return buffer;
}


/**
 * candl_dependence_print_scop function:
 * This function adds to the .scop file provided as the 'input' the
 * optional tags to represent the dependences 'dependence' of the
 * program. Finally, it prints the updated .scop to the file 'output'.
 *
 *
 */
#ifdef CANDL_SUPPORTS_SCOPLIB
/**
 * Read one dependence from a string.
 *
 */
static
CandlDependence* candl_dependence_read_one_dep(char* str, char** next,
					       CandlProgram* program)
{
  CandlDependence* dep = candl_dependence_malloc();
  CandlMatrix* msource = NULL;
  CandlMatrix* mtarget = NULL;
  char buffer[1024];

  int i, j, k;
  int id;
  int id2;
  /* # Description of dependence xxx */
  for (; *str != '\n'; ++str);
  ++str;

  /* # type */
  for (; *str != '\n'; ++str);
  ++str;

  /* {RAW,RAR,WAR,WAW} #(type) */
  for (i = 0; *str != ' '; ++str, ++i)
    buffer[i] = *str;
  buffer[i] = '\0';
  if (! strcmp(buffer, "RAW"))
    dep->type = CANDL_RAW;
  else if (! strcmp(buffer, "RAR"))
    dep->type = CANDL_RAR;
  else if (! strcmp(buffer, "WAR"))
    dep->type = CANDL_WAR;
  else if (! strcmp(buffer, "WAW"))
    dep->type = CANDL_WAW;
  for (; *str != '\n'; ++str);
  ++str;

  /* # From statement xxx */
  for (; *str != '\n'; ++str);
  ++str;
  /* stmt-id */
  for (i = 0; *str != '\n'; ++str, ++i)
    buffer[i] = *str;
  ++str;
  buffer[i] = '\0';
  id = atoi(buffer);
  for (i = 0; i < program->nb_statements; ++i)
    if (program->statement[i]->label == id)
      {
	dep->source = program->statement[i];
	break;
      }

  /* # To statement xxx */
  for (; *str != '\n'; ++str);
  ++str;
  /* stmt-id */
  for (i = 0; *str != '\n'; ++str, ++i)
    buffer[i] = *str;
  ++str;
  buffer[i] = '\0';
  id = atoi(buffer);
  for (i = 0; i < program->nb_statements; ++i)
    if (program->statement[i]->label == id)
      {
	dep->target = program->statement[i];
	break;
      }

  /* # Depth */
  for (; *str != '\n'; ++str);
  ++str;
  /* depth */
  for (i = 0; *str != '\n'; ++str, ++i)
    buffer[i] = *str;
  ++str;
  buffer[i] = '\0';
  dep->depth = atoi(buffer);

  /* # Array id */
  for (; *str != '\n'; ++str);
  ++str;
  /* array-id */
  for (i = 0; *str != '\n'; ++str, ++i)
    buffer[i] = *str;
  ++str;
  buffer[i] = '\0';
  id = atoi(buffer);
  switch (dep->type)
    {
    case CANDL_RAW:
      msource = dep->source->written;
      mtarget = dep->target->read;
      break;
    case CANDL_WAR:
      msource = dep->source->read;
      mtarget = dep->target->written;
      break;
    case CANDL_WAW:
      msource = dep->source->written;
      mtarget = dep->target->written;
      break;
    case CANDL_RAR:
      msource = dep->source->read;
      mtarget = dep->target->read;
      break;
    default:
      exit(1);
      break;
    }
  for (i = 0; i < msource->NbRows && msource->p[i][0] != id; ++i)
    ;
  if (i < msource->NbRows)
    dep->ref_source = i;
  for (i = 0; i < mtarget->NbRows && mtarget->p[i][0] != id; ++i)
    ;
  if (i < mtarget->NbRows)
    dep->ref_target = i;

  /* # Dependence domain */
  for (; *str != '\n'; ++str);
  ++str;

  /* nb-row nb-col */
  for (i = 0; *str != ' '; ++str, ++i)
    buffer[i] = *str;
  ++str;
  buffer[i] = '\0';
  id = atoi(buffer);
  for (i = 0; *str != '\n'; ++str, ++i)
    buffer[i] = *str;
  ++str;
  buffer[i] = '\0';
  id2 = atoi(buffer);

  dep->domain = candl_matrix_malloc(id, id2);
  /* Read matrix elements. */
  for (j = 0; j < id; ++j)
    {
    for (k = 0; k < id2; ++k)
      {
	while (*str && *str == ' ')
	  str++;
	for (i = 0; *str != '\n' && *str != ' '; ++str, ++i)
	  buffer[i] = *str;
	buffer[i] = '\0';
	++str;
	CANDL_set_si(dep->domain->p[j][k], atoi(buffer));
      }
    if (*(str - 1) != '\n')
      {
	for (; *str != '\n'; ++str);
	++str;
      }
    }
  /* Store the next character in the string to be analyzed. */
  *next = str;

  return dep;
}

/**
 * Retrieve a CandlDependence list from the option tag in the scop.
 *
 */
CandlDependence* candl_dependence_read_from_scop(scoplib_scop_p scop,
						 CandlProgram* program)
{
  CandlDependence* first = NULL;
  CandlDependence* currdep = NULL;

  char* deps = scoplib_scop_tag_content(scop,
					"<dependence-polyhedra>",
					"</dependence-polyhedra>");
  /* No dependence, nothing to do. */
  if (deps == NULL)
    return NULL;

  /* Keep the starting address of the array. */
  char* base = deps;

  int i;
  int depcount;
  /* Get the number of dependences. */
  char buffer_nb[32];
  /* # Number of dependences */
  for (; *deps != '\n'; ++deps);
  ++deps;
  for (i = 0; *deps != '\n'; ++i, ++deps)
    buffer_nb[i] = deps[i];
  buffer_nb[i] = '\0';
  ++deps;
  int nbdeps = atoi (buffer_nb);
  char* next;
  /* For each of them, read 1 and shift of the read size. */
  for (depcount = 0; depcount < nbdeps; ++depcount)
    {
      CandlDependence* adep =
	candl_dependence_read_one_dep(deps, &next, program);
      if (first == NULL)
	currdep = first = adep;
      else
	{
	  currdep->next = adep;
	  currdep = currdep->next;
	}
      deps = next;
    }

  /* Be clean. */
  free(base);

  return first;
}

/**
 * Update the scop option tag with the dependence list.
 *
 */
void candl_dependence_update_scop_with_deps(scoplib_scop_p scop,
					    CandlDependence* dependence)
{
  char* start;
  char* stop;
  char* content;
  char* olddeps = NULL;
  char* newdeps;
  char* newopttags;
  char* curr;
  char* tmp;
  int size = 0;
  int size_newdeps;
  int size_olddeps = 0;
  int size_optiontags;

  start = stop = scop->optiontags;
  /* Get the candl tag, if any. */
  content = scoplib_scop_tag_content(scop, "<candl>", "</candl>");
  if (content)
    {
      /* Get the dependence tag, if any. */
      olddeps = scoplib_scop_tag_content_from_string
	(content, "<dependence-polyhedra>", "</dependence-polyhedra>");
      /* Seek for the correct start/stop characters to insert
	 dependences. */
      if (olddeps)
	{
	  size = size_olddeps = strlen(olddeps);
	  while (start && *start && strncmp(start, olddeps, size))
	    ++start;
	  stop = start + size;
	}
      else
	{
	  size = strlen(content);
	  while (start && *start && strncmp(start, content, size))
	    ++start;
	  stop = start;
	}
    }

  /* Convert the program dependences to dotscop representation. */
  newdeps = candl_program_deps_to_string(dependence);

  /* Compute the new size of the full options tags, and allocate a new
     string. */
  size_newdeps = newdeps ? strlen(newdeps) : 0;
  size_optiontags = scop->optiontags ? strlen(scop->optiontags) : 0;
  if (content == NULL)
    size = strlen("<candl>") + strlen("</candl>") +
      strlen("<dependence-polyhedra>")
      + strlen("</dependence-polyhedra>");
  else if (olddeps == NULL)
    size = strlen("<dependence-polyhedra>") + strlen("</dependence-polyhedra>");
  else
    size = 0;
  newopttags = (char*) malloc((size_newdeps + size_optiontags
			       - size_olddeps + size + 1)
			      * sizeof(char));
  if (newopttags == NULL)
    CANDL_FAIL("Error: memory overflow");
  curr = newopttags;

  /* Copy the beginning of the options. */
  for (tmp = scop->optiontags; tmp != start; ++tmp)
    *(curr++) = *tmp;
  *curr = '\0';
  /* Copy the candl tags, if needed. */
  if (content == NULL)
    {
      strcat(newopttags, "<candl>\n");
      curr += strlen("<candl>\n");
    }
  if (olddeps == NULL)
    {
      strcat(newopttags, "<dependence-polyhedra>\n");
      curr += strlen("<dependence-polyhedra>\n");
    }
  /* Copy the program dependences. */
  for (tmp = newdeps; tmp && *tmp; ++tmp)
      *(curr++) = *tmp;
  *curr = '\0';
  /* Copy the candl tags, if needed. */
  if (olddeps == NULL)
    {
      strcat(curr, "</dependence-polyhedra>\n");
      curr += strlen("</dependence-polyhedra>\n");
    }
  if (content == NULL)
    {
      strcat(curr, "</candl>\n");
      curr += strlen("</candl>\n");
    }
  /* Copy the remainder of the options. */
  for (tmp = stop; tmp && *tmp; ++tmp)
      *(curr++) = *tmp;
  *curr = '\0';

  if (scop->optiontags)
    free(scop->optiontags);
  scop->optiontags = newopttags;

  /* Be clean. */
  if (content)
    free(content);
  if (olddeps)
    free(olddeps);
  if (newdeps)
    free(newdeps);
}

/**
 * Print the scop, containing the list of dependences.
 *
 */
void candl_dependence_print_scop(FILE* input, FILE* output,
				 CandlDependence* dependence)
{
  scoplib_scop_p scop;

  /* Go to the beginning of the file. */
  rewind(input);

  /* Re-read the options tags. */
  scop = scoplib_scop_read(input);

  /* Embed the dependences in the scop option tag. */
  candl_dependence_update_scop_with_deps(scop, dependence);

  /* Dump the .scop. */
  scoplib_scop_print_dot_scop(output, scop);
}
#endif


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/* candl_dependence_free function:
 * This function frees the allocated memory for a CandlDependence structure.
 * - 18/09/2003: first version.
 */
void candl_dependence_free(candl_dependence_p dependence)
{
  candl_dependence_p next;

  while (dependence != NULL)
    {
      next = dependence->next;
      candl_matrix_free(dependence->domain);
      free(dependence);
      dependence = next;
    }
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * candl_dependence_malloc function:
 * This function allocates the memory space for a CandlDependence structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - 07/12/2005: first version.
 */
candl_dependence_p candl_dependence_malloc()
{
  candl_dependence_p dependence;

  /* Memory allocation for the CandlDependence structure. */
  dependence = (candl_dependence_p) malloc(sizeof(CandlDependence));
  if (dependence == NULL)
    CANDL_FAIL(" Error: memory overflow");

  /* We set the various fields with default values. */
  dependence->source     = NULL;
  dependence->target     = NULL;
  dependence->depth      = CANDL_UNSET;
  dependence->type       = CANDL_UNSET;
  dependence->ref_source = CANDL_UNSET;
  dependence->ref_target = CANDL_UNSET;
  dependence->domain     = NULL;
  dependence->next       = NULL;
  dependence->usr	 = NULL;

  return dependence;
}


/**
 * candl_dependence_add function:
 * This function adds a CandlDependence structure (dependence) at a given place
 * (now) of a NULL terminated list of CandlDependence structures. The beginning
 * of this list is (start). This function updates (now) to the end of the loop
 * list (loop), and updates (start) if the added element is the first one -that
 * is when (start) is NULL-.
 * - 18/09/2003: first version.
 */
void candl_dependence_add(candl_dependence_p* start,
			  candl_dependence_p* now,
			  candl_dependence_p dependence)
{
  if (dependence != NULL)
    {
      if (*start == NULL)
	{
	  *start = dependence;
	  *now = *start;
	}
      else
	{
	  (*now)->next = dependence;
	  *now = (*now)->next;
	}

      while ((*now)->next != NULL)
	*now = (*now)->next;
    }
}


/**
 * GCD computation.
 */
static
int
candl_dependence_gcd(int a, int b)
{
  int z = 1;

  if (a < 0)
    a *= -1;
  if (b < 0)
    b *= -1;
  if (a == 0)
    return b;
  if (b == 0)
    return a;
  if (b > a)
    {
      int temp = a;
      a = b;
      b = temp;
    }

  while (z != 0)
    {
      z = a % b;
      a = b;
      b = z;
    }

  return a;
}

/**
 *
 *
 */
static
int candl_dependence_gcd_test_context (CandlMatrix* system, int id)
{
  /* FIXME: implement me! */

  return 1;
}


/**
 * candl_dependence_gcd_test function:
 * This functions performs a GCD test on a dependence polyhedra
 * represented exactly by a set of constraints 'system' organized in
 * such a way:
 * - first lines: iteration domain of 'source'
 * - then: iteration domain of 'target'
 * - then: array access equality(ies)
 * - then (optional): precedence constraint inequality.
 *
 * The function returns false if the dependence is impossible, true
 * otherwise. A series of simple checks (SIV/ZIV/MIV/bounds checking)
 * are also performed before the actual GCD test.
 *
 */
int candl_dependence_gcd_test(CandlStatement* source,
			      CandlStatement* target,
			      CandlMatrix* system,
			      int level)
{
  int i;
  int gcd;
  int id;
  int value;
  int null_iter, null_param, null_cst, pos_iter, neg_iter, strict_pred;

  /* Check that the precedence constraint, if any, is not strict in a
     self-dependence. */
  if (source == target &&
      CANDL_get_si(system->p[system->NbRows - 1][0]) == 1 &&
      CANDL_get_si(system->p[system->NbRows - 1][system->NbColumns - 1]) == -1)
    strict_pred = 1;
  else
    strict_pred = 0;

  /* Inspect the array access function equalities. */
  for (id = source->domain->NbRows + target->domain->NbRows;
       id < system->NbRows && CANDL_get_si(system->p[id][0]) == 0; ++id)
    {
      /* Inspect which parts of the access function equality are null,
	 positive or negative. */
      null_iter = null_param = null_cst = pos_iter = neg_iter = 0;
      for (i = 1; i < source->depth + target->depth + 1 &&
	     CANDL_get_si(system->p[id][i]) == 0; ++i)
	;
      if (i == source->depth + target->depth + 1)
	null_iter = 1;
      else
	for (pos_iter = 1, neg_iter = 1;
	     i < source->depth + target->depth + 1; ++i)
	  {
	    if (CANDL_get_si(system->p[id][i]) < 0)
	      pos_iter = 0;
	    else if (CANDL_get_si(system->p[id][i]) > 0)
	      neg_iter = 0;
	  }
      for (; i < system->NbColumns - 1 && CANDL_get_si(system->p[id][i]) == 0;
	   ++i)
	;
      if (i == system->NbColumns - 1)
	null_param = 1;
      null_cst = ! CANDL_get_si(system->p[id][system->NbColumns - 1]);

      /* Some useful ZIV/SIV/MIV tests. */
      if (null_iter && null_param && !null_cst)
	return 0;
      if (null_iter)
	if (! candl_dependence_gcd_test_context (system, id))
	  return 0;
      if (null_cst || !null_param)
	continue;
/* FIXME: implement the loop bound check. */
/*       /\* A clever test on access bounds. *\/ */
/*       if (null_param && pos_iter &&  */
/* 	  CANDL_get_si(system->p[id][system->NbColumns - 1]) > 0) */
/* 	return 0; */
/*       if (null_param && neg_iter &&  */
/* 	  CANDL_get_si(system->p[id][system->NbColumns - 1]) < 0) */
/* 	return 0; */

      /* Compute GCD test for the array access equality. */
      for (i = 1, gcd = CANDL_get_si(system->p[id][i]);
	   i < source->depth + target->depth; ++i)
	gcd = candl_dependence_gcd(gcd, CANDL_get_si(system->p[id][i + 1]));
      value = CANDL_get_si(system->p[id][system->NbColumns - 1]);
      value = value < 0 ? -value : value;
      if ((gcd == 0 && value != 0) || value % gcd)
	return 0;
    }

  return 1;
}


/**
 * candl_dependence_build_system function:
 * this function builds the constraint system corresponding to a data
 * dependence, for a given statement couple (with iteration domains "source"
 * and "target"), for a given reference couple (the source reference is array
 * "ref_s" in "array_s" and the target reference is the array "ref_t" in
 * "array_t"), at a given depth "depth" and knowing if the source is textually
 * before the target (boolean "before"). The system is built... as always !
 * See the [bastoul and Feautrier, PPL 2005] paper for details !
 * - source is the source iteration domain,
 * - target is the target iteration domain,
 * - array_s is the array list for the source,
 * - array_t is the array list for the target,
 * - ref_s is the position of the source reference in array_s,
 * - ref_s is the position of the target reference in array_t,
 * - depth is the dependence depth,
 * - before is 1 if the source is textually before the target, 0 otherwise,
 * - nb_par is the number of parameters.
 ***
 * - 13/12/2005: first version (extracted from candl_dependence_system).
 * - 23/02/2006: a stupid bug corrected in the subscript equality.
 * - 07/04/2007: fix the precedence condition to respect C. Bastoul PhD thesis
 */
static
CandlMatrix * candl_dependence_build_system(source, target, array_s, array_t,
					    ref_s, ref_t, depth, before, nb_par)
CandlStatement  * source, * target;
CandlMatrix * array_s, * array_t ;
int ref_s, ref_t, depth, before, nb_par ;
{ int i, j, nb_rows, nb_columns, nb_dimensions, constraint, s_dims, t_dims ;
  CandlMatrix * system ;
  Entier temp ;

  CANDL_init(temp) ;

  /* We calculate the number of dimensions of the considered array. */
  nb_dimensions = 1;
  while (((ref_s + nb_dimensions + 1) <= array_s->NbRows) &&
         (array_s->p[ref_s + nb_dimensions][0] == 0))
  nb_dimensions ++ ;

  /* The number of dimensions of the source and target domains. */
  s_dims = source->domain->NbColumns - nb_par - 2 ;
  t_dims = target->domain->NbColumns - nb_par - 2 ;

  /* The number of rows of the system is:
   * - the number of constraints of the source iteration domain +
   * - the number of constraints of the target iteration domain +
   * - the number of dimensions of the considered array (subscript equality) +
   * - the number of precedence constraints (equal to depth).
   */
  nb_rows = source->domain->NbRows + target->domain->NbRows +
    nb_dimensions + depth ;

  /* The number of columns of the system is:
   * - the number of source statement surrounding loops +
   * - the number of target statement surrounding loops +
   * - the number of parameters +
   * - 2 (1 for equality/inequality identifier + 1 for the constant).
   */
  nb_columns = s_dims + t_dims + nb_par + 2 ;

  /* We allocate memory space for the constraint system. */
  system = candl_matrix_malloc(nb_rows,nb_columns) ;

  /* Compute the maximal common depth. */
  int min_depth = 0;
  while (min_depth < source->depth && min_depth < target->depth &&
	 source->index[min_depth] == target->index[min_depth])
    ++min_depth;

  /* We fill the constraint system (note that there is no need to put zeros
   * in the empty zones since candl_matrix_alloc initialized all entries to 0):
   */

  /* 1. The source iteration domain constraints. */
  constraint = 0 ;
  for (i=0;i<source->domain->NbRows;i++)
  { for (j=0;j<=s_dims;j++)
    CANDL_assign(system->p[constraint][j],source->domain->p[i][j]) ;

    for (j=s_dims+1;j<source->domain->NbColumns;j++)
    CANDL_assign(system->p[constraint][j+t_dims],source->domain->p[i][j]) ;

    constraint++ ;
  }

  /* 2. The target iteration domain constraints. */
  for (i = 0; i < target->domain->NbRows; i++)
  { CANDL_assign(system->p[constraint][0],target->domain->p[i][0]) ;

    for (j = 1; j < target->domain->NbColumns; j++)
    CANDL_assign(system->p[constraint][j + s_dims], target->domain->p[i][j]) ;

    constraint++ ;
  }

  int subeq = 0;

  /* 3. The equality of the subscripts. */
  for (i = 0; i < nb_dimensions; i++)
  { /* Source iterator coefficients part. */
    for (j = 1; j <= s_dims; j++)
    CANDL_assign(system->p[constraint][j], array_s->p[ref_s + i][j]) ;

    /* Target iterator coefficients part (negative). */
    for (j = 1; j <= t_dims; j++)
    { CANDL_oppose(temp, array_t->p[ref_t + i][j]) ;
      CANDL_assign(system->p[constraint][j + s_dims], temp) ;
    }

    /* Parameters and constant coefficients part. */
    for (j = 1; j <= nb_par + 1; j++)
    CANDL_subtract(system->p[constraint][j + s_dims + t_dims],
                    array_s->p[ref_s + i][j + s_dims],
                    array_t->p[ref_t + i][j + t_dims]) ;
    constraint++ ;
    subeq++;
  }

  /* 4. The precedence constraints (their number is equal to depth). */
  for (i = 0; i < depth; i++)
  { /* i = i' for all dimension less than depth. */
    CANDL_set_si(system->p[constraint][i + 1], -1) ;
    CANDL_set_si(system->p[constraint][s_dims + i + 1], 1) ;
    if (i == depth - 1)
      {
	/* i <= i' at dimension depth if source is textually before target. */
	CANDL_set_si(system->p[constraint][0], 1) ;
	/* If source is textually after target, this is obviously i < i'. */
	if (before || depth < min_depth)
	  //if (before)
	  CANDL_set_si(system->p[constraint][nb_columns - 1], -1);
      }

    constraint++ ;
  }

  CANDL_clear(temp) ;
  return system ;
}


/**
 * candl_dependence_system function :
 * this function builds a node of the dependence graph: it studies a dependence
 * between a given statement couple, reference couple, type and depth. If a
 * data dependence actually exists, it returns a dependence structure, it
 * returns NULL otherwise.
 * - source is the source statement,
 * - target is the target statement,
 * - context is the program context (contraints on global parameters),
 * - array_s is the array list for the source,
 * - array_t is the array list for the target,
 * - ref_s is the position of the source reference in array_s,
 * - ref_s is the position of the target reference in array_t,
 * - depth is the dependence depth,
 * - type is the dependence type (RAW, WAW, WAR or RAR).
 ***
 * - 18/09/2003: first version.
 * - 09/12/2005: few corrections and polishing.
 */
candl_dependence_p candl_dependence_system(CandlStatement* source,
					   CandlStatement* target,
					   CandlMatrix* context,
					   CandlMatrix* array_s,
					   CandlMatrix* array_t,
					   int ref_s, int ref_t,
					   int type, int depth)
{
  candl_dependence_p dependence = NULL;
  CandlMatrix * system;
  PipOptions * options;
  PipQuast * solution;

  /* First, a trivial case: for two different statements at depth 0, there is
   * a dependence only if the source is textually before the target.
   */
  if ((source != target) && (depth == 0) && (source->label > target->label))
    return NULL;

  /* We build the system of constraints. */
  system = candl_dependence_build_system(source,  target,
					 array_s, array_t, ref_s, ref_t,
					 depth,
					 (source->label >= target->label),
					 context->NbColumns-2);

  /* We start by simple SIV/ZIV/GCD tests. */
  if (! candl_dependence_gcd_test(source, target, system, depth))
    return NULL;

  options = pip_options_init();
  options->Simplify = 1;
  options->Urs_parms = -1;
  options->Urs_unknowns = -1;
  solution = pip_solve(system,context, -1, options);

  if ((solution != NULL) &&
      ((solution->list != NULL) || (solution->condition != NULL)))
    {
      dependence = candl_dependence_malloc();

      /* We set the various fields with corresponding values. */
      dependence->type       = type;
      dependence->depth      = depth;
      dependence->source     = source;
      dependence->target     = target;
      dependence->ref_source = ref_s;
      dependence->ref_target = ref_t;
      dependence->domain     = system;
    }
  else {
    candl_matrix_free(system);
  }
  pip_options_free(options);
  pip_quast_free(solution);

  return dependence;
}


/**
 * candl_dependence_between function :
 * this function builds the dependence list from the statement "source" to
 * statement "target": it will study the dependence for each reference and for
 * each depth, under a particular context (context) and according to some
 * user options.
 * - 18/09/2003: first version.
 * - 07/12/2005: (debug) correction of depth bounds.
 * - 09/12/2005: We may take commutativity into consideration.
 */
candl_dependence_p candl_dependence_between(CandlStatement* source,
					    CandlStatement* target,
					    CandlMatrix* context,
					    CandlOptions* options)
{
  int i, j, k, min_depth, max_depth;
  candl_dependence_p new;
  candl_dependence_p dependence = NULL;
  candl_dependence_p now;

  /* If the statements commute and the user asks to use this information to
   * simplify the dependence graph, we return no dependences.
   */
  if (options->commute && candl_statement_commute(source, target))
    return NULL;

  /* In the case of a self-dependence, the dependence depth can be as low as 1
   * (not 0 because at depth 0 there is no loop, thus there is only one
   * instance of the statement !) and as high as the statement depth.
   * In the case of different statements, the dependence depth can be as low
   * as 0 and as high as the number of shared loops.
   */
  if (source == target)
    {
      min_depth = 1;
      max_depth = source->depth;
    }
  else
    {
      /* Depth 0 is for statements that don't share any loop. */
      min_depth = (source->index[0] == target->index[0]) ? 1 : 0;

      max_depth = 0;
      while ((max_depth < source->depth) &&
	     (max_depth < target->depth) &&
	     (source->index[max_depth] == target->index[max_depth]))
	max_depth++;
    }

  /* Flow and output-dependences analysis. */
  for (j = 0; j < source->written->NbRows; j++)
    if (CANDL_notzero_p(source->written->p[j][0]))
      {
	/* Flow-dependences analysis. */
	if (options->raw)
	  for (k = 0; k < target->read->NbRows; k++)
	    if (CANDL_eq(target->read->p[k][0], source->written->p[j][0]))
	      for (i = min_depth; i <= max_depth; i++)
		{
		  new = candl_dependence_system(source, target, context,
						source->written, target->read,
						j, k, CANDL_RAW, i);
		  candl_dependence_add(&dependence, &now, new);
		}
	/* Output-dependences analysis. */
	if (options->waw)
	  for (k = 0; k < target->written->NbRows; k++)
	    if (CANDL_eq(target->written->p[k][0], source->written->p[j][0]))
	      for (i = min_depth; i <= max_depth; i++)
		{
		  new = candl_dependence_system(source, target, context,
						source->written,
						target->written, j, k,
						CANDL_WAW, i);
		  candl_dependence_add(&dependence, &now, new);
		}
      }

  /* Anti and input-dependences analysis. */
  for (j = 0; j < source->read->NbRows; j++)
    if (source->read->p[j][0] != 0)
      {
	/* Anti-dependences analysis. */
	if (options->war)
	  for (k = 0; k < target->written->NbRows; k++)
	    if (CANDL_eq(target->written->p[k][0], source->read->p[j][0]))
	      for (i = min_depth; i <= max_depth; i++)
		{
		  new = candl_dependence_system(source, target, context,
						source->read, target->written,
						j, k, CANDL_WAR, i);
		  candl_dependence_add(&dependence, &now, new);
		}
	/* Input-dependences analysis. */
	if (options->rar)
	  for (k = 0; k < target->read->NbRows; k++)
	    if (CANDL_eq(target->read->p[k][0], source->read->p[j][0]))
	      for (i = min_depth; i <= max_depth; i++)
		{
		  new = candl_dependence_system(source, target, context,
						source->read, target->read,
						j, k, CANDL_RAR, i);
		  candl_dependence_add(&dependence, &now, new);
		}
      }

  return dependence;
}





/**
 * candl_dependence function:
 * this function builds the dependence graph of a program (program)
 * according to some user options (options).
 * - 18/09/2003: first version.
 */
candl_dependence_p candl_dependence(candl_program_p program,
				    CandlOptions* options)
{
  int i, j;
  candl_dependence_p dependence = NULL;
  candl_dependence_p new = NULL;
  candl_dependence_p now;
  CandlStatement  ** statement;
  CandlMatrix * context;

  statement = program->statement;
  context = program->context;

  if (options->scalar_privatization || options->scalar_expansion)
    candl_dependence_analyze_scalars (program, options);

  for (i = 0; i < program->nb_statements; i++)
    { /* We add self dependence. */
      /* S->S */
      new = candl_dependence_between(statement[i], statement[i],
				     context, options);
      candl_dependence_add(&dependence, &now, new);

      for (j = i + 1; j < program->nb_statements; j++)
	{ /* And dependences with other statements. */
	  /* S1->S2 */
	  new = candl_dependence_between(statement[i], statement[j],
					 context, options);
	  candl_dependence_add(&dependence, &now, new);

	  /* S2->S1 */
	  new = candl_dependence_between(statement[j], statement[i],
					 context, options);
	  candl_dependence_add(&dependence, &now, new);
	}
    }

  /* If scalar analysis is called, remove some useless dependences. */
  if (options->scalar_privatization || options->scalar_expansion ||
      options->scalar_renaming)
    candl_dependence_prune_scalar_waw (program, options, &dependence);

  /* Final treatment for scalar analysis. */
  int check = 0;
  if (options->scalar_renaming)
    check = candl_dependence_scalar_renaming (program, options, &dependence);
  if (! check && options->scalar_privatization)
    candl_dependence_prune_with_privatization (program, options, &dependence);

  /* Compute the last writer */
  if (options->lastwriter)  {
      candl_compute_last_writer(dependence, program);
  }

  return dependence;
}


/******************************************************************************
 *                          Scalar analysis functions                         *
 ******************************************************************************/

/**
 * candl_dependence_var_is_scalar function:
 * This function returns true if the variable indexed by 'var_index'
 * is a scalar in the whole program.
 */
int
candl_dependence_var_is_scalar (candl_program_p program, int var_index)
{
  CandlMatrix* m;
  int i, j, k, cpt;

  for (i = 0; i < program->nb_statements; ++i)
    for (m = program->statement[i]->read, cpt = 0; cpt < 2; ++cpt,
	   m = program->statement[i]->written)
      for (j = 0; j < m->NbRows; ++j)
	if (CANDL_get_si(m->p[j][0]) == var_index)
	  {
	    /* Ensure it is not an array. */
	    if (j < m->NbRows - 1 && CANDL_get_si(m->p[j + 1][0]) == 0)
	      return 0;
	    /* Ensure the access function is '0'. */
	    for (k = 1; k < m->NbColumns; ++k)
	      if (CANDL_get_si(m->p[j][k]) != 0)
		return 0;
	  }

  return 1;
}


/**
 * candl_dependence_expand_scalar function:
 * Expand the variable of index 'scalar_idx' by adding one
 * dimension (x becomes x[0]) to each access matrix refering this
 * variable in the statement list.
 *
 */
static
void
candl_dependence_expand_scalar(CandlStatement** sl,
			       int scalar_idx)
{
  CandlMatrix* m;
  int i, l, n, j;

  /* Iterate on all statements of the list. */
  for (i = 0; sl[i] != NULL; ++i)
    {
      /* Check if the scalar is referenced in the 'read' access
	 function. */
      for (j = 0; j < sl[i]->read->NbRows &&
	     CANDL_get_si(sl[i]->read->p[j][0]) != scalar_idx; ++j)
	;
      /* It is. */
      if (j < sl[i]->read->NbRows)
	{
	  /* Add a row to the 'read' matrix, just after the reference
	     to 'scalar_idx'. */
	  m = candl_matrix_malloc (sl[i]->read->NbRows +1,
				   sl[i]->read->NbColumns);
	  for (l = 0; l <= j; ++l)
	    for (n = 0; n < m->NbColumns; ++n)
	      CANDL_assign(m->p[l][n], sl[i]->read->p[l][n]);
	  for (++l; l < m->NbRows; ++l)
	    for (n = 0; n < m->NbColumns; ++n)
	      CANDL_assign(m->p[l][n], sl[i]->read->p[l - 1][n]);
	  for (n = 0; n < m->NbColumns; ++n)
	    CANDL_set_si(m->p[j + 1][n], 0);
	  candl_matrix_free (sl[i]->read);
	  sl[i]->read = m;
	}

      /* Same for 'written' access function. */
      for (j = 0; j < sl[i]->written->NbRows &&
	     CANDL_get_si(sl[i]->written->p[j][0]) != scalar_idx;++j)
	;
      if (j < sl[i]->written->NbRows)
	{
	  m = candl_matrix_malloc (sl[i]->written->NbRows +1,
				   sl[i]->written->NbColumns);
	  for (l = 0; l <= j; ++l)
	    for (n = 0; n < m->NbColumns; ++n)
	      CANDL_assign(m->p[l][n], sl[i]->written->p[l][n]);
	  for (++l; l < m->NbRows; ++l)
	    for (n = 0; n < m->NbColumns; ++n)
	      CANDL_assign(m->p[l][n], sl[i]->written->p[l - 1][n]);
	  for (n = 0; n < m->NbColumns; ++n)
	    CANDL_set_si(m->p[j + 1][n], 0);
	  candl_matrix_free (sl[i]->written);
	  sl[i]->written = m;
	}
    }
}


/**
 * candl_dependence_refvar_chain function:
 * This function returns a chain of statements as a feshly allocated
 * array of pointer on statements, of all statements reading or
 * writing the variable 'var_index', surrounded by the 'level' common
 * loops of 'dom'.  Output is a NULL-terminated array.
 */
CandlStatement**
candl_dependence_refvar_chain(candl_program_p program, CandlStatement* dom,
			      int var_index, int level)
{
  /* No or empty program -> no chain! */
  if (program == NULL || program->nb_statements == 0)
    return NULL;

  int buffer_size = 64;
  CandlStatement** res =
    (CandlStatement**) malloc(buffer_size * sizeof(CandlStatement*));
  int i, j, count = 0;
  CandlStatement* s;

  /* If no dominator is provided, assume we start with the first statement. */
  if (dom == NULL)
    dom = program->statement[0];
  for (i = 0; i < program->nb_statements && program->statement[i] != dom; ++i)
    ;
  /* The dominator is not in the list of statements. */
  if (i == program->nb_statements)
    return NULL;
  for (; i < program->nb_statements; ++i)
    {
      s = program->statement[i];
      /* We look for exactly 'level' common loops. */
      if (s->depth < level)
	continue;
      /* Ensure it has 'level' common loop(s) with the dominator. */
      for (j = 0; j < level&& s->index[j] == dom->index[j]; ++j)
	;
      if (j < level)
	continue;
      /* Ensure the variable is referenced. */
      if (candl_dependence_var_is_ref (s, var_index) != CANDL_VAR_UNDEF)
	{
	  res[count++] = s;
	  if (count == buffer_size)
	    res = realloc(res, (buffer_size*=2) * sizeof(CandlStatement*));
	}
    }

  res = realloc(res, (count + 1) * sizeof(CandlStatement*));
  res[count] = NULL;

  return res;
}


/**
 * candl_dependence_var_is_ref function:
 * This function checks if a var 'var_index' is referenced (DEF or
 * USE) by the statement.
 */
int
candl_dependence_var_is_ref(CandlStatement* s, int var_index)
{
  int j;
  int ret = CANDL_VAR_UNDEF;

  if (s)
    {
      for (j = 0; s->read && j < s->read->NbRows; ++j)
	if (CANDL_get_si(s->read->p[j][0]) == var_index)
	  {
	    ret = CANDL_VAR_IS_USED;
	    break;
	  }
      for (j = 0; s->written && j < s->written->NbRows; ++j)
	if (CANDL_get_si(s->written->p[j][0]) == var_index)
	  {
	    if (ret == CANDL_VAR_IS_USED)
	      ret = CANDL_VAR_IS_DEF_USED;
	    else
	      ret = CANDL_VAR_IS_DEF;
	    break;
	  }
    }

  return ret;
}


/**
 * candl_dependence_compute_lb function:
 * This function assigns to the Entier 'lb' the lexmin of variable
 * 'col'-1 in the polyhedron 'm'.
 */
static
void
candl_dependence_compute_lb (CandlMatrix* m, Entier* lb, int col)
{
  PipOptions* options;
  PipQuast* solution;
  PipList* l;
  options = pip_options_init ();
  options->Simplify = 1;
  options->Urs_parms = -1;
  options->Urs_unknowns = -1;
  /* Compute lexmin. */
  solution = pip_solve (m, NULL, -1, options);
  if ((solution != NULL) &&
      ((solution->list != NULL) || (solution->condition != NULL)))
    {
      for (l = solution->list; l && col > 0; l = l->next, --col)
	;
      CANDL_assign(*lb, l->vector->the_vector[0]);
    }
  pip_options_free (options);
  pip_quast_free (solution);
}


/**
 * candl_dependence_check_domain_is_included function:
 * This function checks if the 'level'-first iterators of 2 domains
 * are in such a way that s1 is larger or equal to s2, for the
 * considered iterator dimensions.
 *
 */
int
candl_dependence_check_domain_is_included(CandlStatement* s1,
					  CandlStatement* s2,
					  CandlMatrix* context,
					  int level)
{
  int max = level;
  Entier lb; CANDL_init(lb);
  max = s1->depth < max ? s1->depth : max;
  max = s2->depth < max ? s2->depth : max;

  CandlMatrix* m = candl_matrix_malloc(s2->domain->NbRows + s2->depth - max +1,
				       s2->domain->NbColumns);
  int i, j;
  /* Duplicate s2 to the dest matrix. */
  for (i = 0; i < s2->domain->NbRows; ++i)
    for (j = 0; j < s2->domain->NbColumns; ++j)
      CANDL_assign(m->p[i][j], s2->domain->p[i][j]);
  /* Make useless dimensions equal to 1. */
  for (j = 0; j < s2->depth - max; ++j)
    {
      candl_dependence_compute_lb (s2->domain, &lb, j + 1 + max);
      CANDL_assign(m->p[i][m->NbColumns - 1], lb);
      CANDL_set_si(m->p[i++][j + 1 + max], -1);
    }
  /* Iterate on all constraints of s1, and check them. */
  for (i = 0; i < s1->domain->NbRows; ++i)
    {
      /* Skip constraints defining other iterators. */
      for (j = max + 1; j <= s1->depth; ++j)
	if (CANDL_get_si(s1->domain->p[i][j]) != 0)
	  break;
      if (j <= s1->depth)
	continue;
      /* Invert the constraint, and add it to m. */
      for (j = 0; j <= max; ++j)
	{
	  CANDL_assign(m->p[m->NbRows - 1][j], s1->domain->p[i][j]);
	  CANDL_oppose(m->p[m->NbRows - 1][j], m->p[m->NbRows - 1][j]);
	}
      for (j = s1->depth + 1; j < s1->domain->NbColumns; ++j)
	{
	  CANDL_assign(m->p[m->NbRows - 1][j - s1->depth + s2->depth],
		       s1->domain->p[i][j]);
	  CANDL_oppose(m->p[m->NbRows - 1][j - s1->depth + s2->depth],
		       m->p[m->NbRows - 1][j - s1->depth + s2->depth]);
	}
      /* Make the inequality strict. */
      CANDL_decrement(m->p[m->NbRows - 1][m->NbColumns - 1],
		      m->p[m->NbRows - 1][m->NbColumns - 1]);
      if (candl_matrix_check_point (m, context))
	{
	  /* There is a point. dom(s1) - dom(s2) > 0. */
	  CANDL_clear(lb);
	  candl_matrix_free(m);
	  return 0;
	}
    }

  CANDL_clear(lb);
  candl_matrix_free(m);

  return 1;
}


/**
 * candl_dependence_extract_scalar_variables function:
 * This functions returns a -1-terminated array of the program scalar
 * variables.
 */
int*
candl_dependence_extract_scalar_variables (candl_program_p program)
{
  /* FIXME: implement a real buffer. */
  int scalars[1024];
  int checked[1024];
  int i, j, k, idx, cpt;
  int count_s = 0, count_c = 0;
  CandlMatrix* m;

  /* Detect all scalar variables. */
  for (i = 0; i < program->nb_statements; ++i)
    for (m = program->statement[i]->read, cpt = 0; cpt < 2; ++cpt,
	   m = program->statement[i]->written)
      for (j = 0; j < m->NbRows; ++j)
	{
	  idx = CANDL_get_si(m->p[j][0]);
	  if (idx != 0)
	    {
	      for (k = 0; k < count_s && scalars[k] != idx; ++k)
		;
	      if (k == count_s)
		{
		  for (k = 0; k < count_c && checked[k] != idx; ++k)
		    ;
		  if (k == count_c)
		    {
		      if (candl_dependence_var_is_scalar(program, idx))
			scalars[count_s++] = idx;
		      else
			checked[count_c++] = idx;
		    }
		  if (count_s == 1024 || count_c == 1024)
		    CANDL_FAIL("Error: Buffer size too small");
		}
	    }
	}

  /* Rearrange the array to the exact size. */
  int* res = (int*) malloc((count_s + 1) * sizeof(int));
  for (i = 0; i < count_s; ++i)
    res[i] = scalars[i];
  res[i] = -1;

  return res;
}


/**
 * candl_dependence_get_array_refs_in_dep function:
 * This function return the array indices referenced in the
 * dependence.
 */
static
void
candl_dependence_get_array_refs_in_dep (CandlDependence* tmp, int* refs,
					int* reft)
{
  if (! tmp)
    return;
  switch (tmp->type)
    {
    case CANDL_RAR:
      *refs = CANDL_get_si(tmp->source->read->p[tmp->ref_source][0]);
      *reft = CANDL_get_si(tmp->target->read->p[tmp->ref_target][0]);
      break;
    case CANDL_RAW:
      *refs = CANDL_get_si(tmp->source->written->p[tmp->ref_source][0]);
      *reft = CANDL_get_si(tmp->target->read->p[tmp->ref_target][0]);
      break;
    case CANDL_WAR:
      *refs = CANDL_get_si(tmp->source->read->p[tmp->ref_source][0]);
      *reft = CANDL_get_si(tmp->target->written->p[tmp->ref_target][0]);
      break;
    case CANDL_WAW:
      *refs = CANDL_get_si(tmp->source->written->p[tmp->ref_source][0]);
      *reft = CANDL_get_si(tmp->target->written->p[tmp->ref_target][0]);
      break;
    default:
      exit(1);
    }
}


/**
 * candl_dependence_prune_scalar_waw function:
 * This function removes all WAW dependences between the same scalar
 * (they are useless dependences).
 */
void
candl_dependence_prune_scalar_waw (candl_program_p program,
				   CandlOptions* options,
				   CandlDependence** deps)
{
  int* scalars;
  int i;
  int refs, reft;
  CandlDependence* tmp;
  CandlDependence* next;
  CandlDependence* pred = NULL;

  if (options->verbose)
    fprintf (stderr, "[Candl] Scalar Analysis: Remove all WAW between the same"
	     " scalar\n");

  scalars = candl_dependence_extract_scalar_variables (program);

  for (tmp = *deps; tmp; )
    {
      candl_dependence_get_array_refs_in_dep (tmp, &refs, &reft);
      if (refs == reft && tmp->type == CANDL_WAW)
	{
	  for (i = 0; scalars[i] != -1 && scalars[i] != refs; ++i)
	    ;
	  if (scalars[i] != -1)
	    {
	      candl_matrix_free (tmp->domain);
	      next = tmp->next;
	      if (pred == NULL)
		*deps = next;
	      else
		pred->next = next;
	      free (tmp);
	      tmp = next;
	      continue;
	    }
	}
      pred = tmp;
      tmp = tmp->next;
    }

  free(scalars);
}


/**
 * candl_dependence_scalar_renaming function:
 * This function renames scalars in the program. In case scalars have
 * been renamed, the dependence analysis is re-run.
 */
int
candl_dependence_scalar_renaming (candl_program_p program,
				  CandlOptions* options,
				  CandlDependence** deps)
{
  int i, j, k;
  CandlStatement** stmts;
  CandlStatement* defs[1024];
  CandlStatement* uses[1024];
  CandlStatement* current[program->nb_statements];
  int parts[program->nb_statements];
  CandlStatement* s;
  CandlStatement* last_def;
  CandlMatrix* m;
  int defs_c, uses_c;
  int val, cpt, tmp, has_changed = 0;
  int newvar = 0;

  for (i = 0; i < program->nb_statements; ++i)
    current[i] = NULL;

  if (options->verbose)
    fprintf (stderr, "[Candl] Scalar Analysis: Perform scalar renaming\n");

  /* Compute the first free variable index seed. */
  for (i = 0; i < program->nb_statements; ++i)
    for (m = program->statement[i]->read, cpt = 0; cpt < 2;
	 ++cpt, m = program->statement[i]->written)
      for (j = 0; j < m->NbRows; ++j)
	if (CANDL_get_si(m->p[j][0]) >= newvar)
	  newvar = CANDL_get_si(m->p[j][0]) + 1;

  /* Get the list of scalars. */
  int* scalars = candl_dependence_extract_scalar_variables (program);

  /* Iterate on all scalars. */
  for (i = 0; scalars[i] != -1; ++i)
    {
      /* Get all statements referencing the scalar. */
      stmts = candl_dependence_refvar_chain (program, NULL, scalars[i], 0);

      /* If the chain starts by a USE, we can't do anything. */
      if (stmts[0] == NULL || candl_dependence_var_is_ref (stmts[0], scalars[i])
	  != CANDL_VAR_IS_DEF)
	{
	  free (stmts);
	  continue;
	}

      /* Get all defs. */
      for (j = 0, defs_c = 0; stmts[j]; ++j)
	if (candl_dependence_var_is_ref (stmts[j], scalars[i])
	    == CANDL_VAR_IS_DEF)
	  defs[defs_c++] = stmts[j];
      /* Get all uses. */
      for (j = 0, uses_c = 0; stmts[j]; ++j)
	{
	  val = candl_dependence_var_is_ref (stmts[j], scalars[i]);
	  if (val == CANDL_VAR_IS_USED ||  val == CANDL_VAR_IS_DEF_USED)
	    uses[uses_c++] = stmts[j];
	}

      /* Clean buffer. */
      for (j = 0; j < program->nb_statements; ++j)
	current[j] = NULL;

      free (stmts);

      /* Iterate on all DEFs. */
      for (j = 0, last_def = NULL; j < defs_c; ++j)
	{
	  if (last_def == NULL)
	    last_def = defs[j];
	  else
	    {
	      /* Ensure the current DEF covers all iterations covered
		 by the last checked one. */
	      for (k = 0; k < last_def->depth && k < defs[j]->depth &&
		     last_def->index[k] == defs[j]->index[k]; ++k)
		;
	      /* We only need to check when there are common loops. */
	      if (k && ! candl_dependence_check_domain_is_included
		  (last_def, defs[j], program->context, k + 1))
		{
		  current[defs[j]->label] = last_def;
		  continue;
		}
	      else
		last_def = defs[j];
	    }
	  /* Create DEF-USE table. */
	  for (k = 0; k < uses_c; ++k)
	    if (uses[k]->label > defs[j]->label)
		current[uses[k]->label] = defs[j];
	}

      /* Initialize partition table. */
      for (j = 0; j < program->nb_statements; ++j)
	parts[j] = -1;

      /* Create partitions. */
      for (j = 0; j < defs_c; ++j)
	for (k = 0; k < program->nb_statements; ++k)
	  if ((current[k] && current[k] == defs[j]) ||
	      (k == defs[j]->label && current[defs[j]->label] == NULL))
	    parts[k] = j;

      /* Check if it is needed to rename the scalar. */
      for (j = 0, tmp = -1; j < program->nb_statements; ++j)
	if (tmp == -1)
	  {
	    if (parts[j] != -1)
	      tmp = parts[j];
	  }
	else
	  if (parts[j] != -1 && tmp != parts[j])
	    break;

      /* Rename scalar variable. */
      if (j != program->nb_statements)
	for (j = 0, tmp = -1; j < program->nb_statements; ++j)
	  if (parts[j] != -1)
	    {
	      if (tmp == -1)
		tmp = parts[j];
	      else
		if (tmp != parts[j])
		  has_changed = 1;
	      s = program->statement[j];
	      for (m = s->read, cpt = 0; cpt < 2; ++cpt, m = s->written)
		for (k = 0; k < m->NbRows; ++k)
		  if (CANDL_get_si(m->p[k][0]) == scalars[i])
		    {
		      if (options->verbose)
			fprintf (stderr, "[Candl] Scalar analysis: Renamed "
				 "variable %d to %d at statement S%d\n",
				 scalars[i], newvar + parts[j], j);
		      CANDL_set_si(m->p[k][0], newvar + parts[j]);
		    }
	    }
    }

  /* Redo the full dependence analysis, if needed. */
  if (has_changed)
    {
      int bopt = options->scalar_renaming;
      options->scalar_renaming = 0;
      if (options->scalar_privatization)
	free (program->scalars_privatizable);
      candl_dependence_free (*deps);
      *deps = candl_dependence (program, options);
      options->scalar_renaming = bopt;
    }
  free (scalars);

  return has_changed;
}


/**
 * candl_dependence_is_loop_carried function:
 * This function returns true if the dependence 'dep' is loop-carried
 * for loop 'loop_index', false otherwise.
 */
int
candl_dependence_is_loop_carried (candl_program_p program,
				  CandlDependence* dep,
				  int loop_index)
{
  int i, j;

  /* Ensure source and sink share common loop 'loop_index', and that
     dependence depth is consistent with the considered loop. */
  for (i = 0; i < dep->source->depth; ++i)
    if (dep->source->index[i] == loop_index)
      break;
  if (i != dep->depth - 1 || i >= dep->target->depth)
    return 0;
  for (j = 0; j < dep->target->depth; ++j)
    if (dep->target->index[j] == loop_index)
      break;
  if (j != i)
    return 0;

  /* Final check. The dependence exists only because the loop
     iterates. Make the loop not iterate and check if there's still
     dependent iterations. */
  CandlMatrix* m = candl_matrix_malloc(dep->domain->NbRows + 2,
				       dep->domain->NbColumns);
  CANDL_set_si(m->p[m->NbRows - 2][i + 1], -1);
  CANDL_set_si(m->p[m->NbRows - 1][dep->source->depth + 1 + j], -1);
  /* Copy the rest of the matrix. */
  for (i = 0; i < dep->domain->NbRows; ++i)
    for (j = 0; j < dep->domain->NbColumns; ++j)
      CANDL_assign(m->p[i][j], dep->domain->p[i][j]);
  /* Compute real lb of loops. */
  Entier lb; CANDL_init(lb);
  candl_dependence_compute_lb (m, &lb, i + 1);
  CANDL_assign(m->p[m->NbRows - 2][m->NbColumns - 1], lb);
  candl_dependence_compute_lb (m, &lb, dep->source->depth + 1 + j);
  CANDL_assign(m->p[m->NbRows - 1][m->NbColumns - 1], lb);
  CANDL_clear(lb);
  int ret = candl_matrix_check_point(m, program->context);

  /* Be clean. */
  candl_matrix_free(m);

  return !ret;
}


/**
 * candl_dependence_prune_with_privatization function: This function
 * prunes the dependence graph 'deps' by removing loop-carried
 * dependences involving a scalar variable privatizable for that loop.
 */
void
candl_dependence_prune_with_privatization (candl_program_p program,
					   CandlOptions* options,
					   CandlDependence** deps)
{
  CandlDependence* next;
  CandlDependence* tmp;
  CandlDependence* pred = NULL;
  int is_priv;
  int i;
  int loop_idx = 0;
  int refs, reft;


  if (options->verbose)
    fprintf (stderr, "[Candl] Scalar Analysis: Remove loop-carried dependences"
	     " on privatizable scalars\n");

  if (program->nb_statements == 0)
    return;
  /* Perform the scalar analysis, if not done before. */
  if (program->scalars_privatizable == NULL)
    {
      CandlOptions* options = candl_options_malloc ();
      options->scalar_privatization = 1;
      candl_dependence_analyze_scalars (program, options);
      candl_options_free (options);
    }

  for (tmp = *deps; tmp; )
    {
      /* Check if the dependence is involving a privatizable scalar. */
      is_priv = 1;
      candl_dependence_get_array_refs_in_dep (tmp, &refs, &reft);
      for (i = 0; i < tmp->source->depth; ++i)
	if (candl_dependence_scalar_is_privatizable_at (program, refs,
							tmp->source->index[i]))
	  break;
      if (i == tmp->source->depth)
	{
	  for (i = 0; i < tmp->target->depth; ++i)
	    if (candl_dependence_scalar_is_privatizable_at
		(program, reft, tmp->target->index[i]))
	      break;
	  if (i == tmp->target->depth)
	    is_priv = 0;
	  else
	    loop_idx = tmp->target->index[i];
	}
      else
	loop_idx = tmp->source->index[i];
      /* Check if the dependence is loop-carried at loop i. */
      if (is_priv && candl_dependence_is_loop_carried (program, tmp, loop_idx))
	{
	  /* It is, the dependence can be removed. */
	  candl_matrix_free (tmp->domain);
	  next = tmp->next;
	  if (pred == NULL)
	    *deps = next;
	  else
	    pred->next = next;
	  free (tmp);
	  tmp = next;
	  continue;
	}
      /* Go to the next victim. */
      pred = tmp;
      tmp = tmp->next;
    }
}


/**
 * candl_dependence_is_privatizable function:
 * This function checks if a given scalar 'var_index' is privatizable
 * for loop 'loop_index'.
 */
int
candl_dependence_scalar_is_privatizable_at (candl_program_p program,
					    int var_index,
					    int loop_index)
{
  int i;

  /* If the scalar analysis wasn't performed yet, do it. */
  if (program->scalars_privatizable == NULL)
    {
      CandlOptions* options = candl_options_malloc();
      options->scalar_privatization = 1;
      candl_dependence_analyze_scalars (program, options);
      candl_options_free (options);
    }

  /* Check in the array of privatizable scalar variables for the tuple
     (var,loop). */
  for (i = 0; program->scalars_privatizable[i] != -1; i += 2)
    if (program->scalars_privatizable[i] == var_index &&
	program->scalars_privatizable[i + 1] == loop_index)
      return 1;

  return 0;
}


/**
 * candl_dependence_analyze_scalars function:
 * This function checks, for all scalar variables of the program, and
 * all loop levels, if the scalar can be privatized at that level.
 */
int
candl_dependence_analyze_scalars(candl_program_p program,
				 CandlOptions* options)
{
  int* scalars;
  CandlStatement** stmts;
  CandlStatement** fullchain = NULL;
  int i, j, k, l, n;
  CandlMatrix* m;
  int idx, max, is_priv, cpt, offset, was_priv;
  CandlStatement* curlast;
  CandlStatement* last;
  int nb_priv = 0;
  int priv_buff_size = 64;

  /* Initialize the list of privatizable scalars to empty. */
  if (options->scalar_privatization)
    {
      program->scalars_privatizable = (int*) malloc(2 * sizeof(int));
      program->scalars_privatizable[0] = program->scalars_privatizable[1] = -1;
    }

  /* Retrieve all scalar variables. */
  scalars = candl_dependence_extract_scalar_variables (program);

  /* For each of those, detect (at any level) if it can be privatized
     / expanded / renamed. */
  for (i = 0; scalars[i] != -1; ++i)
    {
      idx = scalars[i];
      /* Go to the first statement referencing the scalar. */
      for (j = 0; j < program->nb_statements; ++j)
	if (candl_dependence_var_is_ref (program->statement[j], scalars[i])
	    != CANDL_VAR_UNDEF)
	  break;
      /* A weird error occured. */
      if (j == program->nb_statements)
	continue;

      /* Take all statements referencing the scalar. */
      fullchain = candl_dependence_refvar_chain (program, program->statement[j],
						 scalars[i], 0);
      /* Compute the maximum loop depth of the chain. */
      for (k = 0, max = 0; fullchain[k]; ++k)
	max = max < fullchain[k]->depth ? fullchain[k]->depth : max;
      last = fullchain[k - 1];
      /* Initialize the offset for expansion. */
      offset = 0;
      was_priv = 0;
      /* Iterate on all possible depth for analysis. */
      for (k = 1; k <= max; ++k)
	{
	  CandlStatement* s = fullchain[0];
	  if (was_priv)
	    {
	      ++offset;
	      was_priv = 0;
	    }
	  do
	    {
	      /* Take all statements dominated by s referencing the
		 current scalar variable. */
	      stmts = candl_dependence_refvar_chain (program, s, scalars[i], k);
	      /* No more statement in the chain, exit. */
	      if (stmts[0] == NULL)
		break;

	      int c = 0;
	      is_priv = candl_dependence_var_is_ref (stmts[c], scalars[i])
		== CANDL_VAR_IS_DEF;
	      /* Check for privatization, while the entry of the chain
		 is a DEF. */
	      while (stmts[c] && candl_dependence_var_is_ref
		     (stmts[c], scalars[i]) == CANDL_VAR_IS_DEF)
		{
		  /* From the current DEF node, ensure the rest of the
		     chain covers not more than the iteration domain
		     of the DEF. */
		  for (l = c + 1; stmts[l - 1] && stmts[l]; ++l)
		    /* FIXME: we should deal with
		       def_1->use->def_2->use chains where dom(def_2)
		       > dom(def_1). */
		    if (! candl_dependence_check_domain_is_included
			(stmts[c], stmts[l], program->context, k))
		      {
			/* dom(use) - dom(def) > 0. Check if there is
			   another DEF to test at the entry of the
			   block. */
			if (stmts[c + 1])
			  if (candl_dependence_var_is_ref
			      (stmts[c + 1], scalars[i]) != CANDL_VAR_IS_DEF)
			    /* No. The variable is not privatizable. */
			    is_priv = 0;
			break;
		      }
		  if (! is_priv || ! stmts[l])
		    break;
		  /* The chain dominated by stmts[c] is not
		     privatizable. Go for the next DEF at the
		     beginning of the block, if any. */
		  ++c;
		}
	      if (is_priv)
		{
		  /* Perform the privatization / expansion. */
		  if (options->verbose)
		    fprintf (stderr, "[Candl] Scalar Analysis: The variable %d"
			     " can be privatized at loop %d\n",
			     scalars[i], stmts[0]->index[k - 1]);
		  if (options->scalar_expansion)
		    /* Traverse all statements in the chain. */
		    for (l = c; stmts[l]; ++l)
		      {
			/* It's not the first expansion of the scalar,
			   we need to increase its dimension all along
			   the program. */
			if (offset && !was_priv)
			  candl_dependence_expand_scalar (fullchain,
							  scalars[i]);
			/* Perform scalar expansion in the array
			   access functions. */
			for (cpt = 0, m = stmts[l]->read; cpt < 2;
			     ++cpt, m = stmts[l]->written)
			  for (n = 0; n < m->NbRows; ++n)
			    if (CANDL_get_si(m->p[n][0]) == scalars[i])
			      CANDL_set_si(m->p[n + offset][k], 1);
			was_priv = 1;
		      }
		  if (options->scalar_privatization)
		    {
		      /* Memory management for the array of
			 privatizable scalars. */
		      if (nb_priv == 0)
			{
			  free (program->scalars_privatizable);
			  program->scalars_privatizable =
			    (int*)malloc(priv_buff_size * 2 * sizeof(int));
			  for (l = 0; l < priv_buff_size; ++l)
			    program->scalars_privatizable[l] = -1;
			}
		      if (2 * nb_priv == priv_buff_size)
			{
			  program->scalars_privatizable =
			    realloc(program->scalars_privatizable,
				    (priv_buff_size *= 2) * sizeof(int));
			  for (l = 2 * nb_priv; l < priv_buff_size; ++l)
			    program->scalars_privatizable[l] = -1;
			}
		      /* Memorize the scalar information in the
			 privatizable list. */
		      program->scalars_privatizable[nb_priv++] = scalars[i];
		      program->scalars_privatizable[nb_priv++] =
			stmts[0]->index[k - 1];
		    }
		}
	      /* Go to the next block, if any. */
	      for (l = 0; stmts[l]; ++l)
		;
	      curlast = stmts[l - 1];
	      if (curlast != last)
		{
		  for (l = 0; fullchain[l]; ++l)
		    if (fullchain[l] == curlast)
		      s = fullchain[l + 1];
		}
	      free (stmts);
	    }
	  while (curlast != last);
	}
      free (fullchain);
    }

  return 0;
}


/**
 * candl_num_dependences function:
 * This function returns the number of dependences in the dependence
 * list
 * \param dependence The first dependence of the dependence list.
 **
 */
int
candl_num_dependences(CandlDependence *candl_deps)
{
    CandlDependence *candl_dep = candl_deps;

    int num = 0;
    while (candl_dep != NULL)
    {
        num++;
        candl_dep = candl_dep->next;
    }
    return num;
}


/*
 * Convert a PIP quast to a union of polyhedra (Pip matrices)
 *
 * num: number of Pip matrices returned
 *
 **/
static
PipMatrix **quast_to_polyhedra (PipQuast *quast, int *num,
				int nvar, int npar)
{
    int num1, num2;
    PipMatrix** ep;
    PipMatrix** tp;
    PipMatrix** qp;
    int i, j;
    if (quast == NULL)
      {
        *num = 0;
        return NULL;
      }

    if (quast->condition != NULL)
      {
        tp = quast_to_polyhedra(quast->next_then, &num1, nvar, npar);
        ep = quast_to_polyhedra(quast->next_else, &num2, nvar, npar);

        /* Each of the matrices in the then tree needs to be augmented with
         * the condition */
        for (i = 0; i < num1; i++)
	  {
	    int nrows = tp[i]->NbRows;
	    CANDL_set_si(tp[i]->p[nrows][0], 1);
	    for (j = 1; j < 1 + nvar; j++)
	      CANDL_set_si(tp[i]->p[nrows][j], 0);
	    for (j = 0; j < npar + 1; j++)
	      CANDL_assign(tp[i]->p[nrows][1+nvar+j],
			   quast->condition->the_vector[j]);
	    (tp[i]->NbRows)++;
	  }

        for (i = 0; i < num2; i++)
	  {
	    int nrows = ep[i]->NbRows;
	    /* Inequality */
	    CANDL_set_si(ep[i]->p[nrows][0], 1);
	    for (j = 1; j < 1 + nvar; j++)
	      CANDL_set_si(ep[i]->p[nrows][j], 0);
	    for (j = 0; j < npar + 1; j++)
	      CANDL_set_si(ep[i]->p[nrows][1+nvar+j],
			   -CANDL_get_si(quast->condition->the_vector[j]));
	    (ep[i]->NbRows)++;
	  }

        qp = (PipMatrix **) malloc((num1 + num2) * sizeof(PipMatrix*));
	for (i = 0; i < num1; ++i)
	  qp[i] = tp[i];
	for (i = 0; i < num2; ++i)
	  qp[i + num1] = ep[i];

        *num = num1 + num2;

        return qp;

      }
    else
      {
        /* quast condition is NULL */
        PipMatrix *lwmatrix = pip_matrix_alloc(nvar+npar+1, nvar+npar+2);

        PipList *vecList = quast->list;

        int count=0;
        while (vecList != NULL) {
	  /* Equality */
	  CANDL_set_si(lwmatrix->p[count][0], 0);
	  for (j=0; j<nvar; j++)
	    if (j == count)
	      CANDL_set_si(lwmatrix->p[count][j+1], 1);
	    else
	      CANDL_set_si(lwmatrix->p[count][j+1], 0);

	  for (j=0; j<npar; j++)
	    CANDL_set_si(lwmatrix->p[count][j+1+nvar],
			 -CANDL_get_si(vecList->vector->the_vector[j]));
	  /* Constant portion */
	  if (quast->newparm != NULL)
	    /* Don't handle newparm for now */
	    CANDL_set_si(lwmatrix->p[count][npar+1+nvar],
			 -CANDL_get_si(vecList->vector->the_vector[npar+1]));
	  else
	    CANDL_set_si(lwmatrix->p[count][npar+1+nvar],
			 -CANDL_get_si(vecList->vector->the_vector[npar]));

	  count++;

	  vecList = vecList->next;
        }
        lwmatrix->NbRows = count;

        if (count > 0)
	  *num = 1;
        else
	  *num = 0;

	PipMatrix** ret = (PipMatrix**) malloc(sizeof(PipMatrix*));
	ret[0] = lwmatrix;
        return ret;
      }
}


/*
 * Compute last writer for a given dependence; does not make sense if the
 * supplied dependence is not a RAW or WAW dependence
 *
 * Returns 0 if lastwriter is computed successfully and dep domain updated,
 * returns 1 otherwise
 */
static
int candl_dep_compute_lastwriter (CandlDependence **dep, CandlProgram *prog)
{
    PipQuast *lexmax;
    PipMatrix *new_domain;

    int i, j;

    int npar = prog->context->NbColumns-2;

    PipOptions *pipOptions = pip_options_init();

    /* We do a parametric lexmax on the source iterators
     * keeping the target iterators as parameters */
    pipOptions->Maximize = 1;
    pipOptions->Simplify = 1;
    // pipOptions->Deepest_cut = 1;
    // pipOptions->Urs_unknowns = -1;
    // pipOptions->Urs_parms = -1;

    /* Build a context with equalities /inequalities only on the target
     * variables */
    PipMatrix *context = pip_matrix_alloc((*dep)->domain->NbRows,
					  (*dep)->target->depth + npar + 2);

    int nrows = 0;
    for (i = 0; i < (*dep)->domain->NbRows; i++)
      {
	for (j = 1; j < (*dep)->source->depth+1; j++)
	  {
	    if ((*dep)->domain->p[i][j] != 0)
	      break;
	  }
	if (j == (*dep)->source->depth+1)
	  {
	    /* Include this in the context */
	    CANDL_assign(context->p[nrows][0], (*dep)->domain->p[i][0]);
	    for (j = 1; j < 1 + (*dep)->target->depth + npar + 1; j++)
	      CANDL_assign(context->p[nrows][j],
			   (*dep)->domain->p[i][(*dep)->source->depth+j]);

	    nrows++;
	  }
      }
    context->NbRows = nrows;

    /* Parameteric lexmax */
    lexmax = pip_solve((*dep)->domain, context, -1, pipOptions);

    pip_options_free(pipOptions);

    if (lexmax == NULL)
      {
        printf("WARNING: last writer failed (mostly invalid dependence): %s\n",
               "bailing out safely without modification");
        pip_matrix_print(stdout, (*dep)->domain);
        pip_matrix_print(stdout, context);
        return 1;
      }

    int num;
    PipMatrix **qp = quast_to_polyhedra(lexmax, &num, (*dep)->source->depth,
					(*dep)->target->depth + npar);

    /* Update the dependence domains */
    if (num > 0)
      {
	int it_mat;
	PipMatrix* original_domain = (*dep)->domain;
	for (it_mat = 0; it_mat < num; ++it_mat)
	  {
	    new_domain = pip_matrix_alloc(original_domain->NbRows +
					  qp[it_mat]->NbRows,
					  original_domain->NbColumns);
	    for (i = 0; i < original_domain->NbRows; i++)
	      for (j = 0; j < original_domain->NbColumns; j++)
		CANDL_assign(new_domain->p[i][j], original_domain->p[i][j]);

	    for (i = 0; i < qp[it_mat]->NbRows; i++)
	      for (j = 0; j < original_domain->NbColumns; j++)
		CANDL_assign(new_domain->p[i+original_domain->NbRows][j],
			     qp[it_mat]->p[i][j]);

	    (*dep)->domain = new_domain;
	    /* More than 1 pipmatrix from the quast, we need to insert
	       new dependences to have the union of domains. */
	    if (it_mat < num - 1)
	      {
		CandlDependence* newdep = candl_dependence_malloc ();
		newdep->source = (*dep)->source;
		newdep->target = (*dep)->target;
		newdep->depth = (*dep)->depth;
		newdep->type = (*dep)->type;
		newdep->ref_source = (*dep)->ref_source;
		newdep->ref_target = (*dep)->ref_target;
		newdep->usr = (*dep)->usr;
		newdep->next = (*dep)->next;
		(*dep)->next = newdep;
		*dep = newdep;
	      }
	  }

	pip_matrix_free(original_domain);
	for (i = 0; i < num; ++i)
	  pip_matrix_free(qp[i]);
      }

    if (qp)
      free(qp);
    pip_quast_free(lexmax);
    pip_matrix_free(context);

    return 0;
}


/**
 * Compute the last writer for each RAW, WAW, and RAR dependence. This will
 * modify the dependence polyhedra. Be careful of any references to the old
 * dependence polyhedra. They are freed and new ones allocated.
 */
void candl_compute_last_writer (CandlDependence *dep, CandlProgram *prog)
{
  // int count=0;
    while (dep != NULL)    {
        if (dep->type == CANDL_RAW || dep->type == CANDL_WAW || dep->type == CANDL_RAR)   {
	  // printf("Last writer for dep %d: %d %d\n", count++, dep->source->depth, dep->target->depth);
            // candl_matrix_print(stdout, dep->domain);
            candl_dep_compute_lastwriter(&dep, prog);
            // candl_matrix_print(stdout, dep->domain);
        }
        dep = dep->next;
    }
}
