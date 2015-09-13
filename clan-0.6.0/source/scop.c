
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                  scop.c                               **
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
# include <ctype.h>
# include <string.h>
# include <clan/scop.h>
# include <clan/clan.h>
# include <parser.h>


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/

scoplib_scop_p
clan_parse(FILE *, clan_options_p);

/**
 * clan_scop_extract function:
 * this function is a wrapper clan_parse function that parses a file to
 * extract a SCoP and returns, if successful, a pointer to the scoplib_scop_t
 * structure.
 * \param input   The file to parse (already open).
 * \param options Options for file parsing.
 **
 * - 24/05/2008: First version.
 */
scoplib_scop_p
clan_scop_extract(FILE * input, clan_options_p options)
{
  return clan_parse(input,options);
}


/**
 * clan_scop_compact function:
 * This function scans the SCoP to put the right number of columns to every
 * matrix (during construction we used CLAN_MAX_DEPTH and CLAN_MAX_PARAMETERS
 * to define matrix and vector sizes).
 * \param scop The scop to scan to compact matrices.
 **
 * - 02/05/2008: first version.
 */
void
clan_scop_compact(scoplib_scop_p scop,int clan_max_depth)
{
  clan_statement_compact(scop->statement,scop->nb_parameters);
  scoplib_symbol_table_compact(scop->symbol_table,scop->nb_parameters,clan_max_depth);
}


/**
 * clan_scop_fill_options:
 * This function stores the list of variables id in 'varlist' in the
 * option tag of the scop, enclosed by
 * <local-vars></local-vars>.
 *
 */
void
clan_scop_fill_options(scoplib_scop_p scop, int* localvars, int* liveoutvars)
{
  /* Build the string of ids. */
  int i, size;
  char* tag = NULL;
  char* tag1 = NULL;
  char* tag2 = NULL;
  if (localvars && localvars[0] != -1)
    {
      /* localvars is a -1-terminated array. */
      for (i = 0; localvars[i] != -1; ++i)
	;
      size = i;
      char* ids = (char*)malloc(((size * 5) + 1) * sizeof(char));
      ids[0] = '\0';
      char buffer[16];
      for (i = 0; i < size; ++i)
	{
	  if (i == 0)
	    sprintf(buffer, "%d", localvars[i]);
	  else
	    sprintf(buffer, " %d", localvars[i]);
	  strcat(ids, buffer);
	}
      size = strlen("<local-vars>\n") + strlen (ids) +
	strlen ("</local-vars>\n");
      tag1 = (char*)malloc((size + 1) * sizeof(char));
      strcpy(tag1, "<local-vars>\n");
      strcat(tag1, ids);
      strcat(tag1, "\n");
      strcat(tag1, "</local-vars>\n");
      free(ids);
    }

  if (liveoutvars && liveoutvars[0] != -1)
    {
      /* liveoutvars is a -1-terminated array. */
      for (i = 0; liveoutvars[i] != -1; ++i)
	;
      size = i;
      char* ids = (char*)malloc(((size * 5) + 1) * sizeof(char));
      ids[0] = '\0';
      char buffer[16];
      for (i = 0; i < size; ++i)
	{
	  if (i == 0)
	    sprintf(buffer, "%d", liveoutvars[i]);
	  else
	    sprintf(buffer, " %d", liveoutvars[i]);
	  strcat(ids, buffer);
	}
      size = strlen("<live-out-vars>\n") + strlen (ids) +
	strlen ("</live-out-vars>\n");
      tag2 = (char*)malloc((size + 1) * sizeof(char));
      strcpy(tag2, "<live-out-vars>\n");
      strcat(tag2, ids);
      strcat(tag2, "\n");
      strcat(tag2, "</live-out-vars>\n");
      free(ids);
    }

  /* Concatenate the tags. */
  if (tag1 || tag2)
    {
      if (tag1 == NULL)
	tag = tag2;
      else if (tag2 == NULL)
	tag = tag1;
      else
	{
	  tag = (char*)malloc((strlen(tag1) + strlen(tag2) + 1));
	  strcpy(tag, tag1);
	  strcat(tag, tag2);
	  free(tag1);
	  free(tag2);
	}
    }

  if (scop->optiontags == NULL)
    scop->optiontags = tag;
  else
    {
      char* newtag = (char*)malloc((strlen(tag) + strlen(scop->optiontags) + 2)
				   * sizeof(char));
      strcpy(newtag, scop->optiontags);
      strcat(newtag, "\n");
      strcat(newtag, tag);
      free(scop->optiontags);
      scop->optiontags = newtag;
    }
}
