
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 symbol.c                              **
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
# include <clan/symbol.h>


/*+****************************************************************************
 *                               Global variables                             *
 ******************************************************************************/


int symbol_nb_iterators  = 0; /**< Current number of iterator symbols  */
int symbol_nb_parameters = 0; /**< Current number of parameter symbols */
int symbol_nb_arrays     = 0; /**< Current number of array symbols     */
int symbol_nb_functions  = 0; /**< Current number of function symbols  */


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * clan_symbol_print_structure function:
 * Displays a clan_symbol_t structure (*symbol) into a file (file, possibly
 * stdout) in a way that trends to be understandable without falling in a deep
 * depression or, for the lucky ones, getting a headache... It includes an
 * indentation level (level) in order to work with others print_structure
 * functions.
 * \param file   File where informations are printed.
 * \param symbol The symbol whose information have to be printed.
 * \param level  Number of spaces before printing, for each line.
 **
 * - 01/05/2008: first version.
 */
void
clan_symbol_print_structure(FILE * file, clan_symbol_p symbol, int level)
{
  int i, j, first = 1, number = 1;

  if (symbol != NULL)
  {
    /* Go to the right level. */
    for(j = 0; j < level; j++)
      fprintf(file,"|\t");
    fprintf(file,"+-- clan_symbol_t (node %d)\n",number);
  }
  else
  {
    /* Go to the right level. */
    for(j = 0; j < level; j++)
      fprintf(file,"|\t");
    fprintf(file,"+-- NULL symbol\n");
  }

  while (symbol != NULL)
  { if (!first)
    {
      /* Go to the right level. */
      for (j = 0; j < level; j++)
        fprintf(file,"|\t");
      fprintf(file,"|   clan_symbol_t (node %d)\n",number);
    }
    else
      first = 0;

    /* A blank line. */
    for (j = 0; j <= level+1; j++)
      fprintf(file,"|\t");
    fprintf(file,"\n");

    /* Print the identifier. */
    for (i = 0; i <= level; i++)
      fprintf(file,"|\t");
    if (symbol->identifier != NULL)
      fprintf(file,"+-- Identifier: %s\n",symbol->identifier);
    else
      fprintf(file,"+-- No identifier\n");

    /* A blank line. */
    for(j = 0; j <= level+1; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"\n") ;

    /* Go to the right level and print the type. */
    for (j = 0; j <= level; j++)
    fprintf(file,"|\t") ;
    fprintf(file,"Type: ") ;
    switch (symbol->type)
    { case SCOPLIB_TYPE_ITERATOR : fprintf(file,"Iterator\n");  break;
      case SCOPLIB_TYPE_PARAMETER: fprintf(file,"Parameter\n"); break;
      case SCOPLIB_TYPE_ARRAY    : fprintf(file,"Array\n");     break;
      case SCOPLIB_TYPE_FUNCTION : fprintf(file,"Function\n");  break;
      default : fprintf(file,"Unknown\n") ;
    }

    /* A blank line. */
    for (j = 0; j <= level+1; j++)
      fprintf(file,"|\t");
    fprintf(file,"\n");

    /* Go to the right level and print the rank. */
    for (j = 0; j <= level; j++)
      fprintf(file,"|\t");
    fprintf(file,"Rank: %d\n",symbol->rank);

    /* A blank line. */
    for (j = 0; j <= level+1; j++)
      fprintf(file,"|\t");
    fprintf(file,"\n");

    symbol = symbol->next;
    number++;

    /* Next line. */
    if (symbol != NULL)
    {
      for (j = 0; j <= level; j++)
        fprintf(file,"|\t");
      fprintf(file,"V\n");
    }
  }

  /* The last line. */
  for(j=0; j<=level; j++)
    fprintf(file,"|\t");
  fprintf(file,"\n");
}


/**
 * clan_symbol_print function:
 * This function prints the content of a clan_symbol_t structure (*symbol) into
 * a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param symbol The symbol whose information have to be printed.
 **
 * - 01/05/2008: first version.
 */
void
clan_symbol_print(FILE * file, clan_symbol_p symbol)
{
  clan_symbol_print_structure(file,symbol,0);
}


/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/


/**
 * clan_symbol_malloc function:
 * This function allocates the memory space for a clan_symbol_t structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 **
 * - 01/05/2008: first version.
 */
clan_symbol_p
clan_symbol_malloc()
{
  clan_symbol_p symbol;

  symbol = (clan_symbol_p)malloc(sizeof(clan_symbol_t));
  if (symbol == NULL)
  {
    fprintf(stderr, "[Clan] Memory Overflow.\n");
    exit(1);
  }

  symbol->identifier = NULL;
  symbol->next       = NULL;

  return symbol;
}


/**
 * clan_symbol_free function:
 * This function frees the allocated memory for a clan_symbol_t structure.
 * \param symbol The pointer to the symbol we want to free.
 **
 * - 01/05/2008: first version.
 */
void
clan_symbol_free(clan_symbol_p symbol)
{
  clan_symbol_p next;

  while (symbol != NULL)
  {
    next = symbol->next;
    free(symbol->identifier);
    free(symbol);
    symbol = next;
  }
}


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * clan_symbol_lookup function:
 * This function searches the symbol table for a symbol with the identifier
 * provided as parameter. It returns the pointer to the symbol if it already
 * exists inside the table, NULL otherwise.
 * \param symbol     The first node of the list of symbols.
 * \param identifier The identifier we are looking for.
 **
 * - 01/05/2008: first version.
 */
clan_symbol_p
clan_symbol_lookup(clan_symbol_p symbol, char * identifier)
{
  while (symbol != NULL)
  {
    if (strcmp(symbol->identifier,identifier) == 0)
      return symbol;
    else
      symbol = symbol->next;
  }
  return NULL;
}


/**
 * clan_symbol_add function:
 * This function adds a new clan_symbol_t in the symbol table whose address
 * is provided as a parameter. If the symbol table is empty (NULL), the new
 * node will become its first element. A new node is added only if an
 * existing node with the same identifier does not already exist. It returns
 * the pointer to the symbol table node corresponding to the identifier.
 * \param location   The address of the symbol table.
 * \param identifier The identifier of the symbol we want to add.
 * \param type       The new symbol type
 * \param rank       The new symbol rank (depth if iterator, ignored otherwise)
 **
 * - 01/05/2008: first version.
 */
clan_symbol_p
clan_symbol_add(clan_symbol_p * location, char * identifier, int type, int rank)
{
  clan_symbol_p symbol;

  /* If the identifier is already in the table, do nothing. */
  symbol = clan_symbol_lookup(*location,identifier);
  if (symbol != NULL) {
    return symbol;
  }

  /* Else, we allocate and fill a new clan_symbol_t node. */
  symbol = clan_symbol_malloc();

  symbol->identifier = strdup(identifier);

  /* If the type was unknown (iterator or parameter) we know now that it is
   * a parameter, it would have been already in the table otherwise.
   */
  if (type == SCOPLIB_TYPE_UNKNOWN)
    type = SCOPLIB_TYPE_PARAMETER;
  symbol->type = type;

  switch (symbol->type)
  {
    case SCOPLIB_TYPE_ITERATOR : symbol->rank = rank;
                              symbol_nb_iterators++;
			      break;
    case SCOPLIB_TYPE_PARAMETER: symbol->rank = ++symbol_nb_parameters; break;
    case SCOPLIB_TYPE_ARRAY    : symbol->rank = ++symbol_nb_arrays;     break;
    case SCOPLIB_TYPE_FUNCTION : symbol->rank = ++symbol_nb_functions;  break;
  }

  /* We put the new symbol at the beginning of the table (easier ;-) !). */
  symbol->next = *location;
  *location = symbol;

  return symbol;
}

/**
 * clan_symbol_remove function:
 * Deletes a symbol from the list.
 *
 */
void
clan_symbol_remove(clan_symbol_p* location, clan_symbol_p symbol)
{
  clan_symbol_p s = *location;
  if (s == NULL || symbol == NULL)
    return;
  if (s == symbol)
    *location = symbol->next;
  else
    {
      while (s && s->next != symbol)
	s = s->next;
      if (s)
	{
	  s->next = symbol->next;
	  free(symbol);
	}
    }
}


/**
 * clan_symbol_get_rank function:
 * This function returns the rank of the symbol with identifier "identifier"
 * in the symbol table whose first element is "symbol". If the symbol with
 * the specified identifier is not found, it returns -1.
 * \param symbol     The first node of the list of symbols.
 * \param identifier The identifier we want to know the rank.
 **
 * - 01/05/2008: first version.
 */
int
clan_symbol_get_rank(clan_symbol_p symbol, char * identifier)
{
  while (symbol != NULL)
  {
    if (strcmp(symbol->identifier,identifier) == 0)
      return symbol->rank;
    else
      symbol = symbol->next;
  }
  return -1;
}


/**
 * clan_symbol_get_type function:
 * This function returns the type of the symbol with identifier "identifier"
 * in the symbol table whose first element is "symbol". If the symbol with
 * the specified identifier is not found, it returns -1.
 * \param symbol     The first node of the list of symbols.
 * \param identifier The identifier we want to know the type.
 **
 * - 01/05/2008: first version.
 */
int
clan_symbol_get_type(clan_symbol_p symbol, char * identifier)
{
  while (symbol != NULL)
  {
    if (strcmp(symbol->identifier,identifier) == 0)
      return symbol->type;
    else
      symbol = symbol->next;
  }
  return -1;
}


/**
 * clan_symbol_iterators function:
 * this function builds the array of original iterator names for the
 * scoplib_statement_t structures thanks to the parser current state of
 * parser_depth (depth) and parser_iterators (symbols). "symbols"
 * is an array of references to symbol table entries, one for each
 * loop enclosing the statement.
 * \param symbols Array of iterator symbols for the statement.
 * \param depth   The depth of the statement.
 **
 * - 01/05/2008: First version.
 */
char **
clan_symbol_iterators(clan_symbol_p * symbols, int depth)
{
  int i, length;
  char ** iterators;

  iterators = (char **)malloc(depth * sizeof (char *));

  for (i = 0; i < depth; i++)
  {
    length = strlen((symbols[i])->identifier) + 1;
    iterators[i] = (char *)malloc(length * sizeof(char));
    strcpy(iterators[i],(symbols[i])->identifier);
  }

  return iterators;
}


/**
 * clan_symbol_id_array function:
 * this function builds an array of identifier names of a given type
 * thanks to the informations stored in the symbol
 * table and returns it. The identifiers are sort according to their rank.
 * It also returns the array size to the parameter "size".
 * \param symbol The first element of the symbol table.
 * \param type   The type of interesting elements.
 * \param size   The returned array size (this is a _result_).
 **
 * - 02/05/2008: First version.
 * - 03/05/2008: More generic, no more dedicated to parameters only.
 */
char **
clan_symbol_id_array(clan_symbol_p symbol, int type, int * size)
{
  int i, length, nb_identifiers = 0;
  char ** identifiers;
  clan_symbol_p start;

  /* A first scan of the table to find the number of identifiers of "type". */
  start = symbol;
  while (symbol != NULL)
  {
    if (symbol->type == type)
      nb_identifiers++;
    symbol = symbol->next;
  }

  identifiers = (char **)malloc(nb_identifiers * sizeof(char *));

  /* We scan the table a second time to fill the identifier array
   * Not optimal to act this way but overkills are worse!
   */
  i = 0;
  symbol = start;
  while (symbol != NULL)
  {
    if (symbol->type == type)
    {
      length = strlen(symbol->identifier) + 1;
      identifiers[symbol->rank - 1] = (char *)malloc(length * sizeof(char));
      strcpy(identifiers[symbol->rank - 1],symbol->identifier);
      i++;
    }
    symbol = symbol->next;
  }

  if (size != NULL)
   *size = nb_identifiers;

  return identifiers;
}



