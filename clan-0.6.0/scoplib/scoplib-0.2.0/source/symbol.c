
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan/Scoplib                        **
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
 * Written by Cedric Bastoul,            Cedric.Bastoul@inria.fr              *
 *
 * Symbol support:
 *            Prasanth Chatharasi,       prasanth@iith.ac.in                  *
 *                                                                            *
 ******************************************************************************/


# include <stdlib.h>
# include <stdio.h>
# include <ctype.h>
# include <string.h>
#include <assert.h>

# include <scoplib/symbol.h>
# include <assert.h>

/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/

/**
 * scoplib_symbol_print_structure function:
 * Displays a scoplib_symbol_t structure (*symbol) into a file (file, possibly
 * stdout) in a way that trends to be understandable without falling in a deep
 * depression or, for the lucky ones, getting a headache... It includes an
 * indentation level (level) in order to work with others print_structure
 * functions.
 * \param file   File where informations are printed.
 * \param symbol The symbol whose information have to be printed.
 * \param level  Number of spaces before printing, for each line.
 */
void
scoplib_symbol_print_structure(FILE * file, scoplib_symbol_p symbol, int level)
{
  int i, j, first = 1, number = 1;

  if (symbol != NULL)
  {
    /* Go to the right level. */
    for(j = 0; j < level; j++)
      fprintf(file,"|\t");
    fprintf(file,"+-- scoplib_symbol_t (node %d)\n",number);
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
      fprintf(file,"|   scoplib_symbol_t (node %d)\n",number);
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
 * scoplib_symbol_print function:
 * This function prints the content of a scoplib_symbol_t structure (*symbol) into
 * a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param symbol The symbol whose information have to be printed.
 */
void
scoplib_symbol_print(FILE * file, scoplib_symbol_p symbol)
{
  scoplib_symbol_print_structure(file,symbol,0);
}

/**
 * scoplib_symbol_print_dot_scop function:
 * This function prints the content of a scoplib_symbol_t structure (*symbol) into
 * a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param symbol The symbol whose information have to be printed.
 **
 */
void 
scoplib_symbol_print_dot_scop(FILE* file,scoplib_symbol_p symbol,int nb_parameters,
                              char** parameters,int nb_arrays,char** arrays)
{
  scoplib_symbol_p tempsymbol = symbol;
  int num_of_symbols=0;
  while(tempsymbol!=NULL) {
    num_of_symbols++;
    tempsymbol=tempsymbol->next;
  }
  fprintf(file,"# Number of symbols in the symbol list \n");
  fprintf(file,"%d\n\n",num_of_symbols);
  
  tempsymbol = symbol;
  num_of_symbols = 1;
  
  while(tempsymbol!= NULL) {
  
  int num_of_attributes=1;  
  fprintf(file,"<symbol>\n\n");
  fprintf(file,"# ............................................");
  fprintf(file,"%d.%d Symbol Data Type\n",num_of_symbols,num_of_attributes);
  num_of_attributes++;
  
  if(tempsymbol->data_type != NULL) {
    fprintf(file,"# symbol datatype is  provided\n");
    fprintf(file,"%d\n",1);    
    fprintf(file,"# symbol datatype \n");
    fprintf(file,"%s\n\n",tempsymbol->data_type);    
  }
  else {
    fprintf(file,"# symbol datatype is not provided\n");
    fprintf(file,"%d\n\n",0);
  }
  
  fprintf(file,"# ............................................");
  fprintf(file,"%d.%d Symbol Name\n",num_of_symbols,num_of_attributes);
  num_of_attributes++;
  
  if(tempsymbol->identifier != NULL) {
    fprintf(file,"# symbol name is  provided\n");
    fprintf(file,"%d\n",1);    
    fprintf(file,"# symbol name \n");
    fprintf(file,"%s\n\n",tempsymbol->identifier);    
  }
  else {
    fprintf(file,"# symbol name is not provided\n");
    fprintf(file,"%d\n\n",0);
  }  
  
  fprintf(file,"# ............................................");
  fprintf(file,"%d.%d Symbol Type\n",num_of_symbols,num_of_attributes);
  num_of_attributes++;
  
  fprintf(file,"# symbol type \n");
  switch (tempsymbol->type)
  { case SCOPLIB_TYPE_ITERATOR : fprintf(file,"Iterator\n\n");  break;
    case SCOPLIB_TYPE_PARAMETER: fprintf(file,"Parameter\n\n"); break;
    case SCOPLIB_TYPE_ARRAY    : fprintf(file,"Array\n\n");     break;
    case SCOPLIB_TYPE_FUNCTION : fprintf(file,"Function\n\n");  break;
    default : fprintf(file,"Unknown\n") ;
  }
  
  fprintf(file,"# ............................................");
  fprintf(file,"%d.%d Symbol Flag\n",num_of_symbols,num_of_attributes);
  num_of_attributes++;  
  
  fprintf(file,"# Flag for the symbol\n");
  fprintf(file,"%d\n\n",tempsymbol->flag);
  
  fprintf(file,"# ............................................");
  fprintf(file,"%d.%d Symbol Dimension\n",num_of_symbols,num_of_attributes);
  num_of_attributes++;  
  
  if(tempsymbol->dimensions_bounds != NULL) {
    fprintf(file,"# Dimension Information is  provided\n");
    fprintf(file,"%d\n",1);    
    
          
    scoplib_matrix_print_dot_scop(file,tempsymbol->dimensions_bounds,SCOPLIB_TYPE_SYMBOL_TABLE,
                                  -1,NULL,nb_parameters,parameters,0,NULL); 
                                  
/*    int i =0;
    for(i=0;i< tempsymbol->dimensions_bounds->NbRows;i++) {
      printf("Bound:%s\n",scoplib_symbol_table_get_bound(tempsymbol, i+1, parameters, nb_parameters));
    } */
    fprintf(file,"\n");
  }
  else {
    fprintf(file,"# Dimension Information is not provided\n");
    fprintf(file,"%d\n\n",0);
  }   
  
  
  fprintf(file,"</symbol>\n\n\n");
  num_of_symbols++;
  tempsymbol = tempsymbol->next;
  }
  return ;
}

/*+****************************************************************************
 *                    Memory allocation/deallocation function                 *
 ******************************************************************************/


/**
 * scoplib_symbol_malloc function:
 * This function allocates the memory space for a scoplib_symbol_t structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 */
scoplib_symbol_p
scoplib_symbol_malloc()
{
  scoplib_symbol_p symbol;

  symbol = (scoplib_symbol_p)malloc(sizeof(scoplib_symbol_t));
  if (symbol == NULL)
  {
    fprintf(stderr, "[Clan] Memory Overflow.\n");
    exit(1);
  }

  symbol->identifier        = NULL;
  symbol->next              = NULL;
  symbol->dimensions_bounds = NULL;
  symbol->data_type         = NULL;
  symbol->num_of_dimensions = 0;  
  symbol->flag              = 0;  

  return symbol;
}


/**
 * scoplib_symbol_free function:
 * This function frees the allocated memory for a scoplib_symbol_t structure.
 * \param symbol The pointer to the symbol we want to free.
 */
void
scoplib_symbol_free(scoplib_symbol_p symbol)
{
  scoplib_symbol_p next;

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
 * scoplib_symbol_lookup function:
 * This function searches the symbol table for a symbol with the identifier
 * provided as parameter. It returns the pointer to the symbol if it already
 * exists inside the table, NULL otherwise.
 * \param symbol     The first node of the list of symbols.
 * \param identifier The identifier we are looking for.
 */
scoplib_symbol_p
scoplib_symbol_lookup(scoplib_symbol_p symbol, char * identifier)
{
  scoplib_symbol_p tempsymbol = symbol;
  while (tempsymbol != NULL)
  {
    if (strcmp(tempsymbol->identifier,identifier) == 0) {
      return tempsymbol;
    }
    else
      tempsymbol = tempsymbol->next;
  }
  return NULL;
}


/**
 * scoplib_symbol_add function:
 * This function adds a new scoplib_symbol_t in the symbol table whose address
 * is provided as a parameter. If the symbol table is empty (NULL), the new
 * node will become its first element. A new node is added only if an
 * existing node with the same identifier does not already exist. It returns
 * the pointer to the symbol table node corresponding to the identifier.
 * \param location               The address of the symbol table.
 * \param identifier             The identifier of the symbol we want to add.
 * \param type                   The new symbol type
 * \param dimensions_bounds      The constraints on dimensions
 * \param num_of_dimensions      Total number of dimensions
 */
scoplib_symbol_p
scoplib_symbol_add(scoplib_symbol_p * location, char * identifier, int type ,
                   scoplib_matrix_p dimensions_bounds,int num_of_dimensions)
{
  scoplib_symbol_p symbol;

  /* If the identifier is already in the table, do nothing. */
  symbol = scoplib_symbol_lookup(*location,identifier);
  if (symbol != NULL) { 
    return symbol;
  }

  /* Else, we allocate and fill a new scoplib_symbol_t node. */
  symbol = scoplib_symbol_malloc();

  symbol->identifier = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
  strcpy(symbol->identifier,identifier);

  
  /* If the type was unknown (iterator or parameter) we know now that it is
   * a parameter, it would have been already in the table otherwise.
   */
  if (type == SCOPLIB_TYPE_UNKNOWN)
    type = SCOPLIB_TYPE_PARAMETER;
  symbol->type = type;
  
    
  /*Adding the dimension bounds to the symbol entry */
  
  symbol->dimensions_bounds = dimensions_bounds;   
  symbol->num_of_dimensions = num_of_dimensions;   
  int i = 0;
  for(i = 0; i<num_of_dimensions;i++) {
    dimensions_bounds->p[i][0] = 1;  
  }

  /* We put the new symbol at the beginning of the table (easier ;-) !). */   
  
  symbol->next = *location;
  *location = symbol;
 
  return symbol;
}

scoplib_symbol_p
scoplib_symbol_add_with_structure(scoplib_symbol_p * symbolHead,scoplib_symbol_p symbol) {

  scoplib_symbol_p newsymbol = scoplib_symbol_malloc();
  newsymbol->identifier  = strdup(symbol->identifier);
  newsymbol->data_type   = strdup(symbol->data_type);
  newsymbol->type        = symbol->type;
  newsymbol->flag        = symbol->flag;
  newsymbol->num_of_dimensions = symbol->num_of_dimensions;
  newsymbol->dimensions_bounds = scoplib_matrix_copy(symbol->dimensions_bounds);
  
  newsymbol->next = *symbolHead;
  *symbolHead = newsymbol;
  
  return newsymbol;
  
}

/**
 * scoplib_symbol_remove function:
 * Deletes a symbol from the list.
 *
 */
void
scoplib_symbol_remove(scoplib_symbol_p* location, scoplib_symbol_p symbol)
{
  scoplib_symbol_p s = *location;
  if (s == NULL || symbol == NULL)
    return;
  if (s == symbol)
    *location = symbol->next;
  else {
    while (s && s->next != symbol)
	    s = s->next;
    if (s) {
      s->next = symbol->next;
      free(symbol);
    }
  }
}


/**
 * scoplib_symbol_get_type function:
 * This function returns the type of the symbol with identifier "identifier"
 * in the symbol table whose first element is "symbol". If the symbol with
 * the specified identifier is not found, it returns -1.
 * \param symbol     The first node of the list of symbols.
 * \param identifier The identifier we want to know the type.
 */
int
scoplib_symbol_get_type(scoplib_symbol_p symbol, char * identifier)
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
 * scoplib_symbol_table_compact function:
 * This function scans the symbols list to put the right number of columns
 * to every matrix (during construction we used CLAN_MAX_DEPTH and
 * CLAN_MAX_PARAMETERS to define matrix and vector sizes).
 * \param statement     The first statement to scan to compact matrices.
 * \param nb_parameters The true number of parameters in the SCoP.
 * \param max_depth     The maximum depth
 */
void
scoplib_symbol_table_compact(scoplib_symbol_p symbol, int nb_parameters,int max_depth)
{
  int nb_iterators;
  scoplib_symbol_p tempsymbol = symbol;

  while (tempsymbol != NULL) {
    nb_iterators = tempsymbol->num_of_dimensions;
    if(tempsymbol->dimensions_bounds != NULL ) {
      scoplib_matrix_compact(tempsymbol->dimensions_bounds,nb_iterators,nb_parameters,max_depth);
      int num_of_dimensions = tempsymbol->num_of_dimensions;
      while(num_of_dimensions >= 0) {          
        tempsymbol->dimensions_bounds = scoplib_matrix_remove_column(tempsymbol->dimensions_bounds,0);
        num_of_dimensions--;          
      }
    }
    else {
      if(tempsymbol->type == SCOPLIB_TYPE_ARRAY) {
        tempsymbol->num_of_dimensions = 1;
        tempsymbol->dimensions_bounds = scoplib_matrix_malloc(1,nb_parameters+1);
      }     
    }
    tempsymbol = tempsymbol->next;
  }
} 

/**
 * scoplib_symbol_copy function:
 * This function essentially copies the symbol and returns the copy of symbol
 */

scoplib_symbol_p
scoplib_symbol_copy(scoplib_symbol_p symbol) 
{
  scoplib_symbol_p newsymbol = scoplib_symbol_malloc();
  newsymbol->identifier  = strdup(symbol->identifier);
  if(symbol->data_type  != NULL)
    newsymbol->data_type = strdup(symbol->data_type);
  newsymbol->type        = symbol->type;
  newsymbol->flag        = symbol->flag;
  newsymbol->num_of_dimensions = symbol->num_of_dimensions;
  if(symbol->dimensions_bounds != NULL)
    newsymbol->dimensions_bounds = scoplib_matrix_copy(symbol->dimensions_bounds);  
  return newsymbol;
}

/**
  * scoplib_symbol_table_get_bound function:
  * This function will give bound of given dimension in string format 
  * \param symbol Symbol 
  * \param dimension starts from 1 to num of dimensions
  * \param parameter names array
  * \param number of parameters
  */

char*
scoplib_symbol_table_get_bound(scoplib_symbol_p symbol , int dimension, 
                               char** params, int nb_params) 
{
    assert(symbol->dimensions_bounds != NULL);

    int num_dimensions = symbol->dimensions_bounds->NbRows ;

    assert(dimension <= num_dimensions-1);

    int i=0;
    scoplib_int_t value;
    int set_signs = 0;
    int one_set_signs = 0;
    char* bound = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
    char* temp_bound = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
    bound[0] = '\0';
    temp_bound[0] = '\0';

    for(i=0; i< symbol->dimensions_bounds->NbColumns; i++ ) {
        value = symbol->dimensions_bounds->p[dimension][i] ;

        if(SCOPVAL_notzero_p(value)) {

            if(set_signs == 1 ){
                if(SCOPVAL_pos_p(value))
                    strcat(bound,strdup("+"));        
            }

            if(SCOPVAL_one_p(value)) { one_set_signs = 1 ;}
            else if(SCOPVAL_mone_p(value)) { one_set_signs = 1 ;}
            else {
                SCOPVAL_sprint(temp_bound,SCOPLIB_FMT_TXT,value);
                strcat(bound,temp_bound);
            }

            if( i != symbol->dimensions_bounds->NbColumns -1) {
                if(one_set_signs == 1)
                    sprintf(temp_bound,"%s",params[i]);        
                else
                    sprintf(temp_bound,"*%s",params[i]);        
                strcat(bound,temp_bound);        
            }
            set_signs = 1;
            one_set_signs = 0;
        }    
    }
    return bound;
}
