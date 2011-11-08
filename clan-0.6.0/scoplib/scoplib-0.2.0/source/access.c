
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 access.c                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 14/06/2011                     **
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
 * Written by Cedric Bastoul,         Cedric.Bastoul@inria.fr                 *
 *
 * Access support
 *            Prasanth Chatharasi     prasanth@iith.ac.in                     *
 ******************************************************************************/
 
# include <scoplib/access.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>



/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/
 
/**
 * scoplib_access_print function:
 * This function prints the content of a scoplib_access_t structure
 * into a file (file, possibly stdout).
 * \param file   File where informations are printed.
 * \param access The access whose information have to be printed.
 *
 */ 

void
scoplib_access_print(FILE * file, scoplib_access_p access) {

  scoplib_symbol_print(file,access->symbol);
  scoplib_matrix_print(file,access->matrix);
  
}


/**
 * scoplib_access_list_print function:
 * Displays a scoplib_matrix_list_t structure (a list of matrices) into a
 * file (file, possibly stdout). See scoplib_matrix_print_structure for
 * more details.
 * \param file   File where informations are printed.
 * \param list	 The list of accesses whose information have to be printed.
 */ 
 
 
void
scoplib_access_list_print(FILE * file,scoplib_access_list_p list) {

  scoplib_access_list_p templist = list;
  while(templist!= NULL) {
    scoplib_access_print(file,templist->elt);
    templist = templist->next;
  }
} 
 
 
/**
 * scoplib_access_print_dot_scop function:
 * This function prints the content of a scoplib_access_t structure
 * into a file (file, possibly stdout) in the dot scop format
 * \param file   File where informations are printed.
 * \param access The access whose information have to be printed.
 *
 */ 

void
scoplib_access_print_dot_scop(FILE * file,scoplib_scop_p scop, scoplib_access_p access) {

  scoplib_symbol_print_dot_scop(file,access->symbol,scop->nb_parameters,scop->parameters,
                                scop->nb_arrays,scop->arrays);
  scoplib_matrix_print(file,access->matrix);
  
}


/**
 * scoplib_access_list_print_dot_scop function:
 * Displays a scoplib_matrix_list_t structure (a list of matrices) into a
 * file (file, possibly stdout). See scoplib_matrix_print_structure for
 * more details.
 * \param file   File where informations are printed.
 * \param list	 The list of accesses whose information have to be printed.
 */ 
 
 
void
scoplib_access_list_print_dot_scop(FILE * file,scoplib_scop_p scop,scoplib_access_list_p list) {

  scoplib_access_list_p templist = list;
  while(templist!= NULL) {
    scoplib_access_print_dot_scop(file,scop,templist->elt);
    templist = templist->next;
  }
} 
 
/*+****************************************************************************
 *                   Memory allocation/deallocation functions                 *
 ******************************************************************************/
 
 /**
 * scoplib_access_malloc function:
 * This function allocates the memory space for a scoplib_access_t structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * \param NbRows    The number of row of the matrix to allocate.
 * \param NbColumns The number of columns of the matrix to allocate.
 */
 
scoplib_access_p
scoplib_access_malloc(unsigned NbRows, unsigned NbColumns)
{
  scoplib_access_p access = (scoplib_access_p) malloc(sizeof(scoplib_access_t));
  access->matrix = scoplib_matrix_malloc(NbRows,NbColumns);
  access->symbol = scoplib_symbol_malloc();
  return access;
}  

/**
 * scoplib_access_free function:
 * \param access The pointer to the access we want to free.
 *
 */
void
scoplib_access_free(scoplib_access_p access)
{
  if (access == NULL)
    return;

  scoplib_symbol_free(access->symbol);
  scoplib_matrix_free(access->matrix);
  free(access);
}

/**
 * scoplib_access_list_malloc function:
 * This function allocates the memory space for a scoplib_access_list_t
 * structure and sets its fields with default values. Then it returns
 * a pointer to the allocated space.
 */
scoplib_access_list_p
scoplib_access_list_malloc()
{
  scoplib_access_list_p res =
    (scoplib_access_list_p ) malloc(sizeof(scoplib_access_list_t));

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
 * scoplib_access_list_free function:
 * This function frees the allocated memory for a scoplib_access_list_t
 * structure, and all the matrices stored in the list.
 * \param list The pointer to the access list we want to free.
 */
void
scoplib_access_list_free(scoplib_access_list_p list)
{
  scoplib_access_list_p tmp;

  if (list == NULL)
    return;

  while (list) {
      if (list->elt)
      	scoplib_access_free(list->elt);
      tmp = list->next;
      free(list);
      list = tmp;
  }
}

/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
/**
 * scoplib_access_list_add_access function:
 * This function essentially takes head of the list and the node value that is 
 * access that is going to be added to the list and appends to the list at last.
 * \param list     The address of the head of the list
 * \param access   The new node value
 */ 
 
void
scoplib_access_list_add_access(scoplib_access_list_p* list,scoplib_access_p access) 
{
  if((*list)->elt != NULL ) {
      scoplib_access_list_p newlist = scoplib_access_list_malloc();
      newlist->elt = access;
      newlist->next = NULL;
      scoplib_access_list_p temp_list = (*list);
      while(temp_list->next != NULL)
          temp_list = temp_list->next;
      temp_list->next = newlist;
  }
  else
      (*list)->elt = access;

}



/**
 * scoplib_access_matrix_to_access_format function:
 * This function essentially takes read/ write access matrices and convert them
 * to scoplib_access_list_p objects for easily reading the data into manipulations.
 * \param matrix   The matrix to be transformed 
 * \param arrays   The names of the arrays
 * \param nb_arrys Total number of the arrays
 */

    scoplib_access_list_p
scoplib_access_matrix_to_access_format(scoplib_scop_p scop,
        scoplib_matrix_p original_matrix) 
{

    scoplib_access_list_p access_list = scoplib_access_list_malloc();
    scoplib_access_p access  = (scoplib_access_p) malloc(sizeof(scoplib_access_t));
    scoplib_symbol_p symbol  = NULL;
    scoplib_matrix_p matrix;
    int i;
    int break_point = 0;
    int start_point = 0;

    for(i=0;i<original_matrix->NbRows;i++) {
        if(original_matrix->p[i][0] != 0) {
            /* Adding to the access list */
            if(break_point == 1) {
                int j = start_point;
                access  = (scoplib_access_p) malloc(sizeof(scoplib_access_t));        
                matrix = scoplib_matrix_malloc(i-start_point,original_matrix->NbColumns);
                for(;j<i;j++) {
                    scoplib_vector_p tempvector = scoplib_matrix_get_row(original_matrix,j);
                    scoplib_matrix_add_vector(matrix,tempvector,j-start_point);
                }      
                access->matrix = scoplib_matrix_remove_column(scoplib_matrix_copy(matrix),0);
                access->symbol = symbol;
                scoplib_access_list_add_access(&access_list,access);
                start_point = i;      
                free(matrix);
            }

            char* array_name = strdup((scop->arrays)[(original_matrix->p[i][0])-1]);
            scoplib_symbol_p temp_symbol = scoplib_symbol_lookup(scop->symbol_table,array_name);
            if (temp_symbol != NULL) {
                symbol = scoplib_symbol_copy(temp_symbol);
            }
        }

        /* Setting up the break point */
        if(i+1 < original_matrix->NbRows ) {
            if(original_matrix->p[i+1][0] != 0) {
                break_point = 1;
            }
            else
                break_point = 0;      
        }
        /* Adding the last array accesses in the matrix */
        else {
            i++;
            int j = start_point;
            access  = (scoplib_access_p) malloc(sizeof(scoplib_access_t));       
            matrix = scoplib_matrix_malloc(i-start_point,original_matrix->NbColumns);
            for(;j<i;j++) {
                scoplib_vector_p tempvector = scoplib_matrix_get_row(original_matrix,j);
                scoplib_matrix_add_vector(matrix,tempvector,j-start_point);
            }      
            access->matrix = scoplib_matrix_remove_column(scoplib_matrix_copy(matrix),0);
            access->symbol = symbol;
            scoplib_access_list_add_access(&access_list,access);
            start_point = i;      
            free(matrix);
        }          
    } 

    if(access_list->elt != NULL)
        return access_list;
    else
        return NULL;
}

/**
 * scoplib_access_get_read_access_list function:
 * This function essentially takes read access matrices and convert them
 * to scoplib_access_list_p objects for easily reading the data into manipulations.
 * and then returning the read access list
 * \param scop      The scop structure
 * \param statement The statement in which read access matrix has to be transformed
 */

scoplib_access_list_p
scoplib_access_get_read_access_list(scoplib_scop_p scop, scoplib_statement_p statement) {
    return scoplib_access_matrix_to_access_format(scop,statement->read);
}


/**
 * scoplib_access_get_write_access_list function:
 * This function essentially takes write access matrices and convert them
 * to scoplib_access_list_p objects for easily reading the data into manipulations.
 * and then returning the read access list
 * \param scop      The scop structure
 * \param statement The statement in which write access matrix has to be transformed
 */

scoplib_access_list_p
scoplib_access_get_write_access_list(scoplib_scop_p scop ,scoplib_statement_p statement) {
    return scoplib_access_matrix_to_access_format(scop,statement->write);
}
