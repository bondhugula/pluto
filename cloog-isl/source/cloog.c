
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                             cloog.c                               **
    **-------------------------------------------------------------------**
    **       First version: october 25th 2001, CLooG's birth date !      **
    **-------------------------------------------------------------------**/


/******************************************************************************
 *               CLooG : the Chunky Loop Generator (experimental)             *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2001-2005 Cedric Bastoul                                     *
 *                                                                            *
 * This library is free software; you can redistribute it and/or              *
 * modify it under the terms of the GNU Lesser General Public                 *
 * License as published by the Free Software Foundation; either               *
 * version 2.1 of the License, or (at your option) any later version.         *
 *                                                                            *
 * This library is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU          *
 * Lesser General Public License for more details.                            *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public           *
 * License along with this library; if not, write to the Free Software        *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,                         *
 * Boston, MA  02110-1301  USA                                                *
 *                                                                            *
 * CLooG, the Chunky Loop Generator                                           *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


# include <stdlib.h>
# include <stdio.h>
# include "../include/cloog/cloog.h"


int main(int argv, char * argc[])
{ CloogProgram * program ;
  CloogOptions * options ;
  CloogState *state;
  FILE * input, * output ;
   
  state = cloog_state_malloc();

  /* Options and input/output file setting. */
  cloog_options_read(state, argv, argc, &input, &output, &options);

  /* Reading the program informations. */
  program = cloog_program_read(input,options) ;
  fclose(input) ;
  
  /* Generating and printing the code. */
  program = cloog_program_generate(program,options) ;
  if (options->structure)
  cloog_program_print(stdout,program) ;
  cloog_program_pprint(output,program,options) ;
  cloog_program_free(program) ;

  /* Printing the allocation statistics if asked. */
  if (options->leaks) {
    fprintf(output,"/* Domains    : allocated=%5d, freed=%5d, max=%5d. */\n",
           state->domain_allocated, state->domain_freed, state->domain_max);
    fprintf(output,"/* Loops      : allocated=%5d, freed=%5d, max=%5d. */\n",
           state->loop_allocated, state->loop_freed, state->loop_max);
    fprintf(output,"/* Statements : allocated=%5d, freed=%5d, max=%5d. */\n",
           state->statement_allocated, state->statement_freed, state->statement_max);
    fprintf(output,"/* Blocks     : allocated=%5d, freed=%5d, max=%5d. */\n",
           state->block_allocated, state->block_freed, state->block_max);
  }

  /* Inform the user in case of a problem with the allocation statistics. */
  if ((state->domain_allocated    != state->domain_freed)    ||
      (state->loop_allocated      != state->loop_freed)      ||
      (state->statement_allocated != state->statement_freed) ||
      (state->block_allocated     != state->block_freed))
  {
    cloog_msg(options, CLOOG_INFO,
            "an internal problem has been detected (it should have"
	    " no\n             consequence on the correctness of the output)."
	    " Please send (if\n	     you can) your input file, the first line "
	    "given by typing 'cloog -v'\n	     and your full command "
            "line call to CLooG including options to\n	     <cedric.bastoul"
	    "@inria.fr>. Thank you for your participation to get\n"
	    "	     CLooG better and safer.\n") ;
  }

  cloog_options_free(options) ;
  cloog_state_free(state);
  fclose(output) ;
  return 0;
}

