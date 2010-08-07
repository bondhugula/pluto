
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                    candl.c                              **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.        First version: september 8th 2003               **
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
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * CAnDL, the Chunky Dependence Analyser                                      *
 * Written by Cedric Bastoul, Cedric.Bastoul@inria.fr                         *
 *                                                                            *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <candl/candl.h>
#include <candl/dependence.h>
#include <candl/program.h>
#include <candl/violation.h>
#include <candl/options.h>


int main(int argc, char * argv[])
{
  CandlProgram * program;
  CandlOptions * options;
  CandlDependence * dependence;
  CandlViolation * violation = NULL;
  FILE * input, * output;

  /* Options and input/output file setting. */
  candl_options_read(argc, argv, &input, &output, &options);

  /* Reading the program informations. */
#ifdef CANDL_SUPPORTS_SCOPLIB
  if (options->readscop)
    program = candl_program_read_scop(input);
  else if (options->scoptocandl)
    {
      program = candl_program_read_scop(input);
      candl_program_print_candl_file(output, program);
      candl_program_free(program);
      candl_options_free(options);
      return 0;
    }
  else
#endif
    program = candl_program_read(input);

  /* Calculating dependences. */
  dependence = candl_dependence(program, options);

  /* Calculating legality violations. */
  if (options->violgraph)
    violation = candl_violation(program, dependence, options);

  /* Printing data structures if asked. */
  if (options->structure)
    {
      fprintf(output, "Program:\n");
      candl_program_print(output, program);
      fprintf(output, "\nDependences:\n");
      candl_dependence_print(output, dependence);
      fprintf(output, "\nViolations:\n");
      candl_violation_print(output, violation);
    }

#ifdef CANDL_SUPPORTS_SCOPLIB
  if (options->readscop && options->writescop)
    candl_dependence_print_scop (input, output, dependence);
  else
    {
#endif
      /* Printing dependence graph if asked or if there is no transformation. */
      if (options->depgraph || (program->transformation == NULL))
	{
	  candl_dependence_pprint(output, dependence);
	  if (options->view)
	    candl_dependence_view(dependence);
	}
      /* Printing violation graph if asked and if there is a transformation. */
      if (options->violgraph && (program->transformation != NULL))
	{
	  candl_violation_pprint(output, violation);
	  if (options->view)
	    candl_violation_view(violation);
	}
#ifdef CANDL_SUPPORTS_SCOPLIB
    }
#endif

  /* Being clean. */
  candl_violation_free(violation);
  candl_dependence_free(dependence);
  candl_program_free(program);
  candl_options_free(options);
  pip_close();
  fclose(input);
  fclose(output);

  return 0;
}
