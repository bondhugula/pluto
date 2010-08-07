
   /**------ ( ----------------------------------------------------------**
    **       )\                      CAnDL                               **
    **----- /  ) --------------------------------------------------------**
    **     ( * (                   options.c                             **
    **----  \#/  --------------------------------------------------------**
    **    .-"#'-.         First version: september 8th 2003              **
    **--- |"-.-"| -------------------------------------------------------**
          |     |
          |     |
 ******** |     | *************************************************************
 * CAnDL  '-._,-' the Chunky Analyzer for Dependences in Loops (experimental) *
 ******************************************************************************
 *                                                                            *
 * Copyright (C) 2003 Cedric Bastoul                                          *
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
#include <string.h>
#include <candl/candl.h>
#include <candl/options.h>


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * candl_option_print function:
 * This function prints the content of a CandlOptions structure (program) into
 * a file (foo, possibly stdout).
 * April 19th 2003: first version.
 */
void candl_options_print(FILE * foo, CandlOptions * options)
{
  fprintf(foo, "Options:\n");
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * candl_options_free function:
 * This function frees the allocated memory for a CandlOptions structure.
 * April 19th 2003: first version.
 */
void candl_options_free(CandlOptions * options)
{
  free(options);
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * candl_options_malloc function:
 * This functions allocate the memory space for a CandlOptions structure and
 * fill its fields with the defaults values. It returns a pointer to the
 * allocated CandlOptions structure.
 * April 19th 2003: first version.
 */
CandlOptions * candl_options_malloc(void)
{
  CandlOptions * options;

  /* Memory allocation for the CandlOptions structure. */
  options = (CandlOptions *) malloc(sizeof(CandlOptions));
  if (options == NULL)
  {
    fprintf(stderr, "[Candl]ERROR: memory overflow.\n");
    exit(1);
  }

  /* We set the various fields with default values. */
  /* OPTIONS FOR DEPENDENCE COMPUTATION */
  options->waw = 1;       /* WAW (output) dependences matter. */
  options->raw = 1;       /* RAW (flow)   dependences matter. */
  options->war = 1;       /* WAR (anti)   dependences matter. */
  options->rar = 0;       /* RAR (input)  dependences don't matter. */
  options->commute = 0;   /* Don't use commutativity to simplify dependences.*/
  options->fullcheck = 0; /* Don't compute all violations.*/
  options->depgraph  = 0; /* Don't print the dependence graph.*/
  options->violgraph = 0; /* Don' compute the violation graph.*/
  options->scalar_renaming = 0; /* Don't enable scalar renaming. */
  options->scalar_privatization = 0; /* Don't enable scalar privatization. */
  options->scalar_expansion = 0; /* Don't enable scalar expansion. */
  options->lastwriter = 0; /* Compute the last writer for RAW and WAW dependences */
  options->readscop = 0; /* Don't read a .scop format for the input. */
  options->writescop = 0; /* Don't write a .scop format for the output. */
  options->scoptocandl = 0; /* Don't act as a .scop to candl converter. */
  options->verbose = 0; /* Don't be verbose. */
  /* UNDOCUMENTED OPTIONS FOR THE AUTHOR ONLY */
  options->view = 0;      /* Do not visualize the graph with dot and gv.*/
  options->structure = 0; /* Don't print internal dependence structure. */

  return options;
}


/**
 * candl_options_help function:
 * This function displays the quick help when the user set the option -help
 * while calling candl. Prints are cutted to respect the 509 characters
 * limitation of the ISO C 89 compilers.
 * August 5th 2002: first version.
 */
void candl_options_help()
{ printf(
  "Usage: candl [ options | file ] ...\n"
  "Options for data dependence computation:\n"
  "  -waw <boolean>       Consider WAW (output) dependences (1) or not (0)\n"
  "                       (default setting: 1).\n"
  "  -raw <boolean>       Consider RAW (flow)   dependences (1) or not (0)\n"
  "                       (default setting: 1).\n"
  "  -war <boolean>       Consider WAR (anti)   dependences (1) or not (0)\n"
  "                       (default setting: 1).\n"
  "  -rar <boolean>       Consider RAR (input)  dependences (1) or not (0)\n"
  "                       (default setting: 0).\n");
  printf(
  "  -commute   <boolean> Consider commutativity (1) or not (0)\n"
  "                       (default setting: 0).\n"
  "  -fullcheck <boolean> Compute all legality violation (1) or only one (0)\n"
  "                       (default setting: 0).\n"
  "  -depgraph  <boolean> Ask to print the dependence graph (1) or not (0)\n"
  "                       (default setting: 0 but 1 if no transformation).\n"
  "  -violgraph <boolean> Ask to print the violation graph (1) or not (0)\n"
  "                       (default setting: 1 but 0 if no transformation).\n"
  "  -scalren   <boolean> Ask to enable scalar renaming (1) or not (0)\n"
  "                       (default setting: 0).\n"
  "  -scalpriv  <boolean> Ask to enable scalar privatization (1) or not (0)\n"
  "                       (default setting: 0).\n"
  "  -scalexp   <boolean> Ask to enable scalar expansion (1) or not (0)\n"
  "                       (default setting: 0).\n");
  printf(
  "  -view                Ask to display the graphs (1) or not (0)\n"
  "                       (requires dot -graphviz- and gv tools).\n");
  printf(
  "\nGeneral options:\n"
#ifdef CANDL_SUPPORTS_SCOPLIB
  "  -inscop		  Read a .scop formatted file as the input.\n"
  "  -outscop		  Output a .scop formatted file as the output.\n"
  "  -scoptocandl      	  Output a .candl formatted file from a .scop input.\n"
#endif
  "  -o <output>          Name of the output file; 'stdout' is a special\n"
  "                       value: when used, output is standard output\n"
  "                       (default setting: stdout).\n"
  "  -verbose             Display a verbose output.\n"
  "  -v, --version        Display the version information.\n"
  "  -h, --help           Display this information.\n\n"
  "The special value 'stdin' for 'file' makes Candl to read data on\n"
  "standard input.\n\n"
  "For bug reporting or any suggestions, please send an email to the author\n"
  "<cedric.bastoul@inria.fr>.\n");
}


/**
 * candl_options_version function:
 * This function displays some version informations when the user set the
 * option -version while calling candl. Prints are cutted to respect the 509
 * characters limitation of the ISO C 89 compilers.
 * August 5th 2002: first version.
 */
void candl_options_version()
{ printf("Candl %s %s bits   The Chunky Dependence Analyzer\n",
         CANDL_RELEASE,CANDL_VERSION);
  printf(
  "-----\n"
  "Candl is a dependence analyzer for static control programs, coming from \n"
  "the CHUNKY project: a research tool for data-locality improvement. This \n"
  "program is distributed under the terms of the GNU Lesser General Public\n"
  "License (details at http://www.gnu.org/copyleft/gpl.html).\n"
  "-----\n");
  printf(
  "It would be kind to refer the following paper in any publication "
  "resulting \nfrom the use of this software or its library:\n"
  "@Article{Bas05,\n"
  "author    =  {Cedric Bastoul and Paul Feautrier},\n"
  "title     =  {Adjusting a program transformation for legality},\n"
  "journal   =  {Parallel Processing Letters},\n"
  "year      =  2005,\n"
  "volume    =  15,\n"
  "number    =  1,\n"
  "pages     =  {3--17},\n"
  "month     =  {March},\n"
  "}\n"
  "-----\n"
  "For bug reporting or any suggestions, please send an email to the author\n"
  "<cedric.bastoul@inria.fr>.\n");
}


/**
 * candl_options_set function:
 * This function sets the value of an option thanks to the user's calling line.
 * - option is the value to set,
 * - argc are the elements of the user's calling line,
 * - number is the number of the element corresponding to the considered option,
 *   this function adds 1 to number to pass away the option value.
 * August 5th 2002: first version.
 * June 29th 2003: (debug) lack of argument now detected.
 */
void candl_options_set(int * option, int argc, char ** argv, int * number)
{
  char ** endptr;

  if (*number+1 >= argc)
  {
    fprintf(stderr, "[Candl]ERROR: an option lacks of argument.\n");
    exit(1);
  }

  endptr = NULL;
  *option = strtol(argv[*number+1],endptr,10);
  if (endptr != NULL)
  {
    fprintf(stderr, "[Candl]ERROR: %s option value is not valid.\n",
            argv[*number]);
    exit(1);
  }
  *number = *number + 1;
}


/**
 * candl_options_read function:
 * This functions reads all the options and the input/output files thanks
 * the the user's calling line elements (in argc). It fills a CandlOptions
 * structure and the FILE structure corresponding to input and output files.
 * August 5th 2002: first version.
 * April 19th 2003: now in options.c and support of the CandlOptions structure.
 */
void candl_options_read(int argc, char** argv, FILE** input, FILE** output,
			CandlOptions** options)
{
  int i, infos = 0, input_is_set = 0;

  /* CandlOptions structure allocation and initialization. */
  *options = candl_options_malloc();
  /* The default output is the standard output. */
  *output = stdout;

  for (i = 1; i < argc; i++)
    if (argv[i][0] == '-')
      {
	if (!strcmp(argv[i], "-waw"))
	  candl_options_set(&(*options)->waw, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-raw"))
	  candl_options_set(&(*options)->raw, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-war"))
	  candl_options_set(&(*options)->war, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-rar"))
	  candl_options_set(&(*options)->rar, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-commute"))
	  candl_options_set(&(*options)->commute, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-fullcheck"))
	  candl_options_set(&(*options)->fullcheck, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-depgraph"))
	  candl_options_set(&(*options)->depgraph, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-violgraph"))
	  candl_options_set(&(*options)->violgraph, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-scalren"))
	  candl_options_set(&(*options)->scalar_renaming, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-scalpriv"))
	  candl_options_set(&(*options)->scalar_privatization, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-scalexp"))
	  candl_options_set(&(*options)->scalar_expansion, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-lastwriter"))
	  candl_options_set(&(*options)->lastwriter, argc, argv, &i);
	else
	if (!strcmp(argv[i], "-view"))
	  (*options)->view = 1;
	else
	if (!strcmp(argv[i], "-verbose"))
	  (*options)->verbose = 1;
	else
	if (!strcmp(argv[i], "-inscop"))
	  (*options)->readscop = 1;
	else
	if (!strcmp(argv[i], "-outscop"))
	  (*options)->writescop = 1;
	else
	if (!strcmp(argv[i], "-scoptocandl"))
	  (*options)->scoptocandl = 1;
	else
	if ((!strcmp(argv[i], "-struct")) ||
	    (!strcmp(argv[i], "-structure")))
	  (*options)->structure = 1;
	else
	if ((!strcmp(argv[i], "--help")) || (!strcmp(argv[i], "-h")))
	  {
	    candl_options_help();
	    infos = 1;
	  }
	else
	if ((!strcmp(argv[i], "--version")) || (!strcmp(argv[i], "-v")))
	  {
	    candl_options_version();
	    infos = 1;
	  }
	else
	if (!strcmp(argv[i], "-o"))
	  {
	    if (i+1 >= argc)
	      {
		fprintf(stderr,
			"[Candl]ERROR: no output name for -o option.\n");
		exit(1);
	      }

	    /* stdout is a special value, when used, we set output to standard
	     * output.
	     */
	    if (!strcmp(argv[i+1], "stdout"))
	      *output = stdout;
	    else
	      {
		*output = fopen(argv[i+1], "w");
		if (*output == NULL)
		  {
		    fprintf(stderr,
			    "[Candl]ERROR: can't create output file %s.\n",
			    argv[i+1]);
		    exit(1);
		  }
	      }
	    i++;
	  }
	else
	  fprintf(stderr,  "[Candl]ERROR: unknown %s option.\n", argv[i]);
      }
  else
  {
    if (!input_is_set)
      {
	input_is_set = 1;
	/* stdin is a special value, when used, we set input to
	   standard input. */
	if (!strcmp(argv[i], "stdin"))
	  *input = stdin;
	else
	  {
	    *input = fopen(argv[i], "r");
	    if (*input == NULL)
	      {
		fprintf(stderr,
			"[Candl]ERROR: %s file does not exist.\n", argv[i]);
		exit(1);
	      }
	  }
      }
    else
    {
      fprintf(stderr, "[Candl]ERROR: multiple input files.\n");
      exit(1);
    }
  }
  if (!input_is_set)
    {
      if (!infos)
	fprintf(stderr, "[Candl]ERROR: no input file (-h for help).\n");
      exit(1);
    }
}

