
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                            options.c                              **
    **-------------------------------------------------------------------**
    **                  First version: april 19th 2003                   **
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


#include <stdarg.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include "../include/cloog/cloog.h"


/******************************************************************************
 *                          Error reporting functions                         *
 ******************************************************************************/

void cloog_vmsg(CloogOptions *options, enum cloog_msg_type type,
		const char *msg, va_list ap)
{
  const char *type_msg;

  if (options && options->quiet &&
      (type == CLOOG_WARNING || type == CLOOG_INFO))
    return;

  switch(type) {
  case CLOOG_WARNING:
	type_msg = "WARNING";
	break;
  case CLOOG_INFO:
	type_msg = "INFO";
	break;
  case CLOOG_ERROR:
  default:
	type_msg = "ERROR";
	break;
  }
  fprintf(stderr, "[CLooG] %s: ", type_msg);
  vfprintf(stderr, msg, ap);
}

/**
 * Print message to stderr.
 * @param msg printf format string
 */
void cloog_msg(CloogOptions *options, enum cloog_msg_type type,
		const char *msg, ...)
{
  va_list args;

  va_start(args, msg);
  cloog_vmsg(options, type, msg, args);
  va_end(args);
}

/**
 * Print error message to stderr and exit.
 * @param msg printf format string
 */
void cloog_die(const char *msg, ...)
{
  va_list args;

  va_start(args, msg);
  cloog_vmsg(NULL, CLOOG_ERROR, msg, args);
  va_end(args);
  exit(1);
}

/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_option_print function:
 * This function prints the content of a CloogOptions structure (program) into
 * a file (foo, possibly stdout).
 * - April 19th 2003: first version.
 */
void cloog_options_print(FILE * foo, CloogOptions * options)
{ fprintf(foo,"Options:\n") ;
  fprintf(foo,"OPTIONS FOR LOOP GENERATION\n") ;
  fprintf(foo,"l           = %3d,\n",options->l) ;
  fprintf(foo,"f           = %3d,\n",options->f) ;
  fprintf(foo,"stop        = %3d,\n",options->stop) ;
  fprintf(foo,"strides     = %3d,\n",options->strides) ;
  fprintf(foo,"sh          = %3d,\n",options->sh);
  fprintf(foo,"OPTIONS FOR PRETTY PRINTING\n") ;
  fprintf(foo,"esp         = %3d,\n",options->esp) ;
  fprintf(foo,"fsp         = %3d,\n",options->fsp) ;
  fprintf(foo,"otl         = %3d.\n",options->otl) ;
  fprintf(foo,"block       = %3d.\n",options->block) ;
  fprintf(foo,"compilable  = %3d.\n",options->compilable) ;
  fprintf(foo,"callable    = %3d.\n",options->callable) ;
  fprintf(foo,"UNDOCUMENTED OPTIONS FOR THE AUTHOR ONLY\n") ;
  fprintf(foo,"leaks       = %3d.\n",options->leaks) ;
  fprintf(foo,"nobacktrack = %3d.\n",options->nobacktrack) ;
  fprintf(foo,"override    = %3d.\n",options->override) ;
  fprintf(foo,"structure   = %3d.\n",options->structure) ;
  fprintf(foo,"noscalars   = %3d.\n",options->noscalars) ;
  fprintf(foo,"noblocks    = %3d.\n",options->noblocks) ;
  fprintf(foo,"nosimplify  = %3d.\n",options->nosimplify) ;
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_options_free function:
 * This function frees the allocated memory for a CloogOptions structure.
 * - April 19th 2003: first version.
 */
void cloog_options_free(CloogOptions *options)
{
  free(options);
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * cloog_options_help function:
 * This function displays the quick help when the user set the option -help
 * while calling cloog. Prints are cutted to respect the 509 characters
 * limitation of the ISO C 89 compilers.
 * - August 5th 2002: first version.
 */
void cloog_options_help()
{ printf(
  "Usage: cloog [ options | file ] ...\n"
  "Options for code generation:\n"
  "  -l <depth>            Last loop depth to optimize (-1: infinity)\n"
  "                        (default setting: -1).\n"
  "  -f <depth>            First loop depth to start loop separation (-1: "
  "infinity)\n                        (default setting:  1).\n") ;
  printf(
  "  -stop <depth>         Loop depth to stop code generation (-1: infinity)"
  "\n                        (default setting: -1).\n"
  "  -strides <boolean>    Handle non-unit strides (1) or not (0)\n"
  "                        (default setting:  0).\n") ;
  printf(
  "\nOptions for pretty printing:\n"
  "  -otl <boolean>        Simplify loops running one time (1) or not (0)\n"
  "                        (default setting:  1).\n") ;
  printf(
  "  -esp <boolean>        Allow complex equalities spreading (1) or not (0)\n"
  "                        (default setting:  0).\n");
  printf(
  "  -fsp <level>          First level to begin the spreading\n"
  "                        (default setting:  1).\n"
  "  -block <boolean>      Make a new statement block per iterator in C\n"
  "                        programs (1) or not (0) (default setting: 0).\n") ;
  printf(
  "  -compilable <number>  Compilable code by using preprocessor (not 0) or" 
  "\n                        not (0), number being the value of the parameters"
  "\n                        (default setting:  0).\n"
  "  -callable <boolean>   Testable code by using preprocessor (not 0) or" 
  "\n                        not (0) (default setting:  0).\n");
  printf(
  "\nGeneral options:\n"
  "  -o <output>           Name of the output file; 'stdout' is a special\n"
  "                        value: when used, output is standard output\n"
  "                        (default setting: stdout).\n"
  "  -v, --version         Display the version information (and more).\n"
  "  -q, --quiet           Don't print any informational messages.\n"
  "  -h, --help            Display this information.\n\n") ;
  printf(
  "The special value 'stdin' for 'file' makes CLooG to read data on\n"
  "standard input.\n\n"
  "For bug reporting or any suggestions, please send an email to the author\n"
  "<cedric.bastoul@inria.fr>.\n") ;
}


/**
 * cloog_options_version function:
 * This function displays some version informations when the user set the
 * option -version while calling cloog. Prints are cutted to respect the 509
 * characters limitation of the ISO C 89 compilers.
 * - August 5th 2002: first version.
 */
void cloog_options_version()
{ printf("%s       The Chunky Loop Generator\n", cloog_version());
  printf(
  "-----\n"
  "This is a loop generator for scanning Z-polyhedra. It is based on the "
  "work of\nF. Quillere and C. Bastoul on high level code generation and of "
  "the PolyLib Team\non polyhedral computation. This program is distributed "
  "under the terms of the\nGNU Lesser General Public License "
  "(details at http://www.gnu.org/licenses/lgpl-2.1.html).\n"
  "-----\n") ;
  printf(
  "It would be fair to refer the following paper in any publication "
  "resulting from\nthe use of this software or its library:\n"
  "@InProceedings{Bas04,\n"
  "author    =  {Cedric Bastoul},\n"
  "title     =  {Code Generation in the Polyhedral Model Is Easier Than You "
  "Think},\n"
  "booktitle =  {PACT'13 IEEE International Conference on Parallel "
  "Architecture\n             and Compilation Techniques},\n"
  "pages     =  {7--16},\n"
  "month     =  {september},\n"
  "year      =  2004,\n"
  "address   =  {Juan-les-Pins}\n"
  "}\n"
  "-----\n"
  "For any information, please ask the author at "
  "<cedric.bastoul@inria.fr>.\n") ;
} 


/**
 * cloog_options_set function:
 * This function sets the value of an option thanks to the user's calling line.
 * - option is the value to set,
 * - argc are the elements of the user's calling line,
 * - number is the number of the element corresponding to the considered option,
 *   this function adds 1 to number to pass away the option value.
 **
 * - August 5th 2002: first version.
 * - June 29th 2003: (debug) lack of argument now detected.
 */
void cloog_options_set(int * option, int argv, char ** argc, int * number)
{ char ** endptr ;
  
  if (*number+1 >= argv)
    cloog_die("an option lacks of argument.\n");

  endptr = NULL ;
  *option = strtol(argc[*number+1],endptr,10) ;
  if (endptr != NULL)
    cloog_die("value '%s' for option '%s' is not valid.\n",
	      argc[*number+1], argc[*number]);
  *number = *number + 1 ;
}


/**
 * cloog_options_malloc function:
 * This functions allocate the memory space for a CLoogOptions structure and
 * fill its fields with the defaults values. It returns a pointer to the
 * allocated CloogOptions structure.
 * - April    19th 2003: first version.
 * - November 21th 2005: name changed (before it was cloog_options_init).
 */
CloogOptions *cloog_options_malloc(CloogState *state)
{ CloogOptions * options ;

  /* Memory allocation for the CloogOptions structure. */
  options = (CloogOptions *)malloc(sizeof(CloogOptions)) ;
  if (options == NULL) 
    cloog_die("memory overflow.\n");
  
  options->state = state;

  /* We set the various fields with default values. */
  /* OPTIONS FOR LOOP GENERATION */
  options->l           = -1 ;  /* Last level to optimize: infinity. */
  options->f           =  1 ;  /* First level to optimize: the first. */
  options->stop        = -1 ;  /* Generate all the code. */
  options->strides     =  0 ;  /* Generate a code with unit strides. */
  options->sh	       =  0;   /* Compute actual convex hull. */
  options->name	       = "";
  /* OPTIONS FOR PRETTY PRINTING */
  options->esp         =  1 ;  /* We don't want Equality SPreading.*/
  options->fsp         =  1 ;  /* The First level to SPread is the first. */
  options->otl         =  1 ;  /* We want to fire One Time Loops. */
  options->block       =  0 ;  /* We don't want to force statement blocks. */
  options->compilable  =  0 ;  /* No compilable code. */
  options->callable    =  0 ;  /* No callable code. */
  options->quiet       =  0;   /* Do print informational messages. */
  /* UNDOCUMENTED OPTIONS FOR THE AUTHOR ONLY */
  options->leaks       =  0 ;  /* I don't want to print allocation statistics.*/
  options->nobacktrack =  0 ;  /* No backtrack in Quillere's algorithm.*/
  options->override    =  0 ;  /* I don't want to override CLooG decisions.*/
  options->structure   =  0 ;  /* I don't want to print internal structure.*/
  options->noblocks    =  0 ;  /* I do want to make statement blocks.*/
  options->noscalars   =  0 ;  /* I do want to use scalar dimensions.*/
  options->nosimplify  =  0 ;  /* I do want to simplify polyhedra.*/
  
  return options ;
}



/**
 * cloog_options_read function:
 * This functions reads all the options and the input/output files thanks
 * the the user's calling line elements (in argc). It fills a CloogOptions
 * structure and the FILE structure corresponding to input and output files.
 * - August 5th 2002: first version.
 * - April 19th 2003: now in options.c and support of the CloogOptions structure.
 */
void cloog_options_read(CloogState *state, int argc, char **argv,
			FILE **input, FILE **output, CloogOptions **options)
{ int i, infos=0, input_is_set=0 ;
  
  /* CloogOptions structure allocation and initialization. */
  *options = cloog_options_malloc(state);
  
  /* The default output is the standard output. */
  *output = stdout ;

  for (i=1;i<argc;i++)
  if (argv[i][0] == '-')
  { if (strcmp(argv[i],"-l")   == 0)
    cloog_options_set(&(*options)->l,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-f")   == 0)
    cloog_options_set(&(*options)->f,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-stop")   == 0)
    cloog_options_set(&(*options)->stop,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-strides")   == 0)
    cloog_options_set(&(*options)->strides,argc,argv,&i) ;
    else if (strcmp(argv[i],"-sh")   == 0)
      cloog_options_set(&(*options)->sh,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-otl") == 0)
    cloog_options_set(&(*options)->otl,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-esp") == 0)
    cloog_options_set(&(*options)->esp,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-fsp") == 0)
    cloog_options_set(&(*options)->fsp,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-block") == 0)
    cloog_options_set(&(*options)->block,argc,argv,&i) ;
    else
    if (strcmp(argv[i],"-compilable") == 0)
      cloog_options_set(&(*options)->compilable, argc, argv, &i);
    else if (strcmp(argv[i], "-callable") == 0)
      cloog_options_set(&(*options)->callable, argc, argv, &i);
    else
    if (strcmp(argv[i],"-loopo") == 0) /* Special option for the LooPo team ! */
    { (*options)->esp   = 0 ;
      (*options)->block = 1 ;
    }
    else
    if (strcmp(argv[i],"-bipbip") == 0)/* Special option for the author only !*/
      (*options)->nobacktrack = 1 ;
    else
    if (strcmp(argv[i],"-leaks") == 0)
    (*options)->leaks = 1 ;
    else
    if (strcmp(argv[i],"-nobacktrack") == 0)
    (*options)->nobacktrack = 1 ;
    else
    if (strcmp(argv[i],"-override") == 0)
    (*options)->override = 1 ;
    else
    if (strcmp(argv[i],"-noblocks") == 0)
    (*options)->noblocks = 1 ;
    else
    if (strcmp(argv[i],"-noscalars") == 0)
    (*options)->noscalars = 1 ;
    else
    if (strcmp(argv[i],"-nosimplify") == 0)
    (*options)->nosimplify = 1 ;
    else
    if ((strcmp(argv[i],"-struct") == 0) || (strcmp(argv[i],"-structure") == 0))
    (*options)->structure = 1 ;
    else
    if ((strcmp(argv[i],"--help") == 0) || (strcmp(argv[i],"-h") == 0))
    { cloog_options_help() ;
      infos = 1 ;
    }
    else
    if ((strcmp(argv[i],"--version") == 0) || (strcmp(argv[i],"-v") == 0))
    { cloog_options_version() ;
      infos = 1 ;
    } else if ((strcmp(argv[i],"--quiet") == 0) || (strcmp(argv[i],"-q") == 0))
      (*options)->quiet = 1;
    else
    if (strcmp(argv[i],"-o") == 0)
    { if (i+1 >= argc)
        cloog_die("no output name for -o option.\n");

      /* stdout is a special value, when used, we set output to standard
       * output.
       */
      if (strcmp(argv[i+1],"stdout") == 0)
      *output = stdout ;
      else
      { *output = fopen(argv[i+1],"w") ;
        if (*output == NULL)
          cloog_die("can't create output file %s.\n", argv[i+1]);
      }
      i ++ ;    
    }
    else
      cloog_msg(*options, CLOOG_WARNING, "unknown %s option.\n", argv[i]);
  }
  else
  { if (!input_is_set)
    { input_is_set = 1 ;
      (*options)->name = argv[i] ;
      /* stdin is a special value, when used, we set input to standard input. */
      if (strcmp(argv[i],"stdin") == 0)
      *input = stdin ;
      else
      { *input = fopen(argv[i],"r") ;
        if (*input == NULL)
          cloog_die("%s file does not exist.\n", argv[i]);
      }
    } 
    else
      cloog_die("multiple input files.\n");
  }
  if (!input_is_set)
  { if (!infos)
      cloog_die("no input file (-h for help).\n");
    exit(1) ;
  }
}

