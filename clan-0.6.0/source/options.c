
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                options.c                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 24/05/2008                     **
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
# include <string.h>
# include <clan/options.h>
# include <clan/clan.h>


/*+****************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * clan_option_print function:
 * This function prints the content of a clan_options_t structure (program) into
 * a file (foo, possibly stdout).
 * \param foo     File where informations are printed.
 * \param options Option structure whose information have to be printed.
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
void
clan_options_print(FILE * foo, clan_options_p options)
{
  fprintf(foo,"Options:\n");

  if (options->name != NULL)
    fprintf(foo,"name            = %s,\n",options->name);
  else
    fprintf(foo,"name            = NULL,\n");

  fprintf(foo,"castle          = %3d,\n",options->castle);
  fprintf(foo,"structure       = %3d.\n",options->structure);
  fprintf(foo,"inputscop       = %3d.\n",options->inputscop);
  fprintf(foo,"bounded_context = %3d.\n",options->bounded_context);
  fprintf(foo,"symboltable     = %3d.\n",options->symboltable);
  fprintf(foo,"debug           = %3d.\n",options->debug);
}


/*+****************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * clan_options_free function:
 * This function frees the allocated memory for a clan_options_t structure.
 * \param options Option structure to be freed.
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
void
clan_options_free(clan_options_p options)
{
  free(options);
}


/*+****************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/


/**
 * clan_options_help function:
 * This function displays the quick help when the user set the option -help
 * while calling clan. Prints are cut to respect the 509 characters
 * limitation of the ISO C 89 compilers.
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
void
clan_options_help()
{
  printf(
  "Usage: clan [ options | file ] ...\n");
  printf(
  "\nGeneral options:\n"
  "  -o <output>           Name of the output file; 'stdout' is a special\n"
  "                        value: when used, output is standard output\n"
  "                        (default setting: stdout).\n"
  "  -inputscop            Read a .scop as the input.\n"
  "  -arraystag            Dump the arrays table in the <arrays> tag.\n"
  "  -boundedctxt          Bound all global parameters to be >= -1.\n"
  "  -symboltable          Display the symbol table .\n"
  "  -debug                Display the debug information.\n"
  "  -v, --version         Display the release information (and more).\n"
  "  -h, --help            Display this information.\n\n");
  printf(
  "The special value 'stdin' for 'file' or the special option '-' makes clan\n"
  "to read data on standard input.\n\n"
  "For bug reporting or any suggestions, please send an email to the author\n"
  "Cedric Bastoul <cedric.bastoul@inria.fr>.\n");
}


/**
 * clan_options_version function:
 * This function displays some version informations when the user set the
 * option -version while calling clan. Prints are cut to respect the 509
 * characters limitation of the ISO C 89 compilers.
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
void
clan_options_version()
{
  printf("clan %s %s bits       The Chunky Loop Analyzer\n",
         CLAN_RELEASE,CLAN_VERSION);
  printf(
  "-----\n"
  "This is a polyhedral representation extractor for imperative programs using "
  "a C\ngrammar for control flow and array accesses (this includes C, C++,"
  " Java, C#\nand probably few toy languages too). This program is distributed "
  "under the\nterms of the GNU Lesser General Public License, see details of "
  "the licence at\nhttp://www.gnu.org/copyleft/lgpl.html\n"
  "-----\n");
  printf(
  "It would be fair to refer the following paper in any publication "
  "resulting from\nthe use of this software or its library (it defines SCoPs):\n"
  "@InProceedings{Bas03,\n"
  "author    =  {Cedric Bastoul and Albert Cohen and Sylvain Girbal and\n"
  "              Saurabh Sharma and Olivier Temam},\n"
  "title     =  {Putting Polyhedral Loop Transformations to Work},\n"
  "booktitle =  {LCPC'16 International Workshop on Languages and\n"
  "              Compilers for Parallel Computers, LNCS 2958},\n"
  "pages     =  {209--225},\n"
  "month     =  {october},\n"
  "year      =  2003,\n"
  "address   =  {College Station, Texas}\n"
  "}\n"
  "-----\n");
  printf(
  "For any information, please send an email to the author\n"
  "Cedric Bastoul <cedric.bastoul@inria.fr>.\n");
}


/**
 * clan_options_set function:
 * This function sets the value of an option thanks to the user's calling line.
 * \param option The value to set,
 * \param argv   Number of elements in the user's calling line,
 * \param argc   Elements of the user's calling line,
 * \param number Number of the element corresponding to the considered option,
 *               this function adds 1 to number to pass away the option value.
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
void
clan_options_set(int * option, int argv, char ** argc, int * number)
{
  char ** endptr;

  if (*number+1 >= argv)
  { fprintf(stderr, "[clan]ERROR: an option lacks of argument.\n");
    exit(1);
  }

  endptr = NULL;
  *option = strtol(argc[*number+1],endptr,10);
  if (endptr != NULL)
  { fprintf(stderr, "[clan]ERROR: %s value for %s option is not valid.\n",
            argc[*number+1],argc[*number]);
    exit(1);
  }
  *number = *number + 1;
}


/**
 * clan_options_malloc function:
 * This functions allocate the memory space for a clan_options_t structure and
 * fill its fields with the defaults values. It returns a pointer to the
 * allocated clan_options_t structure.
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
clan_options_p
clan_options_malloc(void)
{
  clan_options_p options;

  /* Memory allocation for the clan_options_t structure. */
  options = (clan_options_p)malloc(sizeof(clan_options_t));
  if (options == NULL)
  { fprintf(stderr, "[clan]ERROR: memory overflow.\n");
    exit(1);
  }

  /* We set the various fields with default values. */
  options->name      = NULL; /* Name of the input file is not set. */
  options->castle    = 1;    /* Do print the Clan McCloog castle in output. */
  options->structure = 0;    /* Don't print internal structure.*/
  options->inputscop = 0;    /* Default input is a source file, not a .scop.*/
  options->arraystag = 0;    /* Don't dump the array list in the
				<arrays> tag. */
  options->bounded_context = 0;/* Don't bound the global parameters. */
  options->symboltable = 0;
  options->debug = 0;
  return options;
}


/**
 * clan_options_read function:
 * This functions reads all the options and the input/output files thanks
 * the the user's calling line elements (in argc). It fills a clan_options_t
 * structure and the FILE structure corresponding to input and output files.
 * \param argv    Number of strings in command line.
 * \param argc    Array of command line strings.
 * \param input   Input  file (modified by the function).
 * \param output  Output file (modified by the function).
 **
 * - 24/05/2008: first version (from CLooG 0.14.1).
 */
clan_options_p
clan_options_read(int argv, char ** argc, FILE ** input, FILE ** output)
{
  int i, infos=0, input_is_set=0;
  clan_options_p options;

  /* clan_options_t structure allocation and initialization. */
  options = clan_options_malloc();

  /* The default output is the standard output. */
  *output = stdout;

  for (i=1; i < argv; i++)
  {
    if (argc[i][0] == '-')
    {
      if (argc[i][1] == '\0')
      {
        /* "-" alone is a special option to set input to standard input. */
        input_is_set = 1;
	*input = stdin;
      }
      else
      if (strcmp(argc[i],"-castle") == 0)
        clan_options_set(&(options)->castle,argv,argc,&i);
      else
      if (strcmp(argc[i],"-structure") == 0)
        options->structure = 1;
      else
      if (strcmp(argc[i],"-inputscop") == 0)
        options->inputscop = 1;
      else
      if (strcmp(argc[i],"-arraystag") == 0)
        options->arraystag = 1;
      else
      if (strcmp(argc[i],"-boundedctxt") == 0)
        options->bounded_context = 1;
      else
      if (strcmp(argc[i],"-symboltable") == 0)
        options->symboltable = 1;        
      else
      if (strcmp(argc[i],"-debug") == 0)
        options->debug = 1;                
      else
      if ((strcmp(argc[i],"--help") == 0) || (strcmp(argc[i],"-h") == 0))
      {
        clan_options_help();
        infos = 1;
      }
      else
      if ((strcmp(argc[i],"--version") == 0) || (strcmp(argc[i],"-v") == 0))
      {
        clan_options_version();
        infos = 1;
      }
      else
      if (strcmp(argc[i],"-o") == 0)
      {
        if (i+1 >= argv)
        {
	  fprintf(stderr, "[clan]ERROR: no output name for -o option.\n");
          exit(1);
        }

        /* stdout is a special value to set output to standard output. */
        if (strcmp(argc[i+1],"stdout") == 0)
          *output = stdout;
        else
        {
	  *output = fopen(argc[i+1],"w");
          if (*output == NULL)
          {
	    fprintf(stderr, "[clan]ERROR: can't create output file %s.\n",
	            argc[i+1]);
            exit(1);
          }
        }
        i++;
      }
      else
        fprintf(stderr, "[clan]WARNING: unknown %s option.\n",argc[i]);
    }
    else
    {
      if (!input_is_set)
      {
        input_is_set = 1;
        options->name = argc[i];
        /* stdin is a special value to set input to standard input. */
        if (strcmp(argc[i],"stdin") == 0)
          *input = stdin;
        else
        {
	  *input = fopen(argc[i],"r");
          if (*input == NULL)
          {
	    fprintf(stderr, "[clan]ERROR: %s file does not exist.\n",argc[i]);
            exit(1);
          }
        }
      }
      else
      {
        fprintf(stderr, "[clan]ERROR: multiple input files.\n");
        exit(1);
      }
    }
  }

  if (!input_is_set)
  {
    if (!infos)
      fprintf(stderr, "[clan]ERROR: no input file (-h for help).\n");
    clan_options_free(options);
    exit(1);
  }

  return options;
}
