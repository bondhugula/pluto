
   /**-------------------------------------------------------------------**
    **                              CLooG                                **
    **-------------------------------------------------------------------**
    **                            program.c                              **
    **-------------------------------------------------------------------**
    **                 First version: october 25th 2001                  **
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
/* CAUTION: the english used for comments is probably the worst you ever read,
 *          please feel free to correct and improve it !
 */


# include <sys/types.h>
# include <sys/time.h>
#include <stdarg.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <ctype.h>
# include <unistd.h>
# include "../include/cloog/cloog.h"
#ifdef CLOOG_RUSAGE
# include <sys/resource.h>
#endif


/******************************************************************************
 *                          Structure display function                        *
 ******************************************************************************/


/**
 * cloog_program_print function:
 * this function is a human-friendly way to display the CloogProgram data
 * structure, it shows all the different fields and includes an indentation
 * level (level) in order to work with others print_structure functions.
 * - July 1st 2005: first version based on the old cloog_program_print function.
 */
void cloog_program_print_structure(file, program, level)
FILE * file ; 
CloogProgram * program ;
int level ;
{ int i, j ;

  /* Go to the right level. */
  for (i=0; i<level; i++)
  fprintf(file,"|\t") ;
  
  fprintf(file,"+-- CloogProgram\n") ;
  
  /* A blank line. */
  for (i=0; i<=level+1; i++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
  
  /* Print the language. */
  for (i=0; i<=level; i++)
  fprintf(file,"|\t") ;
  fprintf(file, "Language: %c\n",program->language) ;
  
  /* A blank line. */
  for (i=0; i<=level+1; i++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;

  /* Print the scattering dimension number. */
  for (i=0; i<=level; i++)
  fprintf(file,"|\t") ;
  fprintf(file,"Scattering dimension number: %d\n",program->nb_scattdims) ;
  
  /* A blank line. */
  for (i=0; i<=level+1; i++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
  
  /* Print the scalar scattering dimension informations. */
  for (i=0; i<=level; i++)
  fprintf(file,"|\t") ;
  if (program->scaldims != NULL)
  { fprintf(file,"Scalar dimensions:") ;
    for (i=0;i<program->nb_scattdims;i++)
    fprintf(file," %d:%d ",i,program->scaldims[i]) ;
    fprintf(file,"\n") ;
  }
  else
  fprintf(file,"No scalar scattering dimensions\n") ;
  
  /* A blank line. */
  for (i=0; i<=level+1; i++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;

  /* Print the parameter and the iterator names. */
  cloog_names_print_structure(file,program->names,level+1) ;
 
  /* A blank line. */
  for (i=0; i<=level+1; i++)
  fprintf(file,"|\t") ;
  fprintf(file,"\n") ;
  
  /* Print the context. */
  cloog_domain_print_structure(file, program->context, level+1, "Context");
    
  /* Print the loop. */
  cloog_loop_print_structure(file,program->loop,level+1) ;

  /* One more time something that is here only for a better look. */
  for (j=0; j<2; j++)
  { for (i=0; i<=level; i++)
    fprintf(file,"|\t") ;
      
    fprintf(file,"\n") ;
  }
}


/**
 * cloog_program_dump_cloog function:
 * This function dumps a CloogProgram structure supposed to be completely
 * filled in a CLooG input file (foo possibly stdout) such as CLooG can
 * rebuild almost exactly the data structure from the input file (the number
 * of scattering functions is lost since they are included inside the
 * iteration domains, this can only lead to a less beautiful pretty printing).
 * WARNING: this function do not respect CloogDomain as an object.
 * - June     27th 2003: first version.
 * - May      15th 2005: (debug) several patches by Kristof Beyls.
 * - November 16th 2005: adaptation for CLooG 0.14.0 structures.
 */
void cloog_program_dump_cloog(FILE * foo, CloogProgram * program)
{
  int i;
  CloogLoop * loop ;

  fprintf(foo,
  "# CLooG -> CLooG\n"
  "# This is an automatic dump of a CLooG input file from a CloogProgram data\n"
  "# structure. WARNING: it is highly dangerous and MAY be correct ONLY if\n"
  "# - it has been dumped before loop generation.\n"
  "# - option -noscalars is used (it removes scalar dimensions otherwise)\n"
  "# - option -l is at least the original scattering dimension number\n"
  "# ASK THE AUTHOR IF YOU *NEED* SOMETHING MORE ROBUST\n") ;

  /* Language. */
  if (program->language == 'c')
  fprintf(foo,"# Language: C\n") ;
  else
  fprintf(foo,"# Language: FORTRAN\n") ;
  fprintf(foo,"%c\n\n",program->language) ;

  /* Context. */
  fprintf(foo, "# Context (%d parameter(s)):\n", program->names->nb_parameters);
  cloog_domain_print_constraints(foo, program->context, 0);
  fprintf(foo,"1 # Parameter name(s)\n") ;
  for (i=0;i<program->names->nb_parameters;i++)
  fprintf(foo,"%s ",program->names->parameters[i]) ;

  /* Statement number. */
  i = 0 ;
  loop = program->loop ;
  while (loop != NULL)
  { i++ ;
    loop = loop->next ;
  }
  fprintf(foo,"\n\n# Statement number:\n%d\n\n",i) ;

  /* Iteration domains. */
  i = 1 ;
  loop = program->loop ;
  while (loop != NULL)
  { /* Name of the domain. */
    fprintf(foo,"# Iteration domain of statement %d.\n",i) ;

    cloog_domain_print_constraints(foo, loop->domain, 1);
    fprintf(foo,"0 0 0 # For future options.\n\n") ;
    
    i++ ;
    loop = loop->next ;
  }
  fprintf(foo,"\n1 # Iterator name(s)\n") ;
  for (i=0;i<program->names->nb_scattering;i++)
    fprintf(foo,"%s ",program->names->scattering[i]);
  for (i=0;i<program->names->nb_iterators;i++)
    fprintf(foo,"%s ",program->names->iterators[i]);
  fprintf(foo,"\n\n") ;

  /* Scattering functions (none since included inside domains). */
  fprintf(foo,"# No scattering functions.\n0\n\n") ;
}


/**
 * cloog_program_print function:
 * This function prints the content of a CloogProgram structure (program) into a
 * file (file, possibly stdout).
 * - July 1st 2005: Now this very old function (probably as old as CLooG) is
 *                  only a frontend to cloog_program_print_structure, with a
 *                  quite better human-readable representation.
 */
void cloog_program_print(FILE * file, CloogProgram * program)
{ cloog_program_print_structure(file,program,0) ;
}


static void print_comment(FILE *file, CloogOptions *options,
			  const char *fmt, ...)
{
  va_list args;

  va_start(args, fmt);
  if (options->language == LANGUAGE_FORTRAN) {
    fprintf(file, "! ");
    vfprintf(file, fmt, args);
    fprintf(file, "\n");
  } else {
    fprintf(file, "/* ");
    vfprintf(file, fmt, args);
    fprintf(file, " */\n");
  }
}

static void print_macros(FILE *file)
{
    fprintf(file, "/* Useful macros. */\n") ;
    fprintf(file,
	"#define floord(n,d) (((n)<0) ? -((-(n)+(d)-1)/(d)) : (n)/(d))\n");
    fprintf(file,
	"#define ceild(n,d)  (((n)<0) ? -((-(n))/(d)) : ((n)+(d)-1)/(d))\n");
    fprintf(file, "#define max(x,y)    ((x) > (y) ? (x) : (y))\n") ; 
    fprintf(file, "#define min(x,y)    ((x) < (y) ? (x) : (y))\n\n") ; 
}

static void print_declarations(FILE *file, int n, char **names)
{
    int i;

    fprintf(file, "  int %s", names[0]); 
    for (i = 1; i < n; i++)
	fprintf(file, ", %s", names[i]); 

    fprintf(file, ";\n");
}

static void print_iterator_declarations(FILE *file, CloogProgram *program,
	CloogOptions *options)
{
    CloogNames *names = program->names;

    if (names->nb_scattering) {
	fprintf(file, "  /* Scattering iterators. */\n");
	print_declarations(file, names->nb_scattering, names->scattering);
    }
    if (names->nb_iterators) {
	fprintf(file, "  /* Original iterators. */\n");
	print_declarations(file, names->nb_iterators, names->iterators);
    }
}

static void print_callable_preamble(FILE *file, CloogProgram *program,
	CloogOptions *options)
{
    int j;
    CloogBlockList *blocklist;
    CloogBlock *block;
    CloogStatement *statement;

    fprintf(file, "extern void hash(int);\n\n");

    print_macros(file);

    for (blocklist = program->blocklist; blocklist; blocklist = blocklist->next) {
	block = blocklist->block;
	for (statement = block->statement; statement; statement = statement->next) {
	    fprintf(file, "#define S%d(", statement->number);
	    if (block->depth > 0) {
		fprintf(file, "%s", program->names->iterators[0]);
		for(j = 1; j < block->depth; j++)
		    fprintf(file, ",%s", program->names->iterators[j]);
	    }
	    fprintf(file,") { hash(%d);", statement->number);
	    for(j = 0; j < block->depth; j++)
		fprintf(file, " hash(%s);", program->names->iterators[j]);
	    fprintf(file, " }\n");
	}
    }
    fprintf(file, "\nvoid test("); 
    if (program->names->nb_parameters > 0) {
	fprintf(file, "int %s", program->names->parameters[0]);
	for(j = 1; j < program->names->nb_parameters; j++)
	    fprintf(file, ", int %s", program->names->parameters[j]);
    }
    fprintf(file, ")\n{\n"); 
    print_iterator_declarations(file, program, options);
}

static void print_callable_postamble(FILE *file, CloogProgram *program)
{
    fprintf(file, "}\n"); 
}

/**
 * cloog_program_pprint function:
 * This function prints the content of a CloogProgram structure (program) into a
 * file (file, possibly stdout), in a C-like language.
 * - June 22nd 2005: Adaptation for GMP.
 */
void cloog_program_pprint(file, program, options)
FILE * file ;
CloogProgram * program ;
CloogOptions * options ;
{ int i, j, nb_scattering, indentation=0 ;
  CloogStatement * statement ;
  CloogBlockList * blocklist ;
  CloogBlock * block ;
  struct clast_stmt *root;
   
  if (program->language == 'f')
    options->language = LANGUAGE_FORTRAN ;
  else
    options->language = LANGUAGE_C ;
 
#ifdef CLOOG_RUSAGE
  print_comment(file, options, "Generated from %s by %s in %.2fs.",
		options->name, cloog_version(), options->time);
#else
  print_comment(file, options, "Generated from %s by %s.",
		options->name, cloog_version());
#endif
#ifdef CLOOG_MEMORY
  print_comment(file, options, "CLooG asked for %d KBytes.", options->memory);
  cloog_msg(CLOOG_INFO, "%.2fs and %dKB used for code generation.\n",
	  options->time,options->memory);
#endif
  
  /* If the option "compilable" is set, we provide the whole stuff to generate
   * a compilable code. This code just do nothing, but now the user can edit
   * the source and set the statement macros and parameters values.
   */
  nb_scattering = program->nb_scattdims ;
  if (options->compilable && (program->language == 'c'))
  { /* The headers. */
    fprintf(file,"/* DON'T FORGET TO USE -lm OPTION TO COMPILE. */\n\n") ;
    fprintf(file,"/* Useful headers. */\n") ;
    fprintf(file,"#include <stdio.h>\n") ;
    fprintf(file,"#include <stdlib.h>\n") ;
    fprintf(file,"#include <math.h>\n\n") ;

    /* The value of parameters. */
    fprintf(file,"/* Parameter value. */\n") ;
    for (i = 1; i <= program->names->nb_parameters; i++)
      fprintf(file, "#define PARVAL%d %d\n", i, options->compilable);
    
    /* The macros. */
    print_macros(file);

    /* The statement macros. */
    fprintf(file,"/* Statement macros (please set). */\n") ;
    blocklist = program->blocklist ;
    while (blocklist != NULL)
    { block = blocklist->block ;
      statement = block->statement ;
      while (statement != NULL)
      { fprintf(file,"#define S%d(",statement->number) ;
        if (block->depth > 0)
        { fprintf(file,"%s",program->names->iterators[0]) ;
          for(j=1;j<block->depth;j++)
          fprintf(file,",%s",program->names->iterators[j]) ;
        }
        fprintf(file,") {total++;") ;
	if (block->depth > 0)
        { fprintf(file," printf(\"S%d \%%d",statement->number) ;
          for(j=1;j<block->depth;j++)
          fprintf(file," \%%d") ;
          
          fprintf(file,"\\n\",%s",program->names->iterators[0]) ;
	  for(j=1;j<block->depth;j++)
          fprintf(file,",%s",program->names->iterators[j]) ;
          fprintf(file,");") ;
        }
        fprintf(file,"}\n") ;
        
	statement = statement->next ;
      }
      blocklist = blocklist->next ;
    }
    
    /* The iterator and parameter declaration. */
    fprintf(file,"\nint main() {\n") ; 
    print_iterator_declarations(file, program, options);
    if (program->names->nb_parameters > 0)
    { fprintf(file,"  /* Parameters. */\n") ;
      fprintf(file, "  int %s=PARVAL1",program->names->parameters[0]);
      for(i=2;i<=program->names->nb_parameters;i++)
        fprintf(file, ", %s=PARVAL%d", program->names->parameters[i-1], i);
      
      fprintf(file,";\n");
    }
    fprintf(file,"  int total=0;\n");
    fprintf(file,"\n") ;
    
    /* And we adapt the identation. */
    indentation += 2 ;
  } else if (options->callable && program->language == 'c') {
    print_callable_preamble(file, program, options);
    indentation += 2;
  }
  
  root = cloog_clast_create(program, options);
  clast_pprint(file, root, indentation, options);
  cloog_clast_free(root);
  
  /* The end of the compilable code in case of 'compilable' option. */
  if (options->compilable && (program->language == 'c'))
  {
    fprintf(file, "\n  printf(\"Number of integral points: %%d.\\n\",total);");
    fprintf(file, "\n  return 0;\n}\n");
  } else if (options->callable && program->language == 'c')
    print_callable_postamble(file, program);
}


/******************************************************************************
 *                         Memory deallocation function                       *
 ******************************************************************************/


/**
 * cloog_program_free function:
 * This function frees the allocated memory for a CloogProgram structure.
 */
void cloog_program_free(CloogProgram * program)
{ cloog_names_free(program->names) ;
  cloog_loop_free(program->loop) ;
  cloog_domain_free(program->context) ;
  cloog_block_list_free(program->blocklist) ;
  if (program->scaldims != NULL)
  free(program->scaldims) ;
  
  free(program) ;
}


/******************************************************************************
 *                               Reading function                             *
 ******************************************************************************/


/**
 * cloog_scattering_list_read
 * Read in a list of scattering functions for the nb_statements
 * domains in loop.
 */
static CloogScatteringList *cloog_scattering_list_read(FILE * foo,
	CloogLoop *loop, int nb_statements, int nb_parameters)
{
    int nb_scat = 0;
    char s[MAX_STRING];
    CloogScatteringList *list = NULL, **next = &list;

    /* We read first the number of scattering functions in the list. */
    do {
	if (!fgets(s, MAX_STRING, foo))
	    break;
    } while ((*s=='#' || *s=='\n') || (sscanf(s, " %d", &nb_scat) < 1));

    if (nb_scat == 0)
	return NULL;

    if (nb_scat > nb_statements)
	cloog_die("too many scattering functions.\n");

    while (nb_scat--) {
	*next = (CloogScatteringList *)malloc(sizeof(CloogScatteringList));
	(*next)->scatt = cloog_domain_read_scattering(loop->domain, foo);
	(*next)->next = NULL;

	next = &(*next)->next;
	loop = loop->next;
    }
    return list;
}


static void cloog_program_construct_block_list(CloogProgram *p)
{
    CloogLoop *loop;
    CloogBlockList **next = &p->blocklist;

    for (loop = p->loop; loop; loop = loop->next) {
	*next = cloog_block_list_alloc(loop->block);
	next = &(*next)->next;
    }
}


/**
 * cloog_program_read function:
 * This function read the informations to put in a CloogProgram structure from
 * a file (file, possibly stdin). It returns a pointer to a CloogProgram
 * structure containing the read informations.
 * - October 25th 2001: first version.
 * - September 9th 2002: - the big reading function is now split in several
 *                         functions (one per read data structure).
 *                       - adaptation to the new file format with naming.
 */
CloogProgram * cloog_program_read(FILE * file, CloogOptions * options)
{
  int i, nb_statements;
  char s[MAX_STRING], language, prefix[2]={'c','\0'};
  CloogLoop * current, * next ;
  CloogScatteringList * scatteringl;
  CloogNames *n;
  CloogProgram * p ;
      
  /* Memory allocation for the CloogProgram structure. */
  p = cloog_program_malloc() ;
  
  /* First of all, we read the language to use. */
  while (fgets(s,MAX_STRING,file) == 0) ;
  while ((*s=='#'||*s=='\n') || (sscanf(s," %c",&language)<1))
  fgets(s,MAX_STRING,file) ;
  p->language = language ;
    
  p->names = n = cloog_names_alloc();

  /* We then read the context data. */
  p->context = cloog_domain_read_context(options->state, file);
  n->nb_parameters = cloog_domain_parameter_dimension(p->context);
  
  /* First part of the CloogNames structure: reading of the parameter names. */
  n->parameters = cloog_names_read_strings(file, n->nb_parameters,
						NULL, FIRST_PARAMETER);
      
  /* We read the statement number. */
  while (fgets(s,MAX_STRING,file) == 0) ;
  while ((*s=='#'||*s=='\n') || (sscanf(s," %d",&nb_statements)<1))
  fgets(s,MAX_STRING,file) ;

  /* Statements and domains reading for each statement. */
  if (nb_statements > 0)
  { /* Reading of the first domain. */
    p->loop = cloog_loop_read(options->state, file, 0, n->nb_parameters);
    
    if (p->loop->domain != NULL)
      n->nb_iterators = cloog_domain_dimension(p->loop->domain);
    else
      n->nb_iterators = 0;
    
    /* And the same for each next domain. */
    current = p->loop ;
    for (i=2;i<=nb_statements;i++) {
      next = cloog_loop_read(options->state, file, i-1, n->nb_parameters);
      if (next->domain != NULL &&
          cloog_domain_dimension(next->domain) > n->nb_iterators)
        n->nb_iterators = cloog_domain_dimension(next->domain);
    
      current->next = next ;
      current = current->next ;
    }     
        
    /* Reading of the iterator names. */
    n->iterators = cloog_names_read_strings(file, n->nb_iterators,
						NULL, FIRST_ITERATOR);

    /* Reading and putting the scattering data in program structure. */
    scatteringl = cloog_scattering_list_read(file, p->loop, nb_statements,
						n->nb_parameters);
    
    if (scatteringl != NULL)
    { if (cloog_scattering_list_lazy_same(scatteringl))
        cloog_msg(options, CLOOG_WARNING,
			"some scattering functions are similar.\n");
      
      p->nb_scattdims = cloog_scattering_dimension(scatteringl->scatt,
							    p->loop->domain);
      n->nb_scattering = p->nb_scattdims;
      n->scattering = cloog_names_read_strings(file, p->nb_scattdims, prefix, -1);
    
      /* The boolean array for scalar dimensions is created and set to 0. */
      p->scaldims = (int *)malloc(p->nb_scattdims*(sizeof(int))) ;
      if (p->scaldims == NULL) 
	cloog_die("memory overflow.\n");
      for (i=0;i<p->nb_scattdims;i++)
      p->scaldims[i] = 0 ;
      
      /* We try to find blocks in the input problem to reduce complexity. */
      if (!options->noblocks)
	cloog_program_block(p, scatteringl, options);
      if (!options->noscalars)
	cloog_program_extract_scalars(p, scatteringl, options);
      
      cloog_program_scatter(p, scatteringl, options);
      cloog_scattering_list_free(scatteringl);

      if (!options->noblocks)
	p->loop = cloog_loop_block(p->loop, p->scaldims, p->nb_scattdims);
    }
    else
    { p->nb_scattdims = 0 ;
      p->scaldims  = NULL ;
    }
  
    cloog_names_scalarize(p->names,p->nb_scattdims,p->scaldims) ;

    cloog_program_construct_block_list(p);
  }
  else
  { p->loop      = NULL ;
    p->blocklist = NULL ;
    p->scaldims  = NULL ;
  }
   
  return(p) ;
}


/******************************************************************************
 *                            Processing functions                            *
 ******************************************************************************/
 
 
/**
 * cloog_program_malloc function:
 * This function allocates the memory space for a CloogProgram structure and
 * sets its fields with default values. Then it returns a pointer to the
 * allocated space.
 * - November 21th 2005: first version.
 */
CloogProgram * cloog_program_malloc()
{ CloogProgram * program ;
  
  /* Memory allocation for the CloogProgram structure. */
  program = (CloogProgram *)malloc(sizeof(CloogProgram)) ;
  if (program == NULL) 
    cloog_die("memory overflow.\n");
  
  /* We set the various fields with default values. */
  program->language     = 'c' ;
  program->nb_scattdims = 0 ;
  program->context      = NULL ;
  program->loop         = NULL ;
  program->names        = NULL ;
  program->blocklist    = NULL ;
  program->scaldims     = NULL ;
  program->usr          = NULL;
  
  return program ;
}  


/**
 * cloog_program_generate function:
 * This function calls the Quillere algorithm for loop scanning. (see the
 * Quillere paper) and calls the loop simplification function.
 * - depth is the loop depth we want to optimize (guard free as possible),
 *   the first loop depth is 1 and anegative value is the infinity depth.
 * - sep_level is the level number where we want to start loop separation.
 **
 * - October 26th 2001: first version. 
 * - April   19th 2005: some basic fixes and memory usage feature.
 * - April   29th 2005: (bug fix, bug found by DaeGon Kim) see case 2 below.
 */ 
CloogProgram * cloog_program_generate(program, options)
CloogProgram * program ;
CloogOptions * options ;
{
#ifdef CLOOG_RUSAGE
  float time;
  struct rusage start, end ;
#endif
  CloogLoop * loop ;
#ifdef CLOOG_MEMORY
  char status_path[MAX_STRING_VAL] ;
  FILE * status ;
 
  /* We initialize the memory need to 0. */
  options->memory = 0 ;
#endif

  if (options->override)
  {
    cloog_msg(options, CLOOG_WARNING,
    "you are using -override option, be aware that the "
    "generated\n                code may be incorrect.\n") ;
  }
  else
  { /* Playing with options may be dangerous, here are two possible issues :
     * 1. Using -l option less than scattering dimension number may lead to
     *    an illegal target code (since the scattering is not respected), if
     *    it is the case, we set -l depth to the first acceptable value.
     */
    if ((program->nb_scattdims > options->l) && (options->l >= 0))
    {
      cloog_msg(options, CLOOG_WARNING,
      "-l depth is less than the scattering dimension number "
      "(the \n                generated code may be incorrect), it has been "
      "automaticaly set\n                to this value (use option -override "
      "to override).\n") ;
      options->l = program->nb_scattdims ;
    }
      
    /* 2. Using -f option greater than one while -l depth is greater than the
     *    scattering dimension number may lead to iteration duplication (try
     *    test/daegon_lu_osp.cloog with '-f 3' to test) because of the step 4b
     *    of the cloog_loop_generate function, if it is the case, we set -l to
     *    the first acceptable value.
     */
    if (((options->f > 1) || (options->f < 0)) &&
        ((options->l > program->nb_scattdims) || (options->l < 0)))
    {
      cloog_msg(options, CLOOG_WARNING,
      "-f depth is more than one, -l depth has been "
      "automaticaly set\n                to the scattering dimension number "
      "(target code may have\n                duplicated iterations), -l depth "
      "has been automaticaly set to\n                this value (use option "
      "-override to override).\n") ;
      options->l = program->nb_scattdims ;
    }
  }
  
#ifdef CLOOG_RUSAGE
  getrusage(RUSAGE_SELF, &start) ;
#endif
  if (program->loop != NULL)
  { loop = program->loop ;
    
    /* Here we go ! */
    loop = cloog_loop_generate(loop, program->context, 0, 0,
                               program->scaldims,
			       program->nb_scattdims,
			       program->names->nb_parameters,
			       options);
			          
#ifdef CLOOG_MEMORY
    /* We read into the status file of the process how many memory it uses. */
    sprintf(status_path,"/proc/%d/status",getpid()) ;
    status = fopen(status_path, "r") ;
    while (fscanf(status,"%s",status_path) && strcmp(status_path,"VmData:")!=0);
    fscanf(status,"%d",&(options->memory)) ;
    fclose(status) ;
#endif
    
    if ((!options->nosimplify) && (program->loop != NULL))
    loop = cloog_loop_simplify(loop,program->context,1,
                               program->names->nb_parameters);
   
    program->loop = loop ;
  }
    
#ifdef CLOOG_RUSAGE
  getrusage(RUSAGE_SELF, &end) ;
  /* We calculate the time spent in code generation. */
  time =  (end.ru_utime.tv_usec -  start.ru_utime.tv_usec)/(float)(MEGA) ;
  time += (float)(end.ru_utime.tv_sec - start.ru_utime.tv_sec) ;
  options->time = time ;
#endif
  
  return program ;
}


/**
 * cloog_program_block function:
 * this function gives a last chance to the lazy user to consider statement
 * blocks instead of some statement lists where the whole list may be
 * considered as a single statement from a code generation point of view.
 * For instance two statements with the same iteration domain and the same
 * scattering functions may be considered as a block. This function is lazy
 * and can only find very simple forms of trivial blocks (see
 * cloog_domain_lazy_block function for more details). The useless loops and
 * scattering functions are removed and freed while the statement list of
 * according blocks are filled.
 * - program is the whole program structure (befaore applying scattering),
 * - scattering is the list of scattering functions.
 **
 * - April   30th 2005: first attempt.
 * - June 10-11th 2005: first working version.
 */
void cloog_program_block(CloogProgram *program,
	CloogScatteringList *scattering, CloogOptions *options)
{ int blocked_reference=0, blocked=0, nb_blocked=0 ;
  CloogLoop * reference, * start, * loop ;
  CloogScatteringList * scatt_reference, * scatt_loop, * scatt_start;
  
  if ((program->loop == NULL) || (program->loop->next == NULL))
  return ;
  
  /* The process will use three variables for the linked list :
   * - 'start' is the starting point of a new block,
   * - 'reference' is the node of the block used for the block checking,
   * - 'loop' is the candidate to be inserted inside the block.
   * At the beginning of the process, the linked lists are as follow:
   *         O------>O------>O------>O------>NULL
   *         |       |
   *       start    loop
   *     reference
   */

  reference       = program->loop ;
  start           = program->loop ;
  loop            = reference->next ;
  scatt_reference = scattering ;
  scatt_start     = scattering ;
  scatt_loop      = scattering->next ;
   
  while (loop != NULL)
  { if (cloog_domain_lazy_equal(reference->domain,loop->domain) &&
        cloog_scattering_lazy_block(scatt_reference->scatt, scatt_loop->scatt,
	                        scattering,program->nb_scattdims))
    { /* If we find a block we update the links:
       *     +---------------+
       *     |               v
       *     O       O------>O------>O------>NULL
       *     |       |
       *   start    loop
       * reference
       */
      blocked = 1 ;
      nb_blocked ++ ;
      cloog_block_merge(start->block,loop->block); /* merge frees loop->block */
      loop->block = NULL ;
      start->next = loop->next ;
      scatt_start->next = scatt_loop->next ;
    }
    else
    { /* If we didn't find a block, the next start of a block is updated:
       *     O------>O------>O------>O------>NULL
       *     |       |
       * reference start
       *           loop
       */
      blocked= 0 ;
      start = loop ;
      scatt_start = scatt_loop ;
    }

    /* If the reference node has been included into a block, we can free it. */
    if (blocked_reference)
    { reference->next = NULL ;
      cloog_loop_free(reference) ;
      cloog_scattering_free(scatt_reference->scatt);
      free(scatt_reference) ;
    }
    
    /* The reference and the loop are now updated for the next try, the
     * starting position depends on the previous step.
     *       O   ?   O------>O------>O------>NULL
     *               |       |
     *           reference loop
     */
    reference       = loop ;
    loop            = loop->next ;
    scatt_reference = scatt_loop ;
    scatt_loop      = scatt_loop->next ;
    
    /* We mark the new reference as being blocked or not, if will be freed
     * during the next while loop execution.
     */
    if (blocked)
    blocked_reference = 1 ;
    else
    blocked_reference = 0 ;
  }
  
  /* We free the last blocked reference if any (since in the while loop it was
   * freed during the next loop execution, it was not possible to free the
   * last one inside).
   */
  if (blocked_reference)
  { reference->next = NULL ;
    cloog_loop_free(reference) ;
    cloog_scattering_free(scatt_reference->scatt);
    free(scatt_reference) ;
  }
  
  if (nb_blocked != 0)
    cloog_msg(options, CLOOG_INFO, "%d domains have been blocked.\n", nb_blocked);
}


/**
 * cloog_program_extract_scalars function:
 * this functions finds and removes the dimensions of the scattering functions
 * when they are scalar (i.e. of the shape "dim + scalar = 0") for all
 * scattering functions. The reason is that the processing of such dimensions
 * is trivial and do not need neither a row and a column in the matrix
 * representation of the domain (this will save memory) neither the full
 * Quillere processing (this will save time). The scalar dimensions data are
 * dispatched in the CloogProgram structure (the boolean vector scaldims will
 * say which original dimensions are scalar or not) and to the CloogBlock
 * structures (each one has a scaldims vector that contains the scalar values).
 * - June 14th 2005: first developments.
 * - June 30th 2005: first version.
 */ 
void cloog_program_extract_scalars(CloogProgram *program,
	CloogScatteringList *scattering, CloogOptions *options)
{ int i, j, scalar, current, nb_scaldims=0 ;
  CloogScatteringList *start;
  CloogScattering *old;
  CloogLoop *loop;
  CloogBlock * block ;

  start = scattering ;
    
  for (i=0;i<program->nb_scattdims;i++)
  { scalar = 1 ;
    scattering = start ;
    while (scattering != NULL)
    { if (!cloog_scattering_lazy_isscalar(scattering->scatt, i, NULL))
      { scalar = 0 ;
        break ;
      }
      scattering = scattering->next ;
    }
    
    if (scalar)
    { nb_scaldims ++ ;
      program->scaldims[i] = 1 ;
    }
  }
  
  /* If there are no scalar dimensions, we can continue directly. */
  if (!nb_scaldims)
  return ;

  /* Otherwise, in each block, we have to put the number of scalar dimensions,
   * and to allocate the memory for the scalar values.
   */
  for (loop = program->loop; loop; loop = loop->next) {
    block = loop->block;
    block->nb_scaldims = nb_scaldims ;
    block->scaldims = (cloog_int_t *)malloc(nb_scaldims*sizeof(cloog_int_t));
    for (i=0;i<nb_scaldims;i++)
    cloog_int_init(block->scaldims[i]);
  }
  
  /* Then we have to fill these scalar values, so we can erase those dimensions
   * from the scattering functions. It's easier to begin with the last one,
   * since there would be an offset otherwise (if we remove the i^th dimension,
   * then the next one is not the (i+1)^th but still the i^th...).
   */
  current = nb_scaldims - 1 ;
  for (i=program->nb_scattdims-1;i>=0;i--)
  if (program->scaldims[i])
  {
    scattering = start ;
    for (loop = program->loop; loop; loop = loop->next) {
      block = loop->block;
      if (!cloog_scattering_lazy_isscalar(scattering->scatt, i,
						&block->scaldims[current])) {
	/* We should have found a scalar value: if not, there is an error. */
	cloog_die("dimension %d is not scalar as expected.\n", i);
      }
      scattering = scattering->next ;
    } 
  
    scattering = start ;
    while (scattering != NULL) {
      old = scattering->scatt;
      scattering->scatt = cloog_scattering_erase_dimension(old, i);
      cloog_scattering_free(old);
      scattering = scattering->next ;
    }
    current-- ;
  }
  
  /* We postprocess the scaldims array in such a way that each entry is how
   * many scalar dimensions follows + 1 (the current one). This will make 
   * some other processing easier (e.g. knowledge of some offsets).
   */
  for (i=0;i<program->nb_scattdims-1;i++)
  { if (program->scaldims[i])
    { j = i + 1 ;
      while ((j < program->nb_scattdims) && program->scaldims[j])
      { program->scaldims[i] ++ ;
        j ++ ;
      }
    }
  }
  
  if (nb_scaldims != 0)
    cloog_msg(options, CLOOG_INFO, "%d dimensions (over %d) are scalar.\n",
          nb_scaldims,program->nb_scattdims) ;
}


/**
 * cloog_program_scatter function:
 * This function adds the scattering (scheduling) informations in a program.
 * If names is NULL, this function create names itself such that the i^th
 * name is ci.
 * - November 6th 2001: first version. 
 */
void cloog_program_scatter(CloogProgram *program,
			CloogScatteringList *scattering, CloogOptions *options)
{ int scattering_dim, scattering_dim2, not_enough_constraints=0 ;
  CloogLoop * loop ;
  
  if ((program != NULL) && (scattering != NULL))
  { loop = program->loop ;
    
    /* We compute the scattering dimension and check it is >=0. */
    scattering_dim = cloog_scattering_dimension(scattering->scatt, loop->domain);
    if (scattering_dim < 0)
      cloog_die("scattering has not enough dimensions.\n");
    if (!cloog_scattering_fully_specified(scattering->scatt, loop->domain))
    not_enough_constraints ++ ;
         
    /* The scattering dimension may have been modified by scalar extraction. */
    scattering_dim = cloog_scattering_dimension(scattering->scatt, loop->domain);

    /* Finally we scatter all loops. */
    cloog_loop_scatter(loop, scattering->scatt);
    loop = loop->next ;
    scattering = scattering->next ;    
    
    while ((loop != NULL) && (scattering != NULL))
    { scattering_dim2 = cloog_scattering_dimension(scattering->scatt,
								loop->domain);
      if (scattering_dim2 != scattering_dim)
        cloog_die("scattering dimensions are not the same.\n") ;
      if (!cloog_scattering_fully_specified(scattering->scatt, loop->domain))
      not_enough_constraints ++ ;
      
      cloog_loop_scatter(loop, scattering->scatt);
      loop = loop->next ;
      scattering = scattering->next ;
    }
    if ((loop != NULL) || (scattering != NULL))
      cloog_msg(options, CLOOG_WARNING,
                    "there is not a scattering for each statement.\n");
    
    if (not_enough_constraints)
      cloog_msg(options, CLOOG_WARNING, "not enough constraints for "
                    "%d scattering function(s).\n",not_enough_constraints) ;
  }
}
