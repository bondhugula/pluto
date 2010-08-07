/*
  $Id: errors.c,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $

  Exception management.
  See "arithmetic_errors.h".

  $Log: errors.c,v $
  Revision 1.1.1.1  2010/01/25 07:06:34  uday
  initial import

  Revision 1.1.1.1  2008/05/05 12:39:54  bondhugu
  Import of polylib v. 5.22.3

  Revision 1.6  2006/03/15 19:59:37  verdoolaege
  arith: add some missing consts

  Revision 1.5  2004/08/31 18:01:56  verdoolaege
  remove warning

  Revision 1.4  2004/02/11 10:19:54  verdoolaege
  add const qualifier

  Revision 1.3  2004/02/08 21:53:27  kienhuis
  Update from Fabien Coelho, via Bart Kienhuis

  I've updated here in the C3/Linear library the arithmetic_error
  package that I developped (with others) to handle exceptions in C.
  It adds a simple callback feature which is needed for pips here.
  If you do not use it, it should not harm;-)

  Revision 1.27  2003/09/04 09:40:37  coelho
  init added.
  verbosity mask added.

  Revision 1.26  2003/09/03 14:05:20  coelho
  push/pop callbacks called.

  Revision 1.20  2003/08/18 09:55:09  coelho
  get_exception_name added...

  Revision 1.19  2003/06/13 13:59:07  coelho
  const out.

  Revision 1.18  2003/06/13 13:54:47  coelho
  hop.

  Revision 1.17  2002/04/02 08:44:54  coelho
  timeout_error ajoute.

  Revision 1.16  2000/10/27 13:26:03  ancourt
  exception_thrown -> linear_number_of_exception_thrown

  Revision 1.15  2000/07/27 15:21:55  coelho
  message++

  Revision 1.14  2000/07/27 14:59:59  coelho
  trace added.

  Revision 1.13  2000/07/26 08:41:23  coelho
  the_last_just_thrown_exception management added.

  Revision 1.12  1998/10/26 18:48:34  coelho
  message++.

  Revision 1.11  1998/10/26 14:38:06  coelho
  constants back in.

  Revision 1.10  1998/10/26 14:35:47  coelho
  constants in .h.

  Revision 1.9  1998/10/24 15:36:22  coelho
  better exception error messages...

  Revision 1.8  1998/10/24 15:19:17  coelho
  more verbose throw...

  Revision 1.7  1998/10/24 14:31:13  coelho
  new stack implementation.
  checks for matching catch/uncatch.
  debug mode.

  Revision 1.6  1998/10/24 09:31:06  coelho
  RCS headers.
  'const' tried.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetique.h"

/* global constants to designate exceptions.
   to be put in the type field bellow.
   cproto 4.6 does not line 'const'...
*/
unsigned int overflow_error = 1;
unsigned int simplex_arithmetic_error = 2;
unsigned int user_exception_error = 4;
unsigned int parser_exception_error = 8;
unsigned int timeout_error = 16;

/* catch all */
unsigned int any_exception_error = ~0;


/*
 * On a few systems, type boolean and/or its values FALSE, TRUE may appear
 * in standard header files.  Or you may have conflicts with application-
 * specific header files that you want to include together with these files.
 * Defining HAVE_BOOLEAN before including jpeglib.h should make it work.
 */

#ifndef HAVE_BOOLEAN
typedef int boolean;
#endif
#ifndef FALSE                   /* in case these macros already exist */
#define FALSE   0               /* values of boolean */
#endif
#ifndef TRUE
#define TRUE    1
#endif

char * get_exception_name(unsigned int exception)
{
  if (exception==overflow_error)
    return "overflow_error exception";
  if (exception==simplex_arithmetic_error)
    return "simplex_arithmetic_error exception";
  if (exception==user_exception_error)
    return "user_exception_error exception";
  if (exception==parser_exception_error)
    return "parser_exception_error exception";
  if (exception==timeout_error)
    return "timeout_error exception";
  if (exception==any_exception_error)
    return "all exceptions mask";

  return "unknown or mixed exception";
}

/* keep track of last thrown exception for RETHROW()
 */
unsigned int the_last_just_thrown_exception = 0;

/* whether to run in debug mode (that is to trace catch/uncatch/throw)
 */
static int linear_exception_debug_mode = FALSE;
static unsigned int linear_exception_verbose = 1 | 2 | 16 ;

/* A structure for the exception stack.
 */
typedef struct 
{
  /* exception type.
   */
  int what;

  /* where to jump to.
   */
  jmp_buf where;

  /* location of the CATCH to be matched against the UNCATCH.
   */
  const char * function;
  const char * file;
  int    line;
} 
  linear_exception_holder;

/* exception stack.
   maximum extension.
   current index (next available bucket)
 */
#define MAX_STACKED_CONTEXTS 64
static linear_exception_holder exception_stack[MAX_STACKED_CONTEXTS];
static int exception_index = 0;

/* callbacks...
 */
static exception_callback_t push_callback = NULL;
static exception_callback_t pop_callback = NULL;

void set_exception_callbacks(exception_callback_t push, 
			     exception_callback_t pop)
{
  if (push_callback!=NULL || pop_callback!=NULL)
  {
    fprintf(stderr, "exception callbacks already defined! (%p, %p)\n",
	    push_callback, pop_callback);
    abort();
  }

  push_callback = push;
  pop_callback = pop;
}

/* total number of exceptions thrown, for statistics.
 */
int linear_number_of_exception_thrown = 0;

/* dump stack
 */
void dump_exception_stack_to_file(FILE * f)
{
  int i;
  fprintf(f, "[dump_exception_stack_to_file] size=%d\n", exception_index);
  for (i=0; i<exception_index; i++)
  {
    fprintf(f, 
	    "%d: [%s:%d in %s (%d)]\n",
	    i, 
	    exception_stack[i].file,
	    exception_stack[i].line,
	    exception_stack[i].function,
	    exception_stack[i].what);
  }
  fprintf(f, "\n");
}

void dump_exception_stack()
{
  dump_exception_stack_to_file(stderr);
}

#define exception_debug_message(type) 				  \
    fprintf(stderr, "%s[%s:%d %s (%d)/%d]\n", 			  \
	    type, file, line, function, what, exception_index) 

#define exception_debug_trace(type) \
    if (linear_exception_debug_mode) { exception_debug_message(type); }


/* push a what exception on stack.
 */
jmp_buf * 
push_exception_on_stack(
    int what,
    const char * function,
    const char * file,
    int line)
{
  exception_debug_trace("PUSH ");

  if (exception_index==MAX_STACKED_CONTEXTS)
  {
    exception_debug_message("push");
    fprintf(stderr, "exception stack overflow\n");
    dump_exception_stack();
    abort();
  }

  if (push_callback) push_callback(file, function, line);

  the_last_just_thrown_exception = 0;

  exception_stack[exception_index].what = what;
  exception_stack[exception_index].function = function;
  exception_stack[exception_index].file = file;
  exception_stack[exception_index].line = line;

  return & exception_stack[exception_index++].where;
}

#if !defined(same_string_p)
#define same_string_p(s1, s2) (strcmp((s1),(s2))==0)
#endif

/* pop a what exception.
   check for any mismatch!
 */
void
pop_exception_from_stack(
    int what,
    const char * function,
    const char * file,
    int line)
{  
  exception_debug_trace("POP  ");

  if (exception_index==0)
  {
    exception_debug_message("pop");
    fprintf(stderr, "exception stack underflow\n");
    dump_exception_stack();
    abort();
  }

  if (pop_callback) pop_callback(file, function, line);

  exception_index--;
  the_last_just_thrown_exception = 0;

  if ((exception_stack[exception_index].what != what) ||
      !same_string_p(exception_stack[exception_index].file, file) ||
      !same_string_p(exception_stack[exception_index].function, function))
  {
    exception_debug_message("pop");
    fprintf(stderr, 
	    "exception stack mismatch at depth=%d:\n"
	    "   CATCH: %s:%d in %s (%d)\n"
	    " UNCATCH: %s:%d in %s (%d)\n",
	    exception_index,
	    exception_stack[exception_index].file,
	    exception_stack[exception_index].line,
	    exception_stack[exception_index].function,
	    exception_stack[exception_index].what,
	    file, line, function, what);
    dump_exception_stack();
    abort();
  }
}

/* throws an exception of a given type by searching for 
   the specified 'what' in the current exception stack.
*/
void throw_exception(
    int what,
    const char * function,
    const char * file,
    int line)
{
  int i;
  
  exception_debug_trace("THROW");

  the_last_just_thrown_exception = what; /* for rethrow */

  for (i=exception_index-1; i>=0; i--)
  {
    if (pop_callback) 
      /* pop with push parameters! */
      pop_callback(exception_stack[i].file,
		   exception_stack[i].function,
		   exception_stack[i].line);

    if (exception_stack[i].what & what) 
    {
      exception_index = i;
      linear_number_of_exception_thrown++;

      if (linear_exception_debug_mode)
	fprintf(stderr, "---->[%s:%d %s (%d)/%d]\n", 
		exception_stack[i].file,
		exception_stack[i].line,
		exception_stack[i].function,
		exception_stack[i].what,
		i);
  
      /* trace some exceptions... */
      if (linear_exception_verbose & what)
	fprintf(stderr, "exception %d/%d: %s(%s:%d) -> %s(%s:%d)\n",
		what, exception_stack[i].what,
		function, file, line,
		exception_stack[i].function, 
		exception_stack[i].file,
		exception_stack[i].line);

      longjmp(exception_stack[i].where, 0);
    }
  }

  /* error. */
  exception_debug_message("throw");
  fprintf(stderr,
	  "exception not found in stack:\n"
	  "an exception was THROWN without a proper matching CATCH\n");
  dump_exception_stack();
  abort();
}

void linear_initialize_exception_stack(
  unsigned int verbose_exceptions,
  exception_callback_t push, 
  exception_callback_t pop)
{
  linear_exception_verbose = verbose_exceptions;
  set_exception_callbacks(push, pop);
}
