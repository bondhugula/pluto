/* 
 * $Id: arithmetic_errors.h,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 *
 * managing arithmetic errors...
 * detecting and managing arithmetic errors on Values should be
 * systematic. These macros gives a C++ look and feel to this
 * management. 
 *
 * (c) CA et FC, Sept 1997
 *
 * $Log: arithmetic_errors.h,v $
 * Revision 1.1.1.1  2010/01/25 07:06:34  uday
 * initial import
 *
 * Revision 1.1.1.1  2008/05/05 12:39:54  bondhugu
 * Import of polylib v. 5.22.3
 *
 * Revision 1.4  2006/03/15 19:59:37  verdoolaege
 * arith: add some missing consts
 *
 * Revision 1.3  2004/02/08 21:53:27  kienhuis
 * Update from Fabien Coelho, via Bart Kienhuis
 *
 * I've updated here in the C3/Linear library the arithmetic_error
 * package that I developped (with others) to handle exceptions in C.
 * It adds a simple callback feature which is needed for pips here.
 * If you do not use it, it should not harm;-)
 *
 * Revision 1.34  2003/09/03 13:59:46  coelho
 * ++
 *
 * Revision 1.33  2003/09/03 13:35:34  coelho
 * no more callback.
 *
 * Revision 1.32  2003/08/18 14:55:38  coelho
 * callback fix.
 *
 * Revision 1.31  2003/08/18 14:16:45  coelho
 * NULL callback added.
 *
 * Revision 1.30  2003/06/13 13:59:55  coelho
 * hop.
 *
 * Revision 1.29  2000/07/27 15:01:55  coelho
 * hop.
 *
 * Revision 1.28  2000/07/26 09:11:58  coelho
 * hop.
 *
 * Revision 1.27  2000/07/26 09:07:32  coelho
 * *** empty log message ***
 *
 * Revision 1.26  2000/07/26 09:06:32  coelho
 * the_last_just_thrown_exception declared.
 *
 * Revision 1.25  2000/07/26 08:41:40  coelho
 * RETHROW added.
 *
 * Revision 1.24  1998/10/26 14:37:48  coelho
 * constants moved out.
 *
 * Revision 1.23  1998/10/26 14:36:13  coelho
 * constants explicitely defined in .h.
 *
 * Revision 1.22  1998/10/24 15:18:26  coelho
 * THROW macro updated to tell its source.
 *
 * Revision 1.21  1998/10/24 14:33:08  coelho
 * parser exception added.
 *
 * Revision 1.20  1998/10/24 14:32:45  coelho
 * simpler macros.
 *
 * Revision 1.19  1998/10/24 09:22:47  coelho
 * size update.
 *
 * Revision 1.18  1998/10/24 09:21:45  coelho
 * const added to constants.
 *
 */

#if !defined(linear_arithmetic_error_included)
#define linear_arithmetic_error_included

#include <setjmp.h>

typedef void (*exception_callback_t)(const char *, const char *, int);

/*
const unsigned int overflow_error = 1;
const unsigned int simplex_arithmetic_error = 2;
const unsigned int user_exception_error = 4;
const unsigned int parser_exception_error = 8;
const unsigned int any_exception_error = ~0;
*/

/* use gnu cpp '__FUNCTION__' extension if possible.
 */
#if defined(__GNUC__)
#define __CURRENT_FUNCTION_NAME__ __FUNCTION__
#else
#define __CURRENT_FUNCTION_NAME__ "<unknown>"
#endif

/* 'const' out because of cproto 4.6. FC 13/06/2003 */
#define EXCEPTION extern unsigned int

#define THROW(what) \
   (throw_exception(what, __CURRENT_FUNCTION_NAME__, __FILE__, __LINE__))

#define CATCH(what) 							\
   if (setjmp(*push_exception_on_stack(what, __CURRENT_FUNCTION_NAME__,	\
				     __FILE__, __LINE__)))

#define UNCATCH(what)						\
     (pop_exception_from_stack(what, __CURRENT_FUNCTION_NAME__,	\
			       __FILE__, __LINE__))

#define TRY else

extern unsigned int the_last_just_thrown_exception;
#define RETHROW() THROW(the_last_just_thrown_exception)

#endif /* linear_arithmetic_error_included */

/* end of it.
 */
