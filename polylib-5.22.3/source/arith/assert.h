/* Version "abort" de l'assert de /usr/include/assert.h 
 * Il est installe dans Linear de maniere a masquer /usr/include/assert.h
 *
 * Il faut faire un include de <stdio.h> pour l'utiliser.
 */

# ifndef NDEBUG
# define _assert(ex)	{if (!(ex)){(void)fprintf(stderr,"Assertion failed: file \"%s\", line %d\n", __FILE__, __LINE__);(void) abort();}}
# define assert(ex)	_assert(ex)
# else
# define _assert(ex)
# define assert(ex)
# endif
