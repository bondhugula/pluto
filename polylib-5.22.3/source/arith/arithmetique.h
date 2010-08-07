/* header file built by cproto */
#ifndef arithmetique_header_included
#define arithmetique_header_included

/** package arithmetique
 *
 * $Id: arithmetique.h,v 1.1.1.1 2010/01/25 07:06:34 uday Exp $
 *
 * Francois Irigoin, mai 1989
 *
 * Modifications
 *  - rewrite of DIVIDE which was wrong (Remi Triolet, Francois Irigoin, 
 *    april 90)
 *  - simplification of POSITIVE_DIVIDE by suppressing one modulo
 *  - B.Meister : added addmul, operation existing in gmp and quite useful 
 *    (05-2005)
 */

/* We would like linear to be generic about the "integer" type used
 * to represent integer values. Thus Value is defined here. It should
 * be changed to "int" "long" or "long long". In an ideal world,
 * any source modification should be limited to this package.
 *
 * Indeed, we cannot switch easily to bignums that need constructors 
 * dans destructors... That would lead to too many modifications...
 * C++ would make things easier and cleaner...
 *
 * Fabien COELHO
 */

#include <stdio.h>
#include <limits.h>   /* Included for getting constants: INT_MAX, etc.. */

#ifdef GNUMP
#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#ifndef mp_get_memory_functions
#if defined(__cplusplus)
extern "C" {
#endif
void mp_get_memory_functions(
		void *(**alloc_func_ptr) (size_t),
		void *(**realloc_func_ptr) (void *, size_t, size_t),
		void (**free_func_ptr) (void *, size_t));
#if defined(__cplusplus)
}
#endif
#endif
#endif 

#ifdef CLN
#include <sstream>
#define WANT_OBFUSCATING_OPERATORS
#include <cln/cln.h>
#endif

/* 
   #        ####   #    #   ####           #        ####   #    #   ####
   #       #    #  ##   #  #    #          #       #    #  ##   #  #    #
   #       #    #  # #  #  #               #       #    #  # #  #  #
   #       #    #  #  # #  #  ###          #       #    #  #  # #  #  ###
   #       #    #  #   ##  #    #          #       #    #  #   ##  #    #
   ######   ####   #    #   ####           ######   ####   #    #   ####
   
*/

/* 
 * Constants like LONG_LONG_MAX are not defined with ansi options, so they are
 * defined here. 
 */  

#ifndef LONG_LONG_MAX

/* would fix on solaris:
 * #define LONG_LONG_MAX LLONG_MAX
 * #define LONG_LONG_MIN LLONG_MIN
 */

#ifndef __LONG_LONG_MAX__
#define __LONG_LONG_MAX__ 9223372036854775807LL
#endif
#undef LONG_LONG_MAX
#define LONG_LONG_MAX __LONG_LONG_MAX__
#undef LONG_LONG_MIN
#define LONG_LONG_MIN (-LONG_LONG_MAX-1)
#undef ULONG_LONG_MAX
#define ULONG_LONG_MAX (LONG_LONG_MAX * 2ULL + 1)
#endif

#if defined(LINEAR_VALUE_IS_LONGLONG)

#define LINEAR_VALUE_STRING "long long int"
typedef long long int Value;
#if defined(WIN32) && !defined(unix)
    /* Mingw or Windows need an incompatible format string. */
#   define VALUE_FMT "%I64d"
#else
#   define VALUE_FMT "%lld"
#endif
#define VALUE_CONST(val) (val##LL) 

/* 
 * CAUTION! 'VALUE_MIN' is defined as 'LONG_LONG_MIN +1' so as to preserve the
 * symmetry (-min==max) and to have a NAN value. FC 
 */ 

#define VALUE_NAN LONG_LONG_MIN
#define VALUE_MIN (LONG_LONG_MIN+1LL)
#define VALUE_MAX LONG_LONG_MAX
#define VALUE_SQRT_MIN long_to_value(LONG_MIN) 
#define VALUE_SQRT_MAX long_to_value(LONG_MAX)
#define VALUE_ZERO (0LL)
#define VALUE_ONE  (1LL)
#define VALUE_MONE (-1LL)

#define VALUE_TO_LONG(val) \
    ((long)((val)>(Value)LONG_MIN&&(val)<=(Value)LONG_MAX)?\
     (val):(THROW(overflow_error), LONG_MIN))

#define VALUE_TO_INT(val) \
    ((int)((val)>(Value)INT_MIN&&(val)<=(Value)INT_MAX)?\
     (val):(THROW(overflow_error), INT_MIN))

#define VALUE_TO_DOUBLE(val) ((double)(val))

/* #define VALUE_TO_FLOAT(val) ((float)(val)): Doesn't seem to work with gcc */
#define VALUE_TO_FLOAT(val) ((float)((int)(val)))

/* end LINEAR_VALUE_IS_LONGLONG */

/* 
 
   #        ####   #    #   ####
   #       #    #  ##   #  #    #
   #       #    #  # #  #  #
   #       #    #  #  # #  #  ###
   #       #    #  #   ##  #    #
   ######   ####   #    #   ####
 
*/

#elif defined(LINEAR_VALUE_IS_LONG)

#define LINEAR_VALUE_STRING "long int"
typedef long Value;
#define VALUE_FMT "%ld"
#define VALUE_CONST(val) (val##L)
#define VALUE_NAN LONG_MIN
#define VALUE_MIN (LONG_MIN+1L)
#define VALUE_MAX LONG_MAX
#define VALUE_SQRT_MIN int_to_value(INT_MIN)
#define VALUE_SQRT_MAX int_to_value(INT_MAX)
#define VALUE_ZERO 0L
#define VALUE_ONE  1L
#define VALUE_MONE -1L
#define VALUE_TO_LONG(val) (val)
#define VALUE_TO_INT(val) ((int)(val))
#define VALUE_TO_FLOAT(val) ((float)(val))
#define VALUE_TO_DOUBLE(val) ((double)(val))

/* end LINEAR_VALUE_IS_LONG */

/* 
   ######  #        ####     ##     #####
   #       #       #    #   #  #      #
   #####   #       #    #  #    #     #
   #       #       #    #  ######     #
   #       #       #    #  #    #     #
   #       ######   ####   #    #     #
 
*/

/*
#elif defined(LINEAR_VALUE_IS_FLOAT)

#define LINEAR_VALUE_STRING "float"
typedef float Value;
#define VALUE_FMT "%f"
#define VALUE_CONST(val) (val)
#define VALUE_MIN FLOAT_MIN
#define VALUE_MAX FLOAT_MAX
#define VALUE_ZERO 0.0
#define VALUE_ONE  1.0
#define VALUE_MONE -1.0
#define VALUE_TO_LONG(val) ((long)(val))
#define VALUE_TO_INT(val) ((int)(val))
#define VALUE_TO_FLOAT(val) ((float)(val))
#define VALUE_TO_DOUBLE(val) ((double)(val))

*/

/* end LINEAR_VALUE_IS_FLOAT */

/*
   ####   #    #    ##    #####           #   #
  #    #  #    #   #  #   #    #           # #
  #       ######  #    #  #    #         #######
  #       #    #  ######  #####            # #
  #    #  #    #  #    #  #   #           #   #
   ####   #    #  #    #  #    #
  
   */

/* Char version is used to detect invalid assignments */

#elif defined(LINEAR_VALUE_IS_CHARS)

#define LINEAR_VALUE_STRING "char"
typedef union { char *s; long l; int i; float f; double d;} Value;
#define VALUE_FMT "%s"
#define VALUE_CONST(val) ((Value)(val))
#define VALUE_NAN ((Value)(long)0xdadeebee)
#define VALUE_MIN ((Value)(long)0xdeadbeef)
#define VALUE_MAX ((Value)(long)0xfeedabee)
#define VALUE_ZERO ((Value)0)
#define VALUE_ONE  ((Value)1)
#define VALUE_MONE ((Value)-1)
#define VALUE_TO_LONG(val) (val.l)
#define VALUE_TO_INT(val) (val.i)
#define VALUE_TO_FLOAT(val) (val.f)
#define VALUE_TO_DOUBLE(val) (val.d)

/* end LINEAR_VALUE_IS_CHARS */

/*
    #    #    #   #####
    #    ##   #     #
    #    # #  #     #
    #    #  # #     #
    #    #   ##     #
    #    #    #     #

*/

#elif defined(LINEAR_VALUE_IS_INT)

#define LINEAR_VALUE_STRING "int"
typedef int Value;
#define VALUE_FMT "%d"
#define VALUE_CONST(val) (val)
#define VALUE_NAN INT_MIN
#define VALUE_MIN (INT_MIN+1)
#define VALUE_MAX INT_MAX
#define VALUE_ZERO  0
#define VALUE_ONE   1
#define VALUE_MONE -1
#define VALUE_TO_LONG(val) ((long)(val))
#define VALUE_TO_INT(val) ((int)(val))
#define VALUE_TO_FLOAT(val) ((float)(val))
#define VALUE_TO_DOUBLE(val) ((double)(val))

/* end LINEAR_VALUE_IS_INT */

#elif defined(GNUMP)

#define LINEAR_VALUE_STRING "gmp"
typedef mpz_t Value;
#define VALUE_FMT "%s"

/* don't use these, use value_set_si instead ! */
#undef VALUE_ZERO
#undef VALUE_ONE
#undef VALUE_MONE
#define VALUE_TO_LONG(val) (mpz_get_si(val))
#define VALUE_TO_INT(val) ((int)mpz_get_si(val))
#define VALUE_TO_FLOAT(val) ((float)((int)mpz_get_si(val)))
#define VALUE_TO_DOUBLE(val) (mpz_get_d(val))

#elif defined(CLN)

#define LINEAR_VALUE_STRING "cln"
typedef cln::cl_I Value;
#define VALUE_FMT "%s"

#define VALUE_TO_INT(val) (cln::cl_I_to_int(val))
#define VALUE_TO_DOUBLE(val) (cln::double_approx(val))

#endif 

/* ***************** MACROS FOR MANIPULATING VALUES ******************** */

#if defined(CLN)

#define value_init(val)        	((val).word = ((cln::cl_uint)cl_FN_tag) << cl_tag_shift)
#define value_assign(v1,v2)    	((v1) = (v2))
#define value_set_si(val,i)    	((val) = (i))    
#define value_set_double(val,d)	((val) = cln::truncate1(cln::cl_R(d)))
#define value_clear(val)       	((val) = 0)
#define value_read(val,str)    	((val) = (str))
#define value_print(Dst,fmt,val)  {std::ostringstream strm; strm << val; \
				   fprintf((Dst),(fmt),strm.str().c_str()); \
				  }
#define value_swap(v1,v2)          {Value tmp; tmp = v2; \
                                    v2 = v1; v1 = tmp;   \
                                   }

/* Boolean operators on 'Value' */

#define value_eq(v1,v2) ((v1)==(v2))
#define value_ne(v1,v2) ((v1)!=(v2))
#define value_gt(v1,v2) ((v1)>(v2))
#define value_ge(v1,v2) ((v1)>=(v2))
#define value_lt(v1,v2) ((v1)<(v2))
#define value_le(v1,v2) ((v1)<=(v2))

#define value_abs_eq(v1,v2) (cln::abs(v1)==cln::abs(v2))
#define value_abs_ne(v1,v2) (cln::abs(v1)!=cln::abs(v2))
#define value_abs_gt(v1,v2) (cln::abs(v1)>cln::abs(v2))
#define value_abs_ge(v1,v2) (cln::abs(v1)>=cln::abs(v2))
#define value_abs_lt(v1,v2) (cln::abs(v1)<cln::abs(v2))
#define value_abs_le(v1,v2) (cln::abs(v1)<=cln::abs(v2))

#define value_sign(val)      (cln::signum(val))
#define value_compare(v1,v2) (cln::compare((v1),(v2)))

#define value_addto(ref,val1,val2) 	((ref) = (val1)+(val2))
#define value_add_int(ref,val,vint)     ((ref) = (val)+(vint))
#define value_addmul(ref, val1, val2)   ((ref) += (val1)*(val2))
#define value_increment(ref,val) 	((ref) = (val)+1)
#define value_multiply(ref,val1,val2)	((ref) = (val1)*(val2))
#define value_subtract(ref,val1,val2) 	((ref) = (val1)-(val2))
#define value_sub_int(ref,val1,val2) 	((ref) = (val1)-(val2))
#define value_decrement(ref,val) 	((ref) = (val)-1)
#define value_division(ref,val1,val2)   ((ref) = cln::truncate1(val1,val2))
#define value_modulus(ref,val1,val2)    ((ref) = cln::truncate2(val1,val2).remainder)
#define value_pdivision(ref,val1,val2)  ((ref) = cln::floor1(val1,val2))
#define value_pmodulus(ref,val1,val2)   ((ref) = cln::floor2(val1,val2).remainder)
#define value_oppose(ref,val)    	((ref) = -(val))
#define value_absolute(ref,val)		((ref) = cln::abs(val))
#define value_minimum(ref,val1,val2)	((ref) = cln::min((val1),(val2)))
#define value_maximum(ref,val1,val2)	((ref) = cln::max((val1),(val2)))
#define value_orto(ref,val1,val2)	((ref) = (val1)|(val2))
#define value_andto(ref,val1,val2)	((ref) = (val1)&(val2))

/* Conditional operations on 'Value' */

#define value_pos_p(val)         ((val) >  0)
#define value_neg_p(val)         ((val) <  0)
#define value_posz_p(val)        ((val) >= 0)
#define value_negz_p(val)        ((val) <= 0)
#define value_zero_p(val)        ((val) == 0)
#define value_notzero_p(val)     ((val) != 0)
#define value_one_p(val)         ((val) == 1)
#define value_notone_p(val)      ((val) != 1)
#define value_mone_p(val)        ((val) == -1)
#define value_notmone_p(val)     ((val) != -1)
#define value_cmp_si(val, n)     (cln::compare(val,n))

#elif defined(GNUMP)

/* Basic macros */

#define value_init(val)        (mpz_init((val)))
#define value_assign(v1,v2)    (mpz_set((v1),(v2)))
#define value_set_si(val,i)    (mpz_set_si((val),(i)))    
#define value_set_double(val,d)(mpz_set_d((val),(d)))
#define value_clear(val)       (mpz_clear((val)))
#define value_read(val,str)    (mpz_set_str((val),(str),10))
#define value_print(Dst,fmt,val)  {char *str; \
				void (*gmp_free) (void *, size_t); \
				str = mpz_get_str(0,10,(val)); \
				fprintf((Dst),(fmt),str); \
				mp_get_memory_functions(NULL, NULL, &gmp_free); \
				(*gmp_free) (str, strlen(str)+1); \
                              }
#define value_swap(val1,val2) (mpz_swap(val1, val2))
                                             
/* Boolean operators on 'Value' */

#define value_eq(v1,v2) (mpz_cmp((v1),(v2)) == 0)
#define value_ne(v1,v2) (mpz_cmp((v1),(v2)) != 0)
#define value_gt(v1,v2) (mpz_cmp((v1),(v2))  > 0)
#define value_ge(v1,v2) (mpz_cmp((v1),(v2)) >= 0)
#define value_lt(v1,v2) (mpz_cmp((v1),(v2))  < 0)
#define value_le(v1,v2) (mpz_cmp((v1),(v2)) <= 0)

#define value_abs_eq(v1,v2) (mpz_cmpabs((v1),(v2)) == 0)
#define value_abs_ne(v1,v2) (mpz_cmpabs((v1),(v2)) != 0)
#define value_abs_gt(v1,v2) (mpz_cmpabs((v1),(v2))  > 0)
#define value_abs_ge(v1,v2) (mpz_cmpabs((v1),(v2)) >= 0)
#define value_abs_lt(v1,v2) (mpz_cmpabs((v1),(v2))  < 0)
#define value_abs_le(v1,v2) (mpz_cmpabs((v1),(v2)) <= 0)

/* Trian operators on 'Value' */

#define value_sign(val)      (mpz_sgn(val))
#define value_compare(v1,v2) (mpz_cmp((v1),(v2)))

/* Binary operations on 'Value' */

#define value_addto(ref,val1,val2)     (mpz_add((ref),(val1),(val2)))
#define value_add_int(ref,val,vint)     (mpz_add_ui((ref),(val),(long)(vint)))
#define value_addmul(ref, val1, val2)   (mpz_addmul((ref), (val1), (val2)))
#define value_increment(ref,val)       (mpz_add_ui((ref),(val),1))
#define value_multiply(ref,val1,val2)  (mpz_mul((ref),(val1),(val2)))
#define value_subtract(ref,val1,val2) (mpz_sub((ref),(val1),(val2)))
#define value_sub_int(ref,val,vint)     (mpz_sub_ui((ref),(val),(long)(vint)))
#define value_decrement(ref,val)       (mpz_sub_ui((ref),(val),1))
#define value_division(ref,val1,val2)  (mpz_tdiv_q((ref),(val1),(val2)))
#define value_modulus(ref,val1,val2)   (mpz_tdiv_r((ref),(val1),(val2)))
#define value_pdivision(ref,val1,val2) (mpz_fdiv_q((ref),(val1),(val2)))
#define value_pmodulus(ref,val1,val2)  (mpz_fdiv_r((ref),(val1),(val2)))
#define value_oppose(ref,val)          (mpz_neg((ref),(val)))
#define value_absolute(ref,val)        (mpz_abs((ref),(val)))
#define value_minimum(ref,val1,val2)   (value_le((val1),(val2)) ?  \
                                        mpz_set((ref),(val1)) :    \
                                        mpz_set((ref),(val2)))  
#define value_maximum(ref,val1,val2)   (value_ge((val1),(val2)) ?  \
                                        mpz_set((ref),(val1)) :    \
                                        mpz_set((ref),(val2)))  
#define value_orto(ref,val1,val2)      (mpz_ior((ref),(val1),(val2)))
#define value_andto(ref,val1,val2)     (mpz_and((ref),(val1),(val2)))

/* Conditional operations on 'Value' */

#define value_pos_p(val)         (mpz_sgn(val) >  0)
#define value_neg_p(val)         (mpz_sgn(val) <  0)
#define value_posz_p(val)        (mpz_sgn(val) >= 0)
#define value_negz_p(val)        (mpz_sgn(val) <= 0)
#define value_zero_p(val)        (mpz_sgn(val) == 0)
#define value_notzero_p(val)     (mpz_sgn(val) != 0)
#define value_one_p(val)         (mpz_cmp_si(val,1) == 0)
#define value_notone_p(val)      (mpz_cmp_si(val,1) != 0)
#define value_mone_p(val)        (mpz_cmp_si(val,-1) ==0)
#define value_notmone_p(val)     (mpz_cmp_si(val,-1) !=0)
#define value_cmp_si(val, n)     (mpz_cmp_si(val,n))

/* ************************************************************************* */

#else /* 'Value' set to longlong|long|float|char *|int */                                     	
/* Basic Macros */    				    

#define value_init(val)            ((val) = 0)
#define value_assign(v1,v2)        ((v1)  = (v2))
#define value_set_si(val,i)        ((val) = (Value)(i))   
#define value_set_double(val,d)    ((val) = (Value)(d)) 
#define value_clear(val)           ((val) = 0)
#define value_read(val,str)        (sscanf((str),VALUE_FMT,&(val)))
#define value_print(Dst,fmt,val)   (fprintf((Dst),(fmt),(val)))
#define value_swap(v1,v2)          {Value tmp; tmp = v2; \
                                    v2 = v1; v1 = tmp;   \
                                   }
/* Cast to 'Value' */

#define int_to_value(i) ((Value)(i))
#define long_to_value(l) ((Value)(l))
#define float_to_value(f) ((Value)(f))
#define double_to_value(d) ((Value)(d))
   
/* Boolean operators on 'Value' */

#define value_eq(v1,v2) ((v1)==(v2))
#define value_ne(v1,v2) ((v1)!=(v2))
#define value_gt(v1,v2) ((v1)>(v2))
#define value_ge(v1,v2) ((v1)>=(v2))
#define value_lt(v1,v2) ((v1)<(v2))
#define value_le(v1,v2) ((v1)<=(v2))

#define value_abs_eq(v1,v2) (value_abs(v1)==value_abs(v2))
#define value_abs_ne(v1,v2) (value_abs(v1)!=value_abs(v2))
#define value_abs_gt(v1,v2) (value_abs(v1)>value_abs(v2))
#define value_abs_ge(v1,v2) (value_abs(v1)>=value_abs(v2))
#define value_abs_lt(v1,v2) (value_abs(v1)<value_abs(v2))
#define value_abs_le(v1,v2) (value_abs(v1)<=value_abs(v2))

/* Trian operators on 'Value' */

#define value_sign(v) (value_eq(v,VALUE_ZERO)?0:value_lt(v,VALUE_ZERO)?-1:1)
#define value_compare(v1,v2) (value_eq(v1,v2)?0:value_lt(v1,v2)?-1:1)

/* Binary operators on 'Value' */

#define value_plus(v1,v2)  		((v1)+(v2))
#define value_div(v1,v2)   		((v1)/(v2))
#define value_mod(v1,v2)   		((v1)%(v2))
#define value_direct_multiply(v1,v2)	((v1)*(v2)) /* direct! */
#define value_minus(v1,v2) 		((v1)-(v2))
#define value_pdiv(v1,v2)  		(DIVIDE((v1),(v2)))
#define value_pmod(v1,v2)  		(MODULO((v1),(v2)))
#define value_min(v1,v2)   		(value_le((v1),(v2))? (v1): (v2))
#define value_max(v1,v2)   		(value_ge((v1),(v2))? (v1): (v2))
#define value_or(v1,v2)  		((v1)|(v2))
#define value_and(v1,v2)  		((v1)&(v2))
#define value_lshift(v1,v2)     	((v1)<<(v2))
#define value_rshift(v1,v2)  	        ((v1)>>(v2))
				  
/* Binary operations on 'Value' */ 

#define value_addto(ref,val1,val2) 	((ref) = (val1)+(val2))
#define value_add_int(ref,val,vint)     ((ref) = (val)+(Value)(vint))
#define value_addmul(ref, val1, val2)   ((ref) += (val1)*(val2))
#define value_increment(ref,val) 	((ref) = (val)+VALUE_ONE)
#define value_direct_product(ref,val1,val2) ((ref) = (val1)*(val2)) /* direct! */
#define value_multiply(ref,val1,val2)	((ref) = value_mult((val1),(val2)))
#define value_subtract(ref,val1,val2) 	((ref) = (val1)-(val2))
#define value_sub_int(ref,val,vint)     ((ref) = (val)-(Value)(vint))
#define value_decrement(ref,val) 	((ref) = (val)-VALUE_ONE)
#define value_division(ref,val1,val2) 	((ref) = (val1)/(val2))
#define value_modulus(ref,val1,val2) 	((ref) = (val1)%(val2))
#define value_pdivision(ref,val1,val2)	((ref) = value_pdiv((val1),(val2)))
#define value_pmodulus(ref,val1,val2)	((ref) = value_pmod((val1),(val2)))
#define value_oppose(ref,val)    	((ref) = value_uminus((val)))
#define value_absolute(ref,val)		((ref) = value_abs((val)))
#define value_minimum(ref,val1,val2)	((ref) = value_min((val1),(val2)))
#define value_maximum(ref,val1,val2)	((ref) = value_max((val1),(val2)))
#define value_orto(ref,val1,val2)	((ref) = (val1)|(val2))
#define value_andto(ref,val1,val2)	((ref) = (val1)&(val2))

/* Unary operators on 'Value' */

#define value_uminus(val)  (-(val))
#define value_not(val)	(~(val))
#define value_abs(val) (value_posz_p(val)? \
    (val) :                                \
    (value_ne((val), VALUE_NAN) ?          \
     value_uminus(val) :                   \
    (THROW (overflow_error), VALUE_NAN )))

/* Conditional operations on 'Value' */

#define value_pos_p(val)      value_gt(val,VALUE_ZERO)
#define value_neg_p(val)      value_lt(val,VALUE_ZERO)
#define value_posz_p(val)     value_ge(val,VALUE_ZERO)
#define value_negz_p(val)     value_le(val,VALUE_ZERO)
#define value_zero_p(val)     value_eq(val,VALUE_ZERO)
#define value_notzero_p(val)  value_ne(val,VALUE_ZERO)
#define value_one_p(val)      value_eq(val,VALUE_ONE)
#define value_notone_p(val)   value_ne(val,VALUE_ONE)
#define value_mone_p(val)     value_eq(val,VALUE_MONE)
#define value_notmone_p(val)  value_ne(val,VALUE_MONE)
#define value_cmp_si(val, n)  (val - VALUE_CONST(n))
#define value_min_p(val)      value_eq(val,VALUE_MIN)
#define value_max_p(val)      value_eq(val,VALUE_MAX)
#define value_notmin_p(val)   value_ne(val,VALUE_MIN)
#define value_notmax_p(val)   value_ne(val,VALUE_MAX)

#endif /* 'Value' set to |longlong|long|float|char *|int */


/* *********************** PROTECTED MULTIPLICATION ********************** */

#include "arithmetic_errors.h"

/* (|v| < MAX / |w|) => v*w is okay
 * I could check ((v*w)/w)==v but a tmp would be useful
 */
#define value_protected_hard_idiv_multiply(v,w,throw)		\
  ((value_zero_p(w) || value_zero_p(v))? VALUE_ZERO:		\
   value_lt(value_abs(v),value_div(VALUE_MAX,value_abs(w)))?	\
   value_direct_multiply(v,w): (throw, VALUE_NAN))

/* is a software idiv is assumed, quick check performed first
 */
#if defined(LINEAR_VALUE_ASSUME_SOFTWARE_IDIV)
#define value_protected_multiply(v,w,throw)				      \
  ((value_le(v,VALUE_SQRT_MAX) && value_le(w,VALUE_SQRT_MAX) &&		      \
   value_ge(v,VALUE_SQRT_MIN) && value_ge(w,VALUE_SQRT_MIN))?		      \
   value_direct_multiply(v,w): value_protected_hard_idiv_multiply(v,w,throw))
#else
#define value_protected_multiply(v,w,throw)		\
   value_protected_hard_idiv_multiply(v,w,throw)
#endif

/* protected versions
 */
#define value_protected_mult(v,w) 				\
    value_protected_multiply(v,w,THROW(overflow_error))
#define value_protected_product(v,w)		\
    v=value_protected_mult(v,w)

/* whether the default is protected or not 
 * this define makes no sense any more... well, doesn't matter. FC.
 */
#if defined(LINEAR_VALUE_PROTECT_MULTIPLY)
#define value_mult(v,w) value_protected_mult(v,w)
#define value_product(v,w) value_protected_product(v,w)
#else

/* I do enforce the protection whatever requested:-)
 * prints out a message and throws the exception, hoping
 * that some valid CATCH waits for it upwards. 
 */
#define value_mult(v,w)							      \
  value_protected_multiply(v,w,						      \
    (fprintf(stderr,"[value_mult] value overflow!\n"),THROW(overflow_error)))
#define value_product(v,w) v=value_mult(v,w)

/* was:
 * #define value_mult(v,w) value_direct_multiply(v,w)
 * #define value_product(v,w) value_direct_product(v,w)
 * could be: protected versions...
 */
#endif

/******************************************************* STATIC VALUE DEBUG */

/* LINEAR_VALUE_IS_CHARS is used for type checking.
 * some operations are not allowed on (char*), thus
 * they are switched to some other operation here...
 */
#if defined(LINEAR_VALUE_IS_CHARS)
#define value_fake_binary(v1,v2) ((Value)((v1).i+(v2).i))
#define value_bool_binary(v1,v2) ((int)((v1).i+(v2).i))
#undef float_to_value
#define float_to_value(f) ((Value)f)
#undef double_to_value
#define double_to_value(f) ((Value)f)
#undef value_uminus
#define value_uminus(v) (v)
#undef value_mult
#define value_mult(v1,v2) value_fake_binary(v1,v2)
#undef value_mod
#define value_mod(v1,v2) value_fake_binary(v1,v2)
#undef value_ge
#define value_ge(v1,v2) value_bool_binary(v1,v2)
#undef value_gt
#define value_gt(v1,v2) value_bool_binary(v1,v2)
#undef value_le
#define value_le(v1,v2) value_bool_binary(v1,v2)
#undef value_lt
#define value_lt(v1,v2) value_bool_binary(v1,v2)
#undef value_ne
#define value_ne(v1,v2) value_bool_binary(v1,v2)
#undef value_eq
#define value_eq(v1,v2) value_bool_binary(v1,v2)
#undef value_plus
#define value_plus(v1,v2) value_fake_binary(v1,v2)
#undef value_minus
#define value_minus(v1,v2) value_fake_binary(v1,v2)
#undef value_pdiv
#define value_pdiv(v1,v2) value_fake_binary(v1,v2)
#undef value_div
#define value_div(v1,v2) value_fake_binary(v1,v2)
#undef value_mod
#define value_mod(v1,v2) value_fake_binary(v1,v2)
#undef value_addto
#define value_addto(v1,v2) value_assign(v1,value_plus(v1,v2))
#undef value_subtract
#define value_subtract(v1,v2) value_addto(v1,v2)
#undef value_product
#define value_product(v1,v2) value_addto(v1,v2)
#undef value_modulus
#define value_modulus(v1,v2) value_addto(v1,v2)
#undef value_division
#define value_division(v1,v2) value_addto(v1,v2)
#undef value_increment
#define value_increment(v) value_addto(v,VALUE_ONE)
#undef value_decrement
#define value_decrement(v) value_addto(v,VALUE_MONE)
#undef value_orto
#define value_orto(ref,val) value_addto(v1,v2)
#undef value_andto
#define value_andto(ref,val) value_addto(v1,v2)	
#undef value_or
#define value_or(v1,v2) value_fake_binary(v1,v2)
#undef value_and
#define value_and(v1,v2) value_fake_binary(v1,v2)
#undef value_lshift
#define value_lshift(v1,v2) value_fake_binary(v1,v2)
#undef value_rshift
#define value_rshift(v1,v2) value_fake_binary(v1,v2)
#endif 

/* for backward compatibility */
#define value_substract(ref,val1,val2) (value_subtract((ref),(val1),(val2)))

/* valeur absolue
 */
#ifndef ABS
#define ABS(x) (((x)>=0) ? (x) : -(x))
#endif

/* minimum et maximum 
 * if they are defined somewhere else, they are very likely 
 * to be defined the same way. Thus the previous def is not overwritten.
 */
#ifndef MIN
#define MIN(x,y) (((x)>=(y))?(y):(x))
#endif
#ifndef MAX
#define MAX(x,y) (((x)>=(y))?(x):(y))
#endif

/* signe d'un entier: -1, 0 ou 1 */
#define SIGN(x) (((x)>0)? 1 : ((x)==0? 0 : -1))

/* division avec reste toujours positif
 * basee sur les equations:
 * a/(-b) = - (a/b)
 * (-a)/b = - ((a+b-1)/b)
 * ou a et b sont des entiers positifs
 */
#define DIVIDE(x,y) ((y)>0? POSITIVE_DIVIDE(x,y) : \
		     -POSITIVE_DIVIDE((x),(-(y))))

/* division avec reste toujours positif quand y est positif: assert(y>=0) */
#define POSITIVE_DIVIDE(x,y) ((x)>0 ? (x)/(y) : - (-(x)+(y)-1)/(y))

/* modulo a resultat toujours positif */
#define MODULO(x,y) ((y)>0 ? POSITIVE_MODULO(x,y) : POSITIVE_MODULO(-x,-y))

/* modulo par rapport a un nombre positif: assert(y>=0)
 *
 * Ce n'est pas la macro la plus efficace que j'aie jamais ecrite: il faut
 * faire, dans le pire des cas, deux appels a la routine .rem, qui n'est
 * surement pas plus cablee que la division ou la multiplication
 */
#define POSITIVE_MODULO(x,y) ((x) > 0 ? (x)%(y) : \
			      ((x)%(y) == 0 ? 0 : ((y)-(-(x))%(y))))
			      
/* errors.c */ 
extern unsigned int overflow_error;
extern unsigned int simplex_arithmetic_error;
extern unsigned int user_exception_error;
extern unsigned int parser_exception_error;
extern unsigned int any_exception_error; 
extern unsigned int the_last_just_thrown_exception;
extern void dump_exception_stack_to_file(FILE * /*f*/);
extern void dump_exception_stack(void);
extern jmp_buf *push_exception_on_stack(int /*what*/, const char * /*function*/, const char * /*file*/, int /*line*/);
extern void pop_exception_from_stack(int /*what*/, const char * /*function*/, const char * /*file*/, int /*line*/);
extern void throw_exception(int /*what*/, const char * /*function*/, const char * /*file*/, int /*line*/);

#endif /* arithmetique_header_included */



