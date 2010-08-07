
   /*+------- <| --------------------------------------------------------**
    **         A                  Clan/Scop                              **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 macros.h                              **
    **-   A   | #   -----------------------------------------------------**
    **   /.\ [""M#         First version: 30/04/2008                     **
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
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 2.1 of the     *
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

#ifndef SCOPLIB_MACROS_H
# define SCOPLIB_MACROS_H


# if defined(SCOPLIB_INT_T_IS_LONGLONG)
#  define SCOPLIB_FMT     "%4lld"
#  define SCOPLIB_FMT_TXT "%lld"
#  define scoplib_int_t long long

# elif defined(SCOPLIB_INT_T_IS_LONG)
#  define SCOPLIB_FMT     "%4ld"
#  define SCOPLIB_FMT_TXT "%ld"
#  define scoplib_int_t long int

# elif defined(SCOPLIB_INT_T_IS_MP)  /* GNUMP */
#include <gmp.h>
#  define SCOPLIB_FMT     "%4s"
#  define SCOPLIB_FMT_TXT "%s"
#  define scoplib_int_t mpz_t

# else
#  error Define SCOPLIB_INT_T_IS_xxx to use this file.

# endif

# define SCOPLIB_DEBUG			0 /* Set to 1 for debug mode,
					     0 otherwise */
# define SCOPLIB_MAX_STRING		2048
# define SCOPLIB_TYPE_ITERATOR		1
# define SCOPLIB_TYPE_PARAMETER		2
# define SCOPLIB_TYPE_ARRAY		3
# define SCOPLIB_TYPE_FUNCTION		4
# define SCOPLIB_TYPE_DOMAIN		6
# define SCOPLIB_TYPE_SCATTERING	7
# define SCOPLIB_TYPE_ACCESS		8
# define SCOPLIB_TYPE_UNKNOWN		9
# define SCOPLIB_FAKE_ARRAY		"fakearray"

# define SCOPLIB_SCOP_PRINT_CASTLE	1
# define SCOPLIB_SCOP_PRINT_ARRAYSTAG	2


/*+****************************************************************************
 *                              SCOP GMP MACROS                               *
 ******************************************************************************/
# ifdef SCOPLIB_INT_T_IS_MP
/* Basic Macros */
#  define SCOPVAL_init(val)                (mpz_init((val)))
#  define SCOPVAL_assign(v1,v2)            (mpz_set((v1),(v2)))
#  define SCOPVAL_set_si(val,i)            (mpz_set_si((val),(i)))
#  define SCOPVAL_get_si(val)              (mpz_get_si((val)))
#  define SCOPVAL_init_set_si(val,i)       (mpz_init_set_si((val),(i)))
#  define SCOPVAL_clear(val)               (mpz_clear((val)))
#  define SCOPVAL_print(Dst,fmt,val)       { char *str; \
                                        str = mpz_get_str(0,10,(val)); \
                                        fprintf((Dst),(fmt),str); free(str); \
                                        }
#  define SCOPVAL_sprint(Dst,fmt,val)      { char * str; \
                                        str = mpz_get_str(0,10,(val)); \
                                        sprintf((Dst),(fmt),str); free(str); \
                                        }

/* Boolean operators on 'scoplib_int_t' */
#  define SCOPVAL_eq(v1,v2)                (mpz_cmp((v1),(v2)) == 0)
#  define SCOPVAL_ne(v1,v2)                (mpz_cmp((v1),(v2)) != 0)

/* Binary operators on 'scoplib_int_t' */
#  define SCOPVAL_increment(ref,val)       (mpz_add_ui((ref),(val),1))
#  define SCOPVAL_addto(ref,val1,val2)     (mpz_add((ref),(val1),(val2)))
#  define SCOPVAL_multo(ref,val1,val2)     (mpz_mul((ref),(val1),(val2)))
#  define SCOPVAL_add_int(ref,val,vint)    (mpz_add_ui((ref),(val),(long)(vint)))
#  define SCOPVAL_subtract(ref,val1,val2)  (mpz_sub((ref),(val1),(val2)))
#  define SCOPVAL_oppose(ref,val)          (mpz_neg((ref),(val)))

/* Conditional operations on 'scoplib_int_t' */
#  define SCOPVAL_pos_p(val)               (mpz_sgn(val) >  0)
#  define SCOPVAL_neg_p(val)               (mpz_sgn(val) <  0)
#  define SCOPVAL_zero_p(val)              (mpz_sgn(val) == 0)
#  define SCOPVAL_notzero_p(val)           (mpz_sgn(val) != 0)
#  define SCOPVAL_one_p(val)               (mpz_cmp_si(val,1)  == 0)
#  define SCOPVAL_mone_p(val)              (mpz_cmp_si(val,-1) == 0)

/*+****************************************************************************
 *                           SCOPVAL BASIC TYPES MACROS                          *
 ******************************************************************************/
# else
/* Basic Macros */
#  define SCOPVAL_init(val)                ((val) = 0)
#  define SCOPVAL_assign(v1,v2)            ((v1)  = (v2))
#  define SCOPVAL_set_si(val,i)            ((val) = (scoplib_int_t)(i))
#  define SCOPVAL_get_si(val)              ((val))
#  define SCOPVAL_init_set_si(val,i)       ((val) = (scoplib_int_t)(i))
#  define SCOPVAL_clear(val)               ((val) = 0)
#  define SCOPVAL_print(Dst,fmt,val)       (fprintf((Dst),(fmt),(val)))
#  define SCOPVAL_sprint(Dst,fmt,val)      (sprintf((Dst),(fmt),(val)))

/* Boolean operators on 'scoplib_int_t' */
#  define SCOPVAL_eq(v1,v2)                ((v1)==(v2))
#  define SCOPVAL_ne(v1,v2)                ((v1)!=(v2))
#  define SCOPVAL_lt(v1,v2)                ((v1)<(v2))
#  define SCOPVAL_gt(v1,v2)                ((v1)>(v2))

/* Binary operators on 'scoplib_int_t' */
#  define SCOPVAL_increment(ref,val)       ((ref) = (val)+(scoplib_int_t)(1))
#  define SCOPVAL_addto(ref,val1,val2)     ((ref) = (val1)+(val2))
#  define SCOPVAL_multo(ref,val1,val2)     ((ref) = (val1)*(val2))
#  define SCOPVAL_add_int(ref,val,vint)    ((ref) = (val)+(scoplib_int_t)(vint))
#  define SCOPVAL_subtract(ref,val1,val2)  ((ref) = (val1)-(val2))
#  define SCOPVAL_oppose(ref,val)          ((ref) = (-(val)))

/* Conditional operations on 'scoplib_int_t' */
#  define SCOPVAL_pos_p(val)               SCOPVAL_gt(val,0)
#  define SCOPVAL_neg_p(val)               SCOPVAL_lt(val,0)
#  define SCOPVAL_zero_p(val)              SCOPVAL_eq(val,0)
#  define SCOPVAL_notzero_p(val)           SCOPVAL_ne(val,0)
#  define SCOPVAL_one_p(val)               SCOPVAL_eq(val,1)
#  define SCOPVAL_mone_p(val)              SCOPVAL_eq(val,-1)

# endif

#endif /* define SCOPLIB_MACROS_H */
