
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.4.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Copy the first part of user declarations.  */

/* Line 189 of yacc.c  */
#line 39 "parser.y"

   #include <stdio.h>
   #include <stdlib.h>
   #include <string.h>
   #include <assert.h>
   #include <clan/macros.h>
   #include <clan/vector.h>
   #include <clan/matrix.h>
   #include <clan/scop.h>
   #include <clan/symbol.h>
   #include <clan/statement.h>
   #include <clan/options.h>

   int yylex(void);
   void yyerror(char *);
   int clan_parse_error = 0; /**< Set to 1 during parsing if
				encountered an error */
   void clan_parser_log(char *);
   scoplib_scop_p clan_parse(FILE *, clan_options_p);

   extern FILE * yyin;                  /**< File to be read by Lex */
   extern char scanner_latest_text[];   /**< Latest text read by Lex */

   /* This is the "parser state", a collection of variables that vary
    * during the parsing and thanks to we can extract all SCoP informations.
    */
   scoplib_scop_p      parser_scop;        /**< SCoP in construction */
   scoplib_statement_p parser_statement;   /**< Statement in construction */
   clan_symbol_p       parser_symbol;      /**< Top of the symbol table */
   int                 parser_recording;   /**< Boolean: do we record or not? */
   char *              parser_record;      /**< What we record
					      (statement body) */
   int                 parser_depth = 0;   /**< Current loop depth */
   int *               parser_scheduling;  /**< Current statement scheduling */
   clan_symbol_p *     parser_iterators;   /**< Current iterator list */
   scoplib_matrix_p    parser_domain;      /**< Current iteration domain */
   int                 parser_nb_cons = 0; /**< Current number of constraints */
   int *               parser_consperdim;  /**< Constraint nb for each
					      dimension */
   int*		       parser_variables_localvars;/**< List of variables
						     in #pragma
						     local-vars */
   int*		       parser_variables_liveout;/**< List of variables
						     in #pragma
						     live-out */
   /* Ugly global variable to keep/read Clan options during parsing. */
   clan_options_p	parser_options = NULL;




/* Line 189 of yacc.c  */
#line 125 "parser.c"

/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     IGNORE = 258,
     IF = 259,
     ELSE = 260,
     FOR = 261,
     PRAGMALOCALVARS = 262,
     PRAGMALIVEOUT = 263,
     MIN = 264,
     MAX = 265,
     CEILD = 266,
     FLOORD = 267,
     REAL = 268,
     ID = 269,
     INTEGER = 270,
     syRPARENTHESIS = 271,
     syLPARENTHESIS = 272,
     syRBRACKET = 273,
     syLBRACKET = 274,
     syRBRACE = 275,
     syLBRACE = 276,
     sySEMICOLON = 277,
     syCOMMA = 278,
     syPOINT = 279,
     syARROW = 280,
     opEQUAL = 281,
     opLEQ = 282,
     opGEQ = 283,
     opLOWER = 284,
     opGREATER = 285,
     opPLUS = 286,
     opMINUS = 287,
     opINCREMENTATION = 288,
     opDECREMENTATION = 289,
     opNOT = 290,
     opMULTIPLY = 291,
     opDIVIDE = 292,
     opMOD = 293,
     opAND = 294,
     opOR = 295,
     opCOMP = 296,
     opASSIGNMENT = 297,
     opPLUSEQUAL = 298,
     opMINUSEQUAL = 299,
     opMULTIPLYEQUAL = 300,
     opDIVIDEEQUAL = 301,
     opMODEQUAL = 302,
     opANDEQUAL = 303,
     opOREQUAL = 304,
     opCOMPEQUAL = 305,
     opLAND = 306,
     opLOR = 307,
     opQMARK = 308,
     opCOLON = 309,
     MAXPRIORITY = 310
   };
#endif
/* Tokens.  */
#define IGNORE 258
#define IF 259
#define ELSE 260
#define FOR 261
#define PRAGMALOCALVARS 262
#define PRAGMALIVEOUT 263
#define MIN 264
#define MAX 265
#define CEILD 266
#define FLOORD 267
#define REAL 268
#define ID 269
#define INTEGER 270
#define syRPARENTHESIS 271
#define syLPARENTHESIS 272
#define syRBRACKET 273
#define syLBRACKET 274
#define syRBRACE 275
#define syLBRACE 276
#define sySEMICOLON 277
#define syCOMMA 278
#define syPOINT 279
#define syARROW 280
#define opEQUAL 281
#define opLEQ 282
#define opGEQ 283
#define opLOWER 284
#define opGREATER 285
#define opPLUS 286
#define opMINUS 287
#define opINCREMENTATION 288
#define opDECREMENTATION 289
#define opNOT 290
#define opMULTIPLY 291
#define opDIVIDE 292
#define opMOD 293
#define opAND 294
#define opOR 295
#define opCOMP 296
#define opASSIGNMENT 297
#define opPLUSEQUAL 298
#define opMINUSEQUAL 299
#define opMULTIPLYEQUAL 300
#define opDIVIDEEQUAL 301
#define opMODEQUAL 302
#define opANDEQUAL 303
#define opOREQUAL 304
#define opCOMPEQUAL 305
#define opLAND 306
#define opLOR 307
#define opQMARK 308
#define opCOLON 309
#define MAXPRIORITY 310




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 90 "parser.y"
 int value;                     /**< An integer value for integers */
         char * symbol;                 /**< A string for identifiers */
         scoplib_vector_p affex;        /**< An affine expression */
         scoplib_matrix_p setex;        /**< A set of affine expressions */
         scoplib_matrix_p rw[2];        /**< Read and write array accesses */
       


/* Line 214 of yacc.c  */
#line 280 "parser.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 292 "parser.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  32
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   569

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  56
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  29
/* YYNRULES -- Number of rules.  */
#define YYNRULES  114
/* YYNRULES -- Number of states.  */
#define YYNSTATES  216

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   310

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     5,     6,     8,    11,    13,    16,    20,
      22,    26,    27,    28,    29,    30,    46,    47,    54,    55,
      58,    61,    64,    67,    70,    76,    80,    82,    89,    91,
      98,   100,   104,   108,   112,   119,   126,   128,   130,   133,
     137,   141,   145,   150,   153,   157,   161,   165,   169,   173,
     177,   181,   183,   185,   187,   189,   191,   193,   195,   197,
     199,   201,   206,   211,   215,   219,   222,   226,   228,   230,
     232,   234,   236,   238,   240,   242,   244,   246,   248,   250,
     252,   254,   257,   259,   263,   266,   270,   276,   278,   281,
     286,   289,   292,   294,   298,   302,   306,   310,   314,   316,
     320,   324,   326,   327,   331,   336,   338,   342,   345,   347,
     349,   351,   353,   355,   357
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      57,     0,    -1,    58,    -1,    -1,    60,    -1,    58,    60,
      -1,     3,    -1,    58,     3,    -1,    20,    58,    21,    -1,
      60,    -1,    20,    58,    21,    -1,    -1,    -1,    -1,    -1,
       6,    16,    82,    61,    42,    69,    62,    22,    72,    63,
      22,    67,    17,    64,    59,    -1,    -1,     4,    16,    72,
      17,    65,    59,    -1,    -1,    66,    75,    -1,     7,    80,
      -1,     8,    80,    -1,    82,    33,    -1,    33,    82,    -1,
      82,    42,    82,    31,    15,    -1,    82,    43,    15,    -1,
      70,    -1,     9,    16,    68,    23,    68,    17,    -1,    70,
      -1,    10,    16,    69,    23,    69,    17,    -1,    71,    -1,
      70,    31,    70,    -1,    70,    32,    70,    -1,    16,    70,
      17,    -1,    11,    16,    70,    23,    71,    17,    -1,    12,
      16,    70,    23,    71,    17,    -1,    15,    -1,    82,    -1,
      32,    15,    -1,    15,    36,    82,    -1,    15,    36,    15,
      -1,    15,    37,    15,    -1,    32,    15,    36,    82,    -1,
      32,    82,    -1,    70,    29,    68,    -1,    70,    30,    69,
      -1,    70,    27,    68,    -1,    70,    28,    69,    -1,    70,
      26,    70,    -1,    16,    72,    17,    -1,    72,    51,    72,
      -1,    43,    -1,    44,    -1,    45,    -1,    46,    -1,    47,
      -1,    48,    -1,    49,    -1,    50,    -1,    33,    -1,    34,
      -1,    78,    42,    77,    22,    -1,    78,    73,    77,    22,
      -1,    78,    74,    22,    -1,    74,    78,    22,    -1,    78,
      22,    -1,    20,    75,    21,    -1,    31,    -1,    32,    -1,
      36,    -1,    37,    -1,    38,    -1,    28,    -1,    30,    -1,
      27,    -1,    29,    -1,    26,    -1,    39,    -1,    40,    -1,
      41,    -1,    84,    -1,    32,    84,    -1,    78,    -1,    77,
      76,    77,    -1,    35,    77,    -1,    16,    77,    17,    -1,
      77,    53,    77,    54,    77,    -1,    82,    -1,    82,    81,
      -1,    82,    16,    80,    17,    -1,    32,    78,    -1,    31,
      78,    -1,    84,    -1,    79,    32,    79,    -1,    79,    31,
      79,    -1,    79,    36,    79,    -1,    79,    37,    79,    -1,
      16,    79,    17,    -1,    78,    -1,    80,    23,    78,    -1,
      80,    23,    79,    -1,    79,    -1,    -1,    18,    70,    19,
      -1,    81,    18,    70,    19,    -1,    14,    -1,    16,    14,
      17,    -1,    39,    14,    -1,    83,    -1,     9,    -1,    10,
      -1,    11,    -1,    12,    -1,    15,    -1,    13,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   140,   140,   169,   177,   178,   179,   180,   181,   192,
     196,   213,   235,   255,   263,   210,   280,   279,   298,   298,
     366,   390,   419,   423,   427,   437,   455,   460,   474,   479,
     493,   497,   503,   509,   513,   518,   535,   542,   551,   558,
     567,   574,   581,   590,   610,   645,   682,   715,   750,   767,
     774,   788,   789,   790,   791,   792,   793,   794,   795,   804,
     805,   818,   831,   845,   858,   871,   879,   892,   893,   894,
     895,   896,   897,   898,   899,   900,   901,   902,   903,   904,
     915,   922,   929,   937,   946,   953,   960,   982,  1008,  1021,
    1029,  1036,  1047,  1048,  1049,  1050,  1051,  1052,  1064,  1071,
    1078,  1085,  1093,  1108,  1116,  1137,  1144,  1151,  1158,  1164,
    1164,  1164,  1164,  1167,  1168
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "IGNORE", "IF", "ELSE", "FOR",
  "PRAGMALOCALVARS", "PRAGMALIVEOUT", "MIN", "MAX", "CEILD", "FLOORD",
  "REAL", "ID", "INTEGER", "syRPARENTHESIS", "syLPARENTHESIS",
  "syRBRACKET", "syLBRACKET", "syRBRACE", "syLBRACE", "sySEMICOLON",
  "syCOMMA", "syPOINT", "syARROW", "opEQUAL", "opLEQ", "opGEQ", "opLOWER",
  "opGREATER", "opPLUS", "opMINUS", "opINCREMENTATION", "opDECREMENTATION",
  "opNOT", "opMULTIPLY", "opDIVIDE", "opMOD", "opAND", "opOR", "opCOMP",
  "opASSIGNMENT", "opPLUSEQUAL", "opMINUSEQUAL", "opMULTIPLYEQUAL",
  "opDIVIDEEQUAL", "opMODEQUAL", "opANDEQUAL", "opOREQUAL", "opCOMPEQUAL",
  "opLAND", "opLOR", "opQMARK", "opCOLON", "MAXPRIORITY", "$accept",
  "program", "instruction_list", "bloc", "instruction", "$@1", "$@2",
  "$@3", "$@4", "$@5", "$@6", "incrementation", "min_affine_expression",
  "max_affine_expression", "affine_expression", "term", "condition",
  "reduction_operator", "unary_operator", "assignment", "binary_operator",
  "expression", "variable", "arithmetic_expression", "variable_list",
  "array_index", "id", "math_func_list", "NUMBER", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    56,    57,    57,    58,    58,    58,    58,    58,    59,
      59,    61,    62,    63,    64,    60,    65,    60,    66,    60,
      60,    60,    67,    67,    67,    67,    68,    68,    69,    69,
      70,    70,    70,    70,    70,    70,    71,    71,    71,    71,
      71,    71,    71,    71,    72,    72,    72,    72,    72,    72,
      72,    73,    73,    73,    73,    73,    73,    73,    73,    74,
      74,    75,    75,    75,    75,    75,    75,    76,    76,    76,
      76,    76,    76,    76,    76,    76,    76,    76,    76,    76,
      77,    77,    77,    77,    77,    77,    77,    78,    78,    78,
      78,    78,    79,    79,    79,    79,    79,    79,    80,    80,
      80,    80,    80,    81,    81,    82,    82,    82,    82,    83,
      83,    83,    83,    84,    84
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     1,     0,     1,     2,     1,     2,     3,     1,
       3,     0,     0,     0,     0,    15,     0,     6,     0,     2,
       2,     2,     2,     2,     5,     3,     1,     6,     1,     6,
       1,     3,     3,     3,     6,     6,     1,     1,     2,     3,
       3,     3,     4,     2,     3,     3,     3,     3,     3,     3,
       3,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     4,     4,     3,     3,     2,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     2,     1,     3,     2,     3,     5,     1,     2,     4,
       2,     2,     1,     3,     3,     3,     3,     3,     1,     3,
       3,     1,     0,     3,     4,     1,     3,     2,     1,     1,
       1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
      18,     6,     0,     0,   102,   102,    18,     0,    18,     4,
       0,     0,     0,   109,   110,   111,   112,   114,   105,   113,
       0,     0,     0,     0,    98,   101,    20,    87,   108,    92,
      21,    18,     1,     7,     5,     0,     0,    59,    60,     0,
      19,     0,   111,   112,    36,     0,     0,     0,    30,     0,
      37,    11,     0,     0,     0,    91,    90,   107,     0,     0,
       0,     0,     0,   102,     0,    88,     8,     0,     0,    65,
       0,    51,    52,    53,    54,    55,    56,    57,    58,     0,
       0,     0,     0,     0,     0,   105,     0,     0,    38,    43,
       0,     0,     0,     0,     0,     0,     0,    16,     0,     0,
     106,    97,    94,    93,    95,    96,    99,   100,     0,     0,
       0,     0,    66,    64,     0,     0,     0,     0,    82,    80,
       0,    63,     0,     0,    40,    39,    41,    33,    49,     0,
      48,   109,    46,    26,   110,    47,    28,    44,    45,    31,
      32,    18,    50,     0,    89,     0,   103,     0,     0,    81,
      84,    61,    76,    74,    72,    75,    73,    67,    68,    69,
      70,    71,    77,    78,    79,     0,     0,    62,     0,     0,
      42,     0,     0,    18,    17,     9,    12,   104,    85,     0,
      83,     0,     0,     0,     0,    18,     0,     0,    34,    35,
       0,     0,    10,     0,    86,     0,     0,    13,    27,    29,
       0,     0,     0,     0,     0,    23,    14,    22,     0,     0,
      18,     0,    25,    15,     0,    24
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     7,     8,   174,     9,    99,   186,   200,   210,   141,
      10,   203,   132,   135,   136,    48,    49,    79,    39,    40,
     166,   117,   118,    25,    26,    65,    50,    28,    29
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -89
static const yytype_int16 yypact[] =
{
     480,   -89,   -14,    13,   286,   286,   181,    15,   142,   -89,
     234,   325,   104,   -89,   -89,   -89,   -89,   -89,   -89,   -89,
     114,   349,   349,    35,   -89,   120,   -10,   154,   -89,   -89,
     -10,   528,   -89,   -89,   -89,    39,   234,   -89,   -89,   349,
     -89,   477,    40,    51,   -17,   357,   436,   537,   -89,    -6,
     -89,   -89,    70,   178,    33,   -89,   -89,   -89,   178,   178,
     178,   178,   286,   286,   383,    86,   -89,    76,    85,   -89,
     246,   -89,   -89,   -89,   -89,   -89,   -89,   -89,   -89,   246,
     102,   383,   383,   454,   116,    70,   530,    -3,   111,   -89,
     383,   391,   402,   391,   402,   383,   383,   -89,   325,   117,
     -89,   -89,    26,    26,   -89,   -89,   -89,   120,    24,   417,
       6,   383,   -89,   -89,   277,   317,   246,   180,   -89,   -89,
     200,   -89,    20,   103,   -89,   -89,   -89,   -89,   -89,   104,
      46,   144,   -89,    46,   163,   -89,    46,   -89,   -89,   -89,
     -89,    78,   -89,   402,   -89,    -1,   -89,    57,   137,   -89,
     -89,   -89,   -89,   -89,   -89,   -89,   -89,   -89,   -89,   -89,
     -89,   -89,   -89,   -89,   -89,   246,   246,   -89,   428,   428,
     -89,   391,   402,   181,   -89,   -89,   -89,   -89,   -89,   476,
     -89,   136,   141,   159,   169,   534,   164,   246,   -89,   -89,
     391,   402,   -89,   325,   -89,   179,   182,   162,   -89,   -89,
     173,   462,   104,   197,    -9,   -89,   -89,   -89,   104,   208,
      78,   184,   -89,   -89,   209,   -89
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -89,   -89,    -2,    25,    -5,   -89,   -89,   -89,   -89,   -89,
     -89,   -89,   -88,   -82,    10,   -75,   -38,   -89,   193,   189,
     -89,   -70,    18,    79,     5,   -89,    -4,   -89,   -43
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -4
static const yytype_int16 yytable[] =
{
      27,    27,    11,    34,    31,   137,    27,    87,    51,   120,
      30,    97,   138,    62,   128,    32,   127,    27,    27,    83,
      84,    47,    24,    24,   207,   146,    34,   119,    41,    12,
      95,    96,    27,   208,   209,    27,   119,    95,    96,    55,
      56,   144,    89,   168,   148,    98,   150,    62,    98,    57,
     101,    95,    96,    52,    41,    86,    81,    68,    27,    27,
     142,   176,    60,    61,    58,    59,    27,    82,   108,    60,
      61,   119,   149,   119,   110,    27,   177,    95,    96,   125,
     106,    24,     2,   183,     3,     4,     5,   100,    95,    96,
     184,   122,   123,   181,   182,   179,   180,   112,   173,    54,
     130,   133,   195,   133,   111,   139,   140,   113,    47,   196,
      27,    27,    27,    13,    14,    15,    16,   194,    18,   145,
      35,   147,   119,   119,   121,   170,   169,    17,    52,    19,
      53,   126,    54,    56,    95,    96,   175,   102,   103,   104,
     105,   107,    -2,    23,   119,    33,     2,   129,     3,     4,
       5,    58,    59,   188,   178,   197,    60,    61,   189,   143,
     171,    27,    27,   152,   153,   154,   155,   156,   157,   158,
      63,   185,    64,   159,   160,   161,   162,   163,   164,   172,
      34,   133,   190,    27,     1,     2,   193,     3,     4,     5,
     165,    17,   191,    19,    53,   201,   198,   204,   205,   199,
     133,     6,   151,    47,   211,   175,   152,   153,   154,   155,
     156,   157,   158,    98,   206,   214,   159,   160,   161,   162,
     163,   164,   167,   212,   215,    67,   152,   153,   154,   155,
     156,   157,   158,   165,    80,   213,   159,   160,   161,   162,
     163,   164,     0,    13,    14,    15,    16,     0,    18,     0,
      35,     0,     0,   165,    36,    13,    14,    15,    16,    17,
      18,    19,   114,     0,     0,    21,    22,    37,    38,     0,
       0,     0,     0,    23,     0,     0,     0,    21,   115,     0,
       0,   116,     0,     0,     0,    23,    13,    14,    15,    16,
      17,    85,    19,   114,     0,    13,    14,    15,    16,    17,
      18,    19,    20,     0,     0,     0,     0,     0,    21,   115,
       0,     0,   116,     0,     0,     0,    23,    21,    22,     0,
       0,     0,     0,     0,     0,    23,    13,    14,    15,    16,
      17,    18,    19,    35,    13,    14,    42,    43,     0,    18,
      44,    45,     0,     0,     0,     0,     0,     0,    21,    22,
       0,     0,     0,     0,     0,     0,    23,    46,    13,    14,
      15,    16,     0,    18,    23,    35,    13,    14,    42,    43,
       0,    85,    44,    45,     0,     0,     0,     0,     0,     0,
      21,    22,     0,     0,     0,     0,     0,     0,    23,    46,
       0,     0,    13,    14,    42,    43,    23,    18,    44,   109,
     131,    14,    42,    43,     0,    18,    44,   109,     0,     0,
       0,    13,   134,    42,    43,    46,    18,    44,   109,     0,
       0,     0,    23,    46,     0,     0,    13,    14,    42,    43,
      23,    85,    44,   109,    46,     0,     0,    13,    14,    15,
      16,    23,    18,    44,    35,    13,    14,    15,    16,    46,
      18,    88,    35,     0,     0,     0,    23,     0,     0,     0,
      46,     0,     0,    13,    14,    15,    16,    23,    18,   124,
      35,    13,    14,    15,    16,    23,    18,     0,    35,     0,
      -3,     0,     0,     1,     2,     0,     3,     4,     5,     0,
       0,     0,     0,    23,     0,   202,     0,     0,     0,    69,
       6,    23,   152,   153,   154,   155,   156,   157,   158,     0,
      37,    38,   159,   160,   161,   162,   163,   164,     0,    70,
      71,    72,    73,    74,    75,    76,    77,    78,     0,   165,
     187,    33,     2,     0,     3,     4,     5,    33,     2,     0,
       3,     4,     5,     0,     0,     0,     0,   127,     0,    66,
       0,     0,     0,     0,     0,   192,    90,    91,    92,    93,
      94,    95,    96,    90,    91,    92,    93,    94,    95,    96
};

static const yytype_int16 yycheck[] =
{
       4,     5,    16,     8,     6,    93,    10,    45,    12,    79,
       5,    17,    94,    23,    17,     0,    17,    21,    22,    36,
      37,    11,     4,     5,    33,    19,    31,    70,    10,    16,
      31,    32,    36,    42,    43,    39,    79,    31,    32,    21,
      22,    17,    46,    23,   114,    51,   116,    23,    51,    14,
      17,    31,    32,    14,    36,    45,    16,    39,    62,    63,
      98,   143,    36,    37,    31,    32,    70,    16,    63,    36,
      37,   114,   115,   116,    64,    79,    19,    31,    32,    83,
      62,    63,     4,   171,     6,     7,     8,    17,    31,    32,
     172,    81,    82,   168,   169,   165,   166,    21,    20,    20,
      90,    91,   190,    93,    18,    95,    96,    22,    98,   191,
     114,   115,   116,     9,    10,    11,    12,   187,    14,   109,
      16,   111,   165,   166,    22,   129,    23,    13,    14,    15,
      16,    15,    53,   115,    31,    32,   141,    58,    59,    60,
      61,    62,     0,    39,   187,     3,     4,    36,     6,     7,
       8,    31,    32,    17,    17,   193,    36,    37,    17,    42,
      16,   165,   166,    26,    27,    28,    29,    30,    31,    32,
      16,   173,    18,    36,    37,    38,    39,    40,    41,    16,
     185,   171,    23,   187,     3,     4,    22,     6,     7,     8,
      53,    13,    23,    15,    16,    22,    17,   201,   202,    17,
     190,    20,    22,   193,   208,   210,    26,    27,    28,    29,
      30,    31,    32,    51,    17,    31,    36,    37,    38,    39,
      40,    41,    22,    15,    15,    36,    26,    27,    28,    29,
      30,    31,    32,    53,    41,   210,    36,    37,    38,    39,
      40,    41,    -1,     9,    10,    11,    12,    -1,    14,    -1,
      16,    -1,    -1,    53,    20,     9,    10,    11,    12,    13,
      14,    15,    16,    -1,    -1,    31,    32,    33,    34,    -1,
      -1,    -1,    -1,    39,    -1,    -1,    -1,    31,    32,    -1,
      -1,    35,    -1,    -1,    -1,    39,     9,    10,    11,    12,
      13,    14,    15,    16,    -1,     9,    10,    11,    12,    13,
      14,    15,    16,    -1,    -1,    -1,    -1,    -1,    31,    32,
      -1,    -1,    35,    -1,    -1,    -1,    39,    31,    32,    -1,
      -1,    -1,    -1,    -1,    -1,    39,     9,    10,    11,    12,
      13,    14,    15,    16,     9,    10,    11,    12,    -1,    14,
      15,    16,    -1,    -1,    -1,    -1,    -1,    -1,    31,    32,
      -1,    -1,    -1,    -1,    -1,    -1,    39,    32,     9,    10,
      11,    12,    -1,    14,    39,    16,     9,    10,    11,    12,
      -1,    14,    15,    16,    -1,    -1,    -1,    -1,    -1,    -1,
      31,    32,    -1,    -1,    -1,    -1,    -1,    -1,    39,    32,
      -1,    -1,     9,    10,    11,    12,    39,    14,    15,    16,
       9,    10,    11,    12,    -1,    14,    15,    16,    -1,    -1,
      -1,     9,    10,    11,    12,    32,    14,    15,    16,    -1,
      -1,    -1,    39,    32,    -1,    -1,     9,    10,    11,    12,
      39,    14,    15,    16,    32,    -1,    -1,     9,    10,    11,
      12,    39,    14,    15,    16,     9,    10,    11,    12,    32,
      14,    15,    16,    -1,    -1,    -1,    39,    -1,    -1,    -1,
      32,    -1,    -1,     9,    10,    11,    12,    39,    14,    15,
      16,     9,    10,    11,    12,    39,    14,    -1,    16,    -1,
       0,    -1,    -1,     3,     4,    -1,     6,     7,     8,    -1,
      -1,    -1,    -1,    39,    -1,    33,    -1,    -1,    -1,    22,
      20,    39,    26,    27,    28,    29,    30,    31,    32,    -1,
      33,    34,    36,    37,    38,    39,    40,    41,    -1,    42,
      43,    44,    45,    46,    47,    48,    49,    50,    -1,    53,
      54,     3,     4,    -1,     6,     7,     8,     3,     4,    -1,
       6,     7,     8,    -1,    -1,    -1,    -1,    17,    -1,    21,
      -1,    -1,    -1,    -1,    -1,    21,    26,    27,    28,    29,
      30,    31,    32,    26,    27,    28,    29,    30,    31,    32
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     4,     6,     7,     8,    20,    57,    58,    60,
      66,    16,    16,     9,    10,    11,    12,    13,    14,    15,
      16,    31,    32,    39,    78,    79,    80,    82,    83,    84,
      80,    58,     0,     3,    60,    16,    20,    33,    34,    74,
      75,    78,    11,    12,    15,    16,    32,    70,    71,    72,
      82,    82,    14,    16,    79,    78,    78,    14,    31,    32,
      36,    37,    23,    16,    18,    81,    21,    75,    78,    22,
      42,    43,    44,    45,    46,    47,    48,    49,    50,    73,
      74,    16,    16,    36,    37,    14,    70,    72,    15,    82,
      26,    27,    28,    29,    30,    31,    32,    17,    51,    61,
      17,    17,    79,    79,    79,    79,    78,    79,    80,    16,
      70,    18,    21,    22,    16,    32,    35,    77,    78,    84,
      77,    22,    70,    70,    15,    82,    15,    17,    17,    36,
      70,     9,    68,    70,    10,    69,    70,    68,    69,    70,
      70,    65,    72,    42,    17,    70,    19,    70,    77,    84,
      77,    22,    26,    27,    28,    29,    30,    31,    32,    36,
      37,    38,    39,    40,    41,    53,    76,    22,    23,    23,
      82,    16,    16,    20,    59,    60,    69,    19,    17,    77,
      77,    71,    71,    68,    69,    58,    62,    54,    17,    17,
      23,    23,    21,    22,    77,    68,    69,    72,    17,    17,
      63,    22,    33,    67,    82,    82,    17,    33,    42,    43,
      64,    82,    15,    59,    31,    15
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}

/* Prevent warnings from -Wmissing-prototypes.  */
#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */


/* The lookahead symbol.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*-------------------------.
| yyparse or yypush_parse.  |
`-------------------------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{


    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yytoken = 0;
  yyss = yyssa;
  yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */
  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:

/* Line 1455 of yacc.c  */
#line 141 "parser.y"
    {
	/* The full program was parsed. Allocate and fill the final
	   .scop structures. */
	int nb_parameters, nb_arrays;

        parser_scop->parameters = clan_symbol_id_array(parser_symbol,
                                                       SCOPLIB_TYPE_PARAMETER,
                                                       &nb_parameters);
        parser_scop->nb_parameters = nb_parameters;
        parser_scop->arrays = clan_symbol_id_array(parser_symbol,
                                                   SCOPLIB_TYPE_ARRAY,
                                                   &nb_arrays);
        parser_scop->nb_arrays = nb_arrays;
	if (parser_options->bounded_context)
	  {
	    parser_scop->context = scoplib_matrix_malloc(nb_parameters,
							 nb_parameters+2);
	    int i;
	    for (i = 0; i < nb_parameters; ++i)
	      {
		SCOPVAL_set_si(parser_scop->context->p[i][0], 1);
		SCOPVAL_set_si(parser_scop->context->p[i][i+1], 1);
		SCOPVAL_set_si(parser_scop->context->p[i][nb_parameters +1], 1);
	      }
	  }
	else
	  parser_scop->context = scoplib_matrix_malloc(0,nb_parameters+2);
      }
    break;

  case 11:

/* Line 1455 of yacc.c  */
#line 213 "parser.y"
    {
        clan_symbol_p symbol;
        symbol = clan_symbol_add(&parser_symbol,(yyvsp[(3) - (3)].symbol),
                                 SCOPLIB_TYPE_ITERATOR,parser_depth+1);
	/* Ensure that the returned symbol was either a new one,
	   either from the same type. */
	if (symbol->type != SCOPLIB_TYPE_ITERATOR)
	  {
	    yyerror("[Clan] Error: the input file is not a SCoP\n"
		    "\t> A loop iterator was previously used as a parameter"
		    "\n");
	    return 0;
	  }
	/* Update the rank, in case a symbol with the same name was
	   already existing. */
	if (symbol->rank != parser_depth + 1)
	  symbol->rank = parser_depth + 1;
        parser_iterators[parser_depth] = symbol;
	/* Memorize the current iterator as a negative constraint prefix */
      }
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 235 "parser.y"
    {
        scoplib_vector_p parser_i_term = clan_vector_term(parser_symbol,1,(yyvsp[(3) - (6)].symbol));
	scoplib_vector_tag_inequality(parser_i_term);
	int i, j;
	for (i = 0; i < (yyvsp[(6) - (6)].setex)->NbRows; ++i)
	  {
	    for (j = 1; j < (yyvsp[(6) - (6)].setex)->NbColumns; ++j)
	      SCOPVAL_oppose((yyvsp[(6) - (6)].setex)->p[i][j],(yyvsp[(6) - (6)].setex)->p[i][j]);
	    scoplib_matrix_add_vector((yyvsp[(6) - (6)].setex),parser_i_term,i);
	  }
	scoplib_matrix_insert_matrix(parser_domain,(yyvsp[(6) - (6)].setex),parser_nb_cons);

        parser_nb_cons += (yyvsp[(6) - (6)].setex)->NbRows;
        parser_consperdim[parser_depth] += (yyvsp[(6) - (6)].setex)->NbRows;
	scoplib_vector_free(parser_i_term);
        free((yyvsp[(3) - (6)].symbol));
	scoplib_matrix_free((yyvsp[(6) - (6)].setex));
      }
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 255 "parser.y"
    {
	scoplib_matrix_insert_matrix(parser_domain,(yyvsp[(9) - (9)].setex),parser_nb_cons);
        parser_nb_cons += (yyvsp[(9) - (9)].setex)->NbRows;
        parser_consperdim[parser_depth] += (yyvsp[(9) - (9)].setex)->NbRows;
      }
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 263 "parser.y"
    {
        parser_depth++;
        parser_scheduling[parser_depth] = 0;
      }
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 268 "parser.y"
    {
        parser_depth--;
        parser_scheduling[parser_depth]++;
        parser_nb_cons -= parser_consperdim[parser_depth];
        parser_consperdim[parser_depth] = 0;
	clan_symbol_remove(&parser_symbol, parser_iterators[parser_depth]);
      }
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 280 "parser.y"
    {
	/* Insert the condition constraint in the current parser domain. */
	scoplib_matrix_insert_matrix(parser_domain,(yyvsp[(3) - (4)].setex),parser_nb_cons);
        parser_nb_cons += (yyvsp[(3) - (4)].setex)->NbRows;
      }
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 286 "parser.y"
    {
        parser_nb_cons -= (yyvsp[(3) - (6)].setex)->NbRows;
	/* Remove the condition constraint from the current parser domain. */
	int i, j;
	for (i = parser_nb_cons; i < parser_domain->NbRows - 1; ++i)
	  for (j = 0; j < parser_domain->NbColumns; ++j)
	    SCOPVAL_assign(parser_domain->p[i][j],parser_domain->p[i+1][j]);
      }
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 298 "parser.y"
    {
        parser_statement = scoplib_statement_malloc();
        parser_record = (char *)malloc(SCOPLIB_MAX_STRING * sizeof(char));
        parser_recording = CLAN_TRUE;
        /* Yacc needs Lex to read the next token to ensure we are starting
         * an assignment. So we keep track of the latest text Lex read
         * and we start the statement body with it.
         */
        strcpy(parser_record,scanner_latest_text);
      }
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 309 "parser.y"
    {
	/* Deal with statements without surrounding loop by adding a
	   fake iterator */
	int old_parser_depth = parser_depth;
	if (parser_depth == 0)
	  {
	    char* fakeiter = strdup("fakeiter");
	    clan_symbol_p symbol = clan_symbol_lookup(parser_symbol, fakeiter);
	    if (symbol)
	      free(fakeiter);
	    else
	      symbol = clan_symbol_add(&parser_symbol,fakeiter,
				       SCOPLIB_TYPE_ITERATOR,parser_depth+1);
	    parser_iterators[parser_depth] = symbol;
	    scoplib_vector_p constraint =
	      scoplib_vector_malloc(parser_domain->NbColumns);
	    SCOPVAL_set_si(constraint->p[1],1);
	    parser_depth++;
	    scoplib_matrix_replace_vector(parser_domain,constraint,parser_nb_cons);
	    parser_nb_cons++;
	    scoplib_vector_free(constraint);
	  }
	/* Construct the statement structure from the parser state */
	parser_statement->domain = scoplib_matrix_list_malloc();
	parser_statement->domain->elt = scoplib_matrix_ncopy(parser_domain,
							  parser_nb_cons);
        parser_statement->schedule = clan_matrix_scheduling(parser_scheduling,
                                                            parser_depth);
        parser_statement->read = (yyvsp[(2) - (2)].rw)[0];
        parser_statement->write = (yyvsp[(2) - (2)].rw)[1];
        parser_statement->body = parser_record;
        parser_statement->nb_iterators = parser_depth;
        parser_statement->iterators = clan_symbol_iterators(parser_iterators,
                                                            parser_depth);
	if (parser_statement->write == NULL)
	  parser_statement->write =
	    scoplib_matrix_malloc(0, parser_domain->NbColumns);
	if (parser_statement->read == NULL)
	  parser_statement->read =
	    scoplib_matrix_malloc(0, parser_domain->NbColumns);
        parser_recording = CLAN_FALSE;
        scoplib_statement_add(&(parser_scop->statement),parser_statement);
	/* We were parsing a statement without iterator. Restore the
	   original state */
	if (old_parser_depth == 0)
	  {
	    --parser_depth;
	    --parser_nb_cons;
	    parser_consperdim[parser_depth] = 0;
	  }
        parser_scheduling[parser_depth]++;
      }
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 367 "parser.y"
    {
	int i, j;
	scoplib_matrix_p m = (yyvsp[(2) - (2)].setex);
	for (i = 0; i <  m->NbRows; ++i)
	  {
	    int id = SCOPVAL_get_si(m->p[i][0]);
	    for (j = 0; parser_variables_localvars[j] != -1 &&
		   parser_variables_localvars[j] != id; ++j)
	      ;
	    if (j == CLAN_MAX_LOCAL_VARIABLES)
	      {
		yyerror("[Clan] Error: maximum number of local variables reached\n");
		return 0;
	      }
	    if (parser_variables_localvars[j] == -1)
	      parser_variables_localvars[j] = id;
	  }
      }
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 391 "parser.y"
    {
	int i, j;
	scoplib_matrix_p m = (yyvsp[(2) - (2)].setex);
	for (i = 0; i <  m->NbRows; ++i)
	  {
	    int id = SCOPVAL_get_si(m->p[i][0]);
	    for (j = 0; parser_variables_liveout[j] != -1 &&
		   parser_variables_liveout[j] != id; ++j)
	      ;
	    if (j == CLAN_MAX_LOCAL_VARIABLES)
	      {
		yyerror("[Clan] Error: maximum number of live-out variables reached\n");
		return 0;
	      }
	    if (parser_variables_liveout[j] == -1)
	      parser_variables_liveout[j] = id;
	  }
      }
    break;

  case 22:

/* Line 1455 of yacc.c  */
#line 420 "parser.y"
    {
        free((yyvsp[(1) - (2)].symbol));
      }
    break;

  case 23:

/* Line 1455 of yacc.c  */
#line 424 "parser.y"
    {
        free((yyvsp[(2) - (2)].symbol));
      }
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 428 "parser.y"
    {
       if ((yyvsp[(5) - (5)].value) != 1)
	 {
	   yyerror("[Clan] Error: loop increment is not 1\n");
	   return 0;
	 }
       free ((yyvsp[(1) - (5)].symbol));
       free ((yyvsp[(3) - (5)].symbol));
     }
    break;

  case 25:

/* Line 1455 of yacc.c  */
#line 438 "parser.y"
    {
       if ((yyvsp[(3) - (3)].value) != 1)
	 {
	   yyerror("[Clan] Error: loop increment is not 1\n");
	   return 0;
	 }
       free ((yyvsp[(1) - (3)].symbol));
     }
    break;

  case 26:

/* Line 1455 of yacc.c  */
#line 456 "parser.y"
    {
	(yyval.setex) = scoplib_matrix_from_vector((yyvsp[(1) - (1)].affex));
        scoplib_vector_free((yyvsp[(1) - (1)].affex));
      }
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 462 "parser.y"
    {
       (yyval.setex) = scoplib_matrix_concat((yyvsp[(3) - (6)].setex), (yyvsp[(5) - (6)].setex));
     }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 475 "parser.y"
    {
	(yyval.setex) = scoplib_matrix_from_vector((yyvsp[(1) - (1)].affex));
        scoplib_vector_free((yyvsp[(1) - (1)].affex));
      }
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 481 "parser.y"
    {
       (yyval.setex) = scoplib_matrix_concat((yyvsp[(3) - (6)].setex), (yyvsp[(5) - (6)].setex));
     }
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 494 "parser.y"
    {
        (yyval.affex) = (yyvsp[(1) - (1)].affex);
      }
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 498 "parser.y"
    {
        (yyval.affex) = scoplib_vector_add((yyvsp[(1) - (3)].affex),(yyvsp[(3) - (3)].affex));
        scoplib_vector_free((yyvsp[(1) - (3)].affex));
        scoplib_vector_free((yyvsp[(3) - (3)].affex));
      }
    break;

  case 32:

/* Line 1455 of yacc.c  */
#line 504 "parser.y"
    {
        (yyval.affex) = scoplib_vector_sub((yyvsp[(1) - (3)].affex),(yyvsp[(3) - (3)].affex));
	scoplib_vector_free((yyvsp[(1) - (3)].affex));
        scoplib_vector_free((yyvsp[(3) - (3)].affex));
      }
    break;

  case 33:

/* Line 1455 of yacc.c  */
#line 510 "parser.y"
    {
        (yyval.affex) = (yyvsp[(2) - (3)].affex);
      }
    break;

  case 34:

/* Line 1455 of yacc.c  */
#line 514 "parser.y"
    {
	SCOPVAL_assign((yyvsp[(3) - (6)].affex)->p[0], (yyvsp[(5) - (6)].affex)->p[(yyvsp[(5) - (6)].affex)->Size - 1]);
	(yyval.affex) = (yyvsp[(3) - (6)].affex);
      }
    break;

  case 35:

/* Line 1455 of yacc.c  */
#line 519 "parser.y"
    {
	SCOPVAL_assign((yyvsp[(3) - (6)].affex)->p[0], (yyvsp[(5) - (6)].affex)->p[(yyvsp[(5) - (6)].affex)->Size - 1]);
	(yyval.affex) = (yyvsp[(3) - (6)].affex);
      }
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 536 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_symbol,(yyvsp[(1) - (1)].value),NULL);
      }
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 543 "parser.y"
    {
        clan_symbol_add(&parser_symbol,(yyvsp[(1) - (1)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_symbol,1,(yyvsp[(1) - (1)].symbol));
        free((yyvsp[(1) - (1)].symbol));
      }
    break;

  case 38:

/* Line 1455 of yacc.c  */
#line 552 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_symbol,-((yyvsp[(2) - (2)].value)),NULL);
      }
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 559 "parser.y"
    {
        clan_symbol_add(&parser_symbol,(yyvsp[(3) - (3)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_symbol,(yyvsp[(1) - (3)].value),(yyvsp[(3) - (3)].symbol));
        free((yyvsp[(3) - (3)].symbol));
      }
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 568 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_symbol, ((yyvsp[(1) - (3)].value)) * ((yyvsp[(3) - (3)].value)), NULL);
      }
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 575 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_symbol, ((yyvsp[(1) - (3)].value)) / ((yyvsp[(3) - (3)].value)), NULL);
      }
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 582 "parser.y"
    {
        clan_symbol_add(&parser_symbol,(yyvsp[(4) - (4)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_symbol,-((yyvsp[(2) - (4)].value)),(yyvsp[(4) - (4)].symbol));
        free((yyvsp[(4) - (4)].symbol));
      }
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 591 "parser.y"
    {
        clan_symbol_add(&parser_symbol,(yyvsp[(2) - (2)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_symbol,-1,(yyvsp[(2) - (2)].symbol));
        free((yyvsp[(2) - (2)].symbol));
      }
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 611 "parser.y"
    {
        /* a<b translates to -a+b-1>=0 */
	int i;
	scoplib_vector_p tmp = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),1);
	scoplib_vector_tag_inequality(tmp);
	for (i = 0; i < (yyvsp[(3) - (3)].setex)->NbRows; ++i)
	  {
	    /* We have parsed a ceild/floord at an earlier stage. */
	    if (SCOPVAL_notzero_p((yyvsp[(3) - (3)].setex)->p[i][0]) && !SCOPVAL_one_p((yyvsp[(3) - (3)].setex)->p[i][0]))
	      {
		scoplib_int_t val; SCOPVAL_init(val);
		SCOPVAL_assign(val, (yyvsp[(3) - (3)].setex)->p[i][0]);
		SCOPVAL_set_si((yyvsp[(3) - (3)].setex)->p[i][0], 0);
		scoplib_vector_p tmp2 = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),0);
		int j;
		for (j = 1; j < (yyvsp[(1) - (3)].affex)->Size; ++j)
		  SCOPVAL_multo(tmp2->p[j], (yyvsp[(1) - (3)].affex)->p[j], val);
		scoplib_vector_p tmp3 = scoplib_vector_add_scalar(tmp2,1);
		scoplib_vector_tag_inequality(tmp3);
		scoplib_matrix_sub_vector((yyvsp[(3) - (3)].setex), tmp3, i);
		scoplib_vector_free(tmp2);
		scoplib_vector_free(tmp3);
		SCOPVAL_clear(val);
	      }
	    else
	      scoplib_matrix_sub_vector((yyvsp[(3) - (3)].setex), tmp, i);
	  }
	scoplib_vector_free((yyvsp[(1) - (3)].affex));
	scoplib_vector_free(tmp);
	(yyval.setex) = (yyvsp[(3) - (3)].setex);
      }
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 646 "parser.y"
    {
        /* a>b translates to a-b-1>=0 */
	int i, j;
	scoplib_vector_p tmp = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),-1);
	scoplib_vector_tag_inequality(tmp);
	for (i = 0; i < (yyvsp[(3) - (3)].setex)->NbRows; ++i)
	  {
	    for (j = 1; j < (yyvsp[(3) - (3)].setex)->NbColumns; ++j)
	      SCOPVAL_oppose((yyvsp[(3) - (3)].setex)->p[i][j],(yyvsp[(3) - (3)].setex)->p[i][j]);
	    /* We have parsed a ceild/floord at an earlier stage. */
	    if (SCOPVAL_notzero_p((yyvsp[(3) - (3)].setex)->p[i][0]) && !SCOPVAL_one_p((yyvsp[(3) - (3)].setex)->p[i][0]))
	      {
		scoplib_int_t val; SCOPVAL_init(val);
		SCOPVAL_assign(val, (yyvsp[(3) - (3)].setex)->p[i][0]);
		SCOPVAL_set_si((yyvsp[(3) - (3)].setex)->p[i][0], 0);
		scoplib_vector_p tmp2 = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),0);
		int j;
		for (j = 1; j < (yyvsp[(1) - (3)].affex)->Size; ++j)
		  SCOPVAL_multo(tmp2->p[j], (yyvsp[(1) - (3)].affex)->p[j], val);
		scoplib_vector_p tmp3 = scoplib_vector_add_scalar(tmp2,-1);
		scoplib_vector_tag_inequality(tmp3);
		scoplib_matrix_add_vector((yyvsp[(3) - (3)].setex), tmp3, i);
		scoplib_vector_free(tmp2);
		scoplib_vector_free(tmp3);
		SCOPVAL_clear(val);
	      }
	    else
	      scoplib_matrix_add_vector((yyvsp[(3) - (3)].setex),tmp,i);
	  }
	scoplib_vector_free((yyvsp[(1) - (3)].affex));
	scoplib_vector_free(tmp);
	(yyval.setex) = (yyvsp[(3) - (3)].setex);
      }
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 683 "parser.y"
    {
        /* a<=b translates to -a+b>=0 */
	int i;
	scoplib_vector_p tmp = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),0);
	scoplib_vector_tag_inequality(tmp);
	for (i = 0; i < (yyvsp[(3) - (3)].setex)->NbRows; ++i)
	  {
	    /* We have parsed a ceild/floord at an earlier stage. */
	    if (SCOPVAL_notzero_p((yyvsp[(3) - (3)].setex)->p[i][0]) && !SCOPVAL_one_p((yyvsp[(3) - (3)].setex)->p[i][0]))
	      {
		scoplib_int_t val; SCOPVAL_init(val);
		SCOPVAL_assign(val, (yyvsp[(3) - (3)].setex)->p[i][0]);
		SCOPVAL_set_si((yyvsp[(3) - (3)].setex)->p[i][0], 0);
		scoplib_vector_p tmp2 = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),0);
		int j;
		for (j = 1; j < (yyvsp[(1) - (3)].affex)->Size; ++j)
		  SCOPVAL_multo(tmp2->p[j], (yyvsp[(1) - (3)].affex)->p[j], val);
		scoplib_vector_tag_inequality(tmp2);
		scoplib_matrix_sub_vector((yyvsp[(3) - (3)].setex), tmp2, i);
		scoplib_vector_free(tmp2);
		SCOPVAL_clear(val);
	      }
	    else
	      scoplib_matrix_sub_vector((yyvsp[(3) - (3)].setex),tmp,i);
	  }
	scoplib_vector_free((yyvsp[(1) - (3)].affex));
	scoplib_vector_free(tmp);
	(yyval.setex) = (yyvsp[(3) - (3)].setex);
      }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 716 "parser.y"
    {
        /* a>=b translates to a-b>=0 */
	int i, j;
	scoplib_vector_p tmp = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),0);
	scoplib_vector_tag_inequality(tmp);
	for (i = 0; i < (yyvsp[(3) - (3)].setex)->NbRows; ++i)
	  {
	    for (j = 1; j < (yyvsp[(3) - (3)].setex)->NbColumns; ++j)
	      SCOPVAL_oppose((yyvsp[(3) - (3)].setex)->p[i][j],(yyvsp[(3) - (3)].setex)->p[i][j]);
	    /* We have parsed a ceild/floord at an earlier stage. */
	    if (SCOPVAL_notzero_p((yyvsp[(3) - (3)].setex)->p[i][0]) && !SCOPVAL_one_p((yyvsp[(3) - (3)].setex)->p[i][0]))
	      {
		scoplib_int_t val; SCOPVAL_init(val);
		SCOPVAL_assign(val, (yyvsp[(3) - (3)].setex)->p[i][0]);
		SCOPVAL_set_si((yyvsp[(3) - (3)].setex)->p[i][0], 0);
		scoplib_vector_p tmp2 = scoplib_vector_add_scalar((yyvsp[(1) - (3)].affex),0);
		int j;
		for (j = 1; j < (yyvsp[(1) - (3)].affex)->Size; ++j)
		  SCOPVAL_multo(tmp2->p[j], (yyvsp[(1) - (3)].affex)->p[j], val);
		scoplib_vector_tag_inequality(tmp2);
		scoplib_matrix_add_vector((yyvsp[(3) - (3)].setex), tmp2, i);
		scoplib_vector_free(tmp2);
		SCOPVAL_clear(val);
	      }
	    else
	      scoplib_matrix_add_vector((yyvsp[(3) - (3)].setex),tmp,i);
	  }
	scoplib_vector_free((yyvsp[(1) - (3)].affex));
	scoplib_vector_free(tmp);
	(yyval.setex) = (yyvsp[(3) - (3)].setex);
      }
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 751 "parser.y"
    {
        /* a==b translates to a-b==0 */
	/* Warning: cases like ceild(M,32) == ceild(N,32) are not handled.
	   Assert if we encounter such a case. */
	assert ((SCOPVAL_zero_p((yyvsp[(1) - (3)].affex)->p[0]) || SCOPVAL_one_p((yyvsp[(1) - (3)].affex)->p[0]))
		&& (SCOPVAL_zero_p((yyvsp[(3) - (3)].affex)->p[0]) || SCOPVAL_one_p((yyvsp[(3) - (3)].affex)->p[0])));
	scoplib_vector_p res = scoplib_vector_sub((yyvsp[(1) - (3)].affex),(yyvsp[(3) - (3)].affex));
	scoplib_vector_tag_equality(res);
	(yyval.setex) = scoplib_matrix_from_vector(res);
	scoplib_vector_free(res);
        scoplib_vector_free((yyvsp[(1) - (3)].affex));
	scoplib_vector_free((yyvsp[(3) - (3)].affex));
      }
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 768 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (3)].setex);
      }
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 775 "parser.y"
    {
       (yyval.setex) = scoplib_matrix_concat((yyvsp[(1) - (3)].setex),(yyvsp[(3) - (3)].setex));
       scoplib_matrix_free((yyvsp[(1) - (3)].setex));
       scoplib_matrix_free((yyvsp[(3) - (3)].setex));
     }
    break;

  case 61:

/* Line 1455 of yacc.c  */
#line 819 "parser.y"
    {
	if ((yyvsp[(1) - (4)].setex) == NULL)
	  {
	    yyerror ("[Clan] Error: changing value of iterator/parameter");
	    return 0;
	  }
        (yyval.rw)[0] = (yyvsp[(3) - (4)].setex);
        (yyval.rw)[1] = (yyvsp[(1) - (4)].setex);
      }
    break;

  case 62:

/* Line 1455 of yacc.c  */
#line 832 "parser.y"
    {
	if ((yyvsp[(1) - (4)].setex) == NULL)
	  {
	    yyerror ("[Clan] Error: changing value of iterator/parameter");
	    return 0;
	  }
        (yyval.rw)[0] = scoplib_matrix_concat((yyvsp[(1) - (4)].setex),(yyvsp[(3) - (4)].setex));
        scoplib_matrix_free((yyvsp[(3) - (4)].setex));
        (yyval.rw)[1] = (yyvsp[(1) - (4)].setex);
      }
    break;

  case 63:

/* Line 1455 of yacc.c  */
#line 846 "parser.y"
    {
	if ((yyvsp[(1) - (3)].setex) == NULL)
	  {
	    yyerror ("[Clan] Error: changing value of iterator/parameter");
	    return 0;
	  }
        (yyval.rw)[0] = (yyvsp[(1) - (3)].setex);
        (yyval.rw)[1] = scoplib_matrix_copy((yyvsp[(1) - (3)].setex));
      }
    break;

  case 64:

/* Line 1455 of yacc.c  */
#line 859 "parser.y"
    {
	if ((yyvsp[(2) - (3)].setex) == NULL)
	  {
	    yyerror ("[Clan] Error: changing value of iterator/parameter");
	    return 0;
	  }
       (yyval.rw)[0] = (yyvsp[(2) - (3)].setex);
       (yyval.rw)[1] = scoplib_matrix_copy((yyvsp[(2) - (3)].setex));
      }
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 872 "parser.y"
    {
       (yyval.rw)[0] = (yyvsp[(1) - (2)].setex);
       (yyval.rw)[1] = NULL;
     }
    break;

  case 66:

/* Line 1455 of yacc.c  */
#line 880 "parser.y"
    {
       (yyval.rw)[0] = (yyvsp[(2) - (3)].rw)[0];
       (yyval.rw)[1] = (yyvsp[(2) - (3)].rw)[1];
     }
    break;

  case 80:

/* Line 1455 of yacc.c  */
#line 916 "parser.y"
    {
        (yyval.setex) = NULL;
      }
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 923 "parser.y"
    {
        (yyval.setex) = NULL;
      }
    break;

  case 82:

/* Line 1455 of yacc.c  */
#line 930 "parser.y"
    {
        (yyval.setex) = (yyvsp[(1) - (1)].setex);
      }
    break;

  case 83:

/* Line 1455 of yacc.c  */
#line 938 "parser.y"
    {
        (yyval.setex) = scoplib_matrix_concat((yyvsp[(1) - (3)].setex),(yyvsp[(3) - (3)].setex));
	scoplib_matrix_free((yyvsp[(1) - (3)].setex));
        scoplib_matrix_free((yyvsp[(3) - (3)].setex));
      }
    break;

  case 84:

/* Line 1455 of yacc.c  */
#line 947 "parser.y"
    {
        (yyval.setex) = (yyvsp[(2) - (2)].setex);
      }
    break;

  case 85:

/* Line 1455 of yacc.c  */
#line 954 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (3)].setex);
      }
    break;

  case 86:

/* Line 1455 of yacc.c  */
#line 961 "parser.y"
    {
	scoplib_matrix_p tmp = scoplib_matrix_concat((yyvsp[(1) - (5)].setex),(yyvsp[(3) - (5)].setex));
        (yyval.setex) = scoplib_matrix_concat(tmp,(yyvsp[(5) - (5)].setex));
	scoplib_matrix_free(tmp);
	scoplib_matrix_free((yyvsp[(1) - (5)].setex));
	scoplib_matrix_free((yyvsp[(3) - (5)].setex));
	scoplib_matrix_free((yyvsp[(5) - (5)].setex));
      }
    break;

  case 87:

/* Line 1455 of yacc.c  */
#line 983 "parser.y"
    {
        int rank;
        scoplib_matrix_p matrix;
	char* s = (char*) (yyvsp[(1) - (1)].symbol);
	clan_symbol_p symbol = clan_symbol_lookup(parser_symbol, s);
	// If the variable is an iterator or a parameter, discard it
	// from the read/write clause.
	if ((symbol && symbol->type == SCOPLIB_TYPE_ITERATOR) ||
	     (symbol && symbol->type == SCOPLIB_TYPE_PARAMETER))
	  (yyval.setex) = NULL;
	else
	  {
	    clan_symbol_add(&parser_symbol, s, SCOPLIB_TYPE_ARRAY,parser_depth);
	    rank = clan_symbol_get_rank(parser_symbol, s);
	    matrix = scoplib_matrix_malloc
	      (1, CLAN_MAX_DEPTH + CLAN_MAX_PARAMETERS + 2);
	    clan_matrix_tag_array(matrix, rank);
	    (yyval.setex) = matrix;
	  }
        free((yyvsp[(1) - (1)].symbol));
      }
    break;

  case 88:

/* Line 1455 of yacc.c  */
#line 1009 "parser.y"
    {
        int rank;
        clan_symbol_add(&parser_symbol,(yyvsp[(1) - (2)].symbol),SCOPLIB_TYPE_ARRAY,parser_depth);
        rank = clan_symbol_get_rank(parser_symbol,(yyvsp[(1) - (2)].symbol));
        clan_matrix_tag_array((yyvsp[(2) - (2)].setex),rank);
        (yyval.setex) = (yyvsp[(2) - (2)].setex);
        free((yyvsp[(1) - (2)].symbol));
      }
    break;

  case 89:

/* Line 1455 of yacc.c  */
#line 1022 "parser.y"
    {
	(yyval.setex) = (yyvsp[(3) - (4)].setex);
	free((yyvsp[(1) - (4)].symbol));
      }
    break;

  case 90:

/* Line 1455 of yacc.c  */
#line 1030 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (2)].setex);
      }
    break;

  case 91:

/* Line 1455 of yacc.c  */
#line 1037 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (2)].setex);
      }
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 1065 "parser.y"
    {
	(yyval.setex) = (yyvsp[(1) - (1)].setex);
      }
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 1072 "parser.y"
    {
	(yyval.setex) = scoplib_matrix_concat((yyvsp[(1) - (3)].setex),(yyvsp[(3) - (3)].setex));
      }
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 1079 "parser.y"
    {
	(yyval.setex) = (yyvsp[(1) - (3)].setex);
      }
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 1086 "parser.y"
    {
	(yyval.setex) = NULL;
      }
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 1093 "parser.y"
    {
	(yyval.setex) = NULL;
      }
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 1109 "parser.y"
    {
        (yyval.setex) = scoplib_matrix_from_vector((yyvsp[(2) - (3)].affex));
        scoplib_vector_free((yyvsp[(2) - (3)].affex));
      }
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 1117 "parser.y"
    {
	if ((yyvsp[(1) - (4)].setex) != NULL)
	  scoplib_matrix_insert_vector((yyvsp[(1) - (4)].setex),(yyvsp[(3) - (4)].affex),(yyvsp[(1) - (4)].setex)->NbRows);
        scoplib_vector_free((yyvsp[(3) - (4)].affex));
        (yyval.setex) = (yyvsp[(1) - (4)].setex);
      }
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 1138 "parser.y"
    {
       (yyval.symbol) = (yyvsp[(1) - (1)].symbol);
     }
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 1145 "parser.y"
    {
       (yyval.symbol) = (yyvsp[(2) - (3)].symbol);
     }
    break;

  case 107:

/* Line 1455 of yacc.c  */
#line 1152 "parser.y"
    {
       (yyval.symbol) = (yyvsp[(2) - (2)].symbol);
     }
    break;

  case 108:

/* Line 1455 of yacc.c  */
#line 1159 "parser.y"
    {
       (yyval.symbol) = NULL;
     }
    break;



/* Line 1455 of yacc.c  */
#line 2752 "parser.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#if !defined(yyoverflow) || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}



/* Line 1675 of yacc.c  */
#line 1171 "parser.y"



void
yyerror(char *s)
{
  fprintf(stderr, "%s\n", s);
  clan_parse_error = 1;
}


/**
 * clan_parser_initialize_state function:
 * this function achieves the initialization of the "parser state": a
 * collection of variables that vary during the parsing and thanks to we
 * can extract all SCoP informations.
 **
 * - 02/05/2008: First version.
 */
void
clan_parser_initialize_state(clan_options_p options)
{
  int i, nb_rows, nb_columns, depth;

  nb_rows    = CLAN_MAX_CONSTRAINTS;
  nb_columns = CLAN_MAX_DEPTH + CLAN_MAX_PARAMETERS + 2;
  depth      = CLAN_MAX_DEPTH;

  parser_scop   = scoplib_scop_malloc();
  parser_domain = scoplib_matrix_malloc(nb_rows,nb_columns);
  parser_symbol = NULL;

  parser_scheduling = (int *)malloc(depth * sizeof(int));
  parser_consperdim = (int *)malloc(depth * sizeof(int));
  for (i = 0; i < depth; i++)
  {
    parser_scheduling[i] = 0;
    parser_consperdim[i] = 0;
  }
  parser_iterators = (clan_symbol_p *)malloc(depth * sizeof(clan_symbol_p));
  parser_variables_localvars =
    (int*)malloc((CLAN_MAX_LOCAL_VARIABLES + 1) * sizeof(int));
  parser_variables_liveout =
    (int*)malloc((CLAN_MAX_LOCAL_VARIABLES + 1) * sizeof(int));
  parser_depth = 0;
  parser_nb_cons = 0;
  /* Reset also the Symbol global variables. */
  extern int symbol_nb_iterators;
  symbol_nb_iterators = 0;
  extern int symbol_nb_parameters;
  symbol_nb_parameters = 0;
  extern int symbol_nb_arrays;
  symbol_nb_arrays = 0;
  extern int symbol_nb_functions;
  symbol_nb_functions = 0;

  for (i = 0; i <= CLAN_MAX_LOCAL_VARIABLES; ++i)
    parser_variables_localvars[i] = -1;
  for (i = 0; i <= CLAN_MAX_LOCAL_VARIABLES; ++i)
    parser_variables_liveout[i] = -1;

  parser_options = options;
}

/**
 * clan_parser_free_state function:
 * this function frees the memory allocated for the "parser state", except
 * for parser_scop, obviously.
 **
 * - 02/05/2008: First version.
 */
void
clan_parser_free_state()
{
  scoplib_matrix_free(parser_domain);
  clan_symbol_free(parser_symbol);
  free(parser_scheduling);
  free(parser_consperdim);
  free(parser_iterators);
  free(parser_variables_localvars);
  free(parser_variables_liveout);
}

/**
 * clan_parse function:
 * this function parses a file to extract a SCoP and returns, if successful,
 * a pointer to the scoplib_scop_t structure.
 * \param input   The file to parse (already open).
 * \param options Options for file parsing.
 **
 * - 01/05/2008: First version.
 */
scoplib_scop_p
clan_parse(FILE * input, clan_options_p options)
{
  yyin = input;

  clan_parser_initialize_state(options);

  yyparse();

  fclose(yyin);
  if (! clan_parse_error)
    {
      if (parser_variables_localvars[0] != -1 ||
	  parser_variables_liveout[0] != -1)
	clan_scop_fill_options(parser_scop, parser_variables_localvars,
			       parser_variables_liveout);
      clan_scop_compact(parser_scop);
    }
  else
    parser_scop = NULL;
  clan_parser_free_state();

  return parser_scop;
}


