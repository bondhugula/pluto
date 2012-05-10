
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
   #include <scoplib/symbol.h>

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
   clan_symbol_p       parser_clan_symbol;      /**< Top of the symbol table */
   scoplib_symbol_p    parser_scop_symbol; /**<Top of the symbol table */
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
   
   char* parser_symbol_datatype;
   
    struct pointer {
            int num_of_references;
    };
            
    struct datatype {
            char* datatype_name;
            struct pointer pointer_details;
            int pre_defined_type;
    };     




/* Line 189 of yacc.c  */
#line 139 "parser.c"

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
     CLAN_PARSER_DECL_START = 258,
     CLAN_PARSER_DECL_END = 259,
     CLAN_PARSER_SCOP_START = 260,
     CLAN_PARSER_SCOP_END = 261,
     IGNORE = 262,
     STRUCT = 263,
     IF = 264,
     ELSE = 265,
     FOR = 266,
     PRAGMALOCALVARS = 267,
     PRAGMALIVEOUT = 268,
     MIN = 269,
     MAX = 270,
     CEILD = 271,
     FLOORD = 272,
     REAL = 273,
     ID = 274,
     INTEGER = 275,
     DATATYPE = 276,
     syRPARENTHESIS = 277,
     syLPARENTHESIS = 278,
     syRBRACKET = 279,
     syLBRACKET = 280,
     syRBRACE = 281,
     syLBRACE = 282,
     sySEMICOLON = 283,
     syCOMMA = 284,
     syPOINT = 285,
     syARROW = 286,
     opEQUAL = 287,
     opLEQ = 288,
     opGEQ = 289,
     opLOWER = 290,
     opGREATER = 291,
     opPLUS = 292,
     opMINUS = 293,
     opINCREMENTATION = 294,
     opDECREMENTATION = 295,
     opNOT = 296,
     opMULTIPLY = 297,
     opDIVIDE = 298,
     opMOD = 299,
     opAND = 300,
     opOR = 301,
     opCOMP = 302,
     opASSIGNMENT = 303,
     opPLUSEQUAL = 304,
     opMINUSEQUAL = 305,
     opMULTIPLYEQUAL = 306,
     opDIVIDEEQUAL = 307,
     opMODEQUAL = 308,
     opANDEQUAL = 309,
     opOREQUAL = 310,
     opCOMPEQUAL = 311,
     opLAND = 312,
     opLOR = 313,
     opQMARK = 314,
     opCOLON = 315,
     MAXPRIORITY = 316
   };
#endif
/* Tokens.  */
#define CLAN_PARSER_DECL_START 258
#define CLAN_PARSER_DECL_END 259
#define CLAN_PARSER_SCOP_START 260
#define CLAN_PARSER_SCOP_END 261
#define IGNORE 262
#define STRUCT 263
#define IF 264
#define ELSE 265
#define FOR 266
#define PRAGMALOCALVARS 267
#define PRAGMALIVEOUT 268
#define MIN 269
#define MAX 270
#define CEILD 271
#define FLOORD 272
#define REAL 273
#define ID 274
#define INTEGER 275
#define DATATYPE 276
#define syRPARENTHESIS 277
#define syLPARENTHESIS 278
#define syRBRACKET 279
#define syLBRACKET 280
#define syRBRACE 281
#define syLBRACE 282
#define sySEMICOLON 283
#define syCOMMA 284
#define syPOINT 285
#define syARROW 286
#define opEQUAL 287
#define opLEQ 288
#define opGEQ 289
#define opLOWER 290
#define opGREATER 291
#define opPLUS 292
#define opMINUS 293
#define opINCREMENTATION 294
#define opDECREMENTATION 295
#define opNOT 296
#define opMULTIPLY 297
#define opDIVIDE 298
#define opMOD 299
#define opAND 300
#define opOR 301
#define opCOMP 302
#define opASSIGNMENT 303
#define opPLUSEQUAL 304
#define opMINUSEQUAL 305
#define opMULTIPLYEQUAL 306
#define opDIVIDEEQUAL 307
#define opMODEQUAL 308
#define opANDEQUAL 309
#define opOREQUAL 310
#define opCOMPEQUAL 311
#define opLAND 312
#define opLOR 313
#define opQMARK 314
#define opCOLON 315
#define MAXPRIORITY 316




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 214 of yacc.c  */
#line 104 "parser.y"
 int value;                     /**< An integer value for integers */
         char * symbol;                 /**< A string for identifiers */
         scoplib_vector_p affex;        /**< An affine expression */
         scoplib_matrix_p setex;        /**< A set of affine expressions */
         scoplib_matrix_p rw[2];        /**< Read and write array accesses */
         scoplib_symbol_p symbol_table;         
       


/* Line 214 of yacc.c  */
#line 307 "parser.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif


/* Copy the second part of user declarations.  */


/* Line 264 of yacc.c  */
#line 319 "parser.c"

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
#define YYFINAL  19
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   557

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  62
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  35
/* YYNRULES -- Number of rules.  */
#define YYNRULES  129
/* YYNRULES -- Number of states.  */
#define YYNSTATES  247

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   316

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
      55,    56,    57,    58,    59,    60,    61
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     8,    13,    16,    17,    19,    20,    22,
      25,    29,    32,    33,    35,    38,    45,    46,    50,    51,
      53,    56,    61,    63,    67,    68,    69,    70,    71,    87,
      88,    95,    96,    99,   102,   105,   108,   111,   117,   121,
     123,   130,   132,   139,   141,   145,   149,   153,   160,   167,
     169,   171,   174,   178,   182,   186,   190,   195,   200,   203,
     207,   211,   215,   219,   223,   227,   231,   233,   235,   237,
     239,   241,   243,   245,   247,   249,   251,   256,   261,   265,
     269,   272,   276,   278,   280,   282,   284,   286,   288,   290,
     292,   294,   296,   298,   300,   302,   304,   307,   309,   313,
     316,   320,   326,   328,   331,   336,   339,   342,   344,   348,
     352,   356,   360,   364,   366,   370,   374,   376,   377,   381,
     386,   388,   392,   395,   397,   399,   401,   403,   405,   407
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      63,     0,    -1,     3,    68,     4,    63,    -1,     5,    64,
       6,    63,    -1,     7,    63,    -1,    -1,    65,    -1,    -1,
      72,    -1,    65,    72,    -1,    26,    65,    27,    -1,    42,
      66,    -1,    -1,    21,    -1,     8,    94,    -1,    67,    66,
      70,    69,    28,    68,    -1,    -1,    29,    70,    69,    -1,
      -1,    94,    -1,    94,    93,    -1,    94,    22,    92,    23,
      -1,    72,    -1,    26,    65,    27,    -1,    -1,    -1,    -1,
      -1,    11,    22,    94,    73,    48,    81,    74,    28,    84,
      75,    28,    79,    23,    76,    71,    -1,    -1,     9,    22,
      84,    23,    77,    71,    -1,    -1,    78,    87,    -1,    12,
      92,    -1,    13,    92,    -1,    94,    39,    -1,    39,    94,
      -1,    94,    48,    94,    37,    20,    -1,    94,    49,    20,
      -1,    82,    -1,    14,    22,    80,    29,    80,    23,    -1,
      82,    -1,    15,    22,    81,    29,    81,    23,    -1,    83,
      -1,    82,    37,    82,    -1,    82,    38,    82,    -1,    22,
      82,    23,    -1,    16,    22,    82,    29,    83,    23,    -1,
      17,    22,    82,    29,    83,    23,    -1,    20,    -1,    94,
      -1,    38,    20,    -1,    20,    42,    94,    -1,    94,    42,
      20,    -1,    20,    42,    20,    -1,    20,    43,    20,    -1,
      38,    20,    42,    94,    -1,    38,    94,    42,    20,    -1,
      38,    94,    -1,    82,    35,    80,    -1,    82,    36,    81,
      -1,    82,    33,    80,    -1,    82,    34,    81,    -1,    82,
      32,    82,    -1,    22,    84,    23,    -1,    84,    57,    84,
      -1,    49,    -1,    50,    -1,    51,    -1,    52,    -1,    53,
      -1,    54,    -1,    55,    -1,    56,    -1,    39,    -1,    40,
      -1,    90,    48,    89,    28,    -1,    90,    85,    89,    28,
      -1,    90,    86,    28,    -1,    86,    90,    28,    -1,    90,
      28,    -1,    26,    87,    27,    -1,    37,    -1,    38,    -1,
      42,    -1,    43,    -1,    44,    -1,    34,    -1,    36,    -1,
      33,    -1,    35,    -1,    32,    -1,    45,    -1,    46,    -1,
      47,    -1,    96,    -1,    38,    96,    -1,    90,    -1,    89,
      88,    89,    -1,    41,    89,    -1,    22,    89,    23,    -1,
      89,    59,    89,    60,    89,    -1,    94,    -1,    94,    93,
      -1,    94,    22,    92,    23,    -1,    38,    90,    -1,    37,
      90,    -1,    96,    -1,    91,    38,    91,    -1,    91,    37,
      91,    -1,    91,    42,    91,    -1,    91,    43,    91,    -1,
      22,    91,    23,    -1,    90,    -1,    92,    29,    90,    -1,
      92,    29,    91,    -1,    91,    -1,    -1,    24,    82,    25,
      -1,    93,    24,    82,    25,    -1,    19,    -1,    22,    19,
      23,    -1,    45,    19,    -1,    95,    -1,    14,    -1,    15,
      -1,    16,    -1,    17,    -1,    20,    -1,    18,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   164,   164,   166,   168,   170,   176,   205,   213,   214,
     215,   220,   222,   228,   230,   236,   238,   244,   246,   261,
     272,   286,   298,   302,   319,   341,   361,   369,   316,   386,
     385,   404,   404,   472,   496,   526,   530,   534,   544,   562,
     567,   581,   586,   600,   604,   610,   616,   620,   625,   642,
     649,   658,   665,   674,   683,   690,   697,   706,   715,   735,
     770,   807,   840,   875,   892,   899,   913,   914,   915,   916,
     917,   918,   919,   920,   929,   930,   943,   956,   970,   983,
     996,  1004,  1017,  1018,  1019,  1020,  1021,  1022,  1023,  1024,
    1025,  1026,  1027,  1028,  1029,  1040,  1047,  1054,  1062,  1071,
    1078,  1085,  1107,  1133,  1146,  1154,  1161,  1172,  1173,  1174,
    1175,  1176,  1177,  1189,  1196,  1203,  1210,  1218,  1233,  1241,
    1262,  1269,  1276,  1283,  1289,  1289,  1289,  1289,  1292,  1293
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "CLAN_PARSER_DECL_START",
  "CLAN_PARSER_DECL_END", "CLAN_PARSER_SCOP_START", "CLAN_PARSER_SCOP_END",
  "IGNORE", "STRUCT", "IF", "ELSE", "FOR", "PRAGMALOCALVARS",
  "PRAGMALIVEOUT", "MIN", "MAX", "CEILD", "FLOORD", "REAL", "ID",
  "INTEGER", "DATATYPE", "syRPARENTHESIS", "syLPARENTHESIS", "syRBRACKET",
  "syLBRACKET", "syRBRACE", "syLBRACE", "sySEMICOLON", "syCOMMA",
  "syPOINT", "syARROW", "opEQUAL", "opLEQ", "opGEQ", "opLOWER",
  "opGREATER", "opPLUS", "opMINUS", "opINCREMENTATION", "opDECREMENTATION",
  "opNOT", "opMULTIPLY", "opDIVIDE", "opMOD", "opAND", "opOR", "opCOMP",
  "opASSIGNMENT", "opPLUSEQUAL", "opMINUSEQUAL", "opMULTIPLYEQUAL",
  "opDIVIDEEQUAL", "opMODEQUAL", "opANDEQUAL", "opOREQUAL", "opCOMPEQUAL",
  "opLAND", "opLOR", "opQMARK", "opCOLON", "MAXPRIORITY", "$accept",
  "program", "scop_instructions", "instruction_list", "pointer",
  "datatype", "declarations", "ctsDeclarations", "variableDeclaration",
  "bloc", "instruction", "$@1", "$@2", "$@3", "$@4", "$@5", "$@6",
  "incrementation", "min_affine_expression", "max_affine_expression",
  "affine_expression", "term", "condition", "reduction_operator",
  "unary_operator", "assignment", "binary_operator", "expression",
  "variable", "arithmetic_expression", "variable_list", "array_index",
  "id", "math_func_list", "NUMBER", 0
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
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    62,    63,    63,    63,    63,    64,    64,    65,    65,
      65,    66,    66,    67,    67,    68,    68,    69,    69,    70,
      70,    70,    71,    71,    73,    74,    75,    76,    72,    77,
      72,    78,    72,    72,    72,    79,    79,    79,    79,    80,
      80,    81,    81,    82,    82,    82,    82,    82,    82,    83,
      83,    83,    83,    83,    83,    83,    83,    83,    83,    84,
      84,    84,    84,    84,    84,    84,    85,    85,    85,    85,
      85,    85,    85,    85,    86,    86,    87,    87,    87,    87,
      87,    87,    88,    88,    88,    88,    88,    88,    88,    88,
      88,    88,    88,    88,    88,    89,    89,    89,    89,    89,
      89,    89,    90,    90,    90,    90,    90,    91,    91,    91,
      91,    91,    91,    92,    92,    92,    92,    92,    93,    93,
      94,    94,    94,    94,    95,    95,    95,    95,    96,    96
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     4,     4,     2,     0,     1,     0,     1,     2,
       3,     2,     0,     1,     2,     6,     0,     3,     0,     1,
       2,     4,     1,     3,     0,     0,     0,     0,    15,     0,
       6,     0,     2,     2,     2,     2,     2,     5,     3,     1,
       6,     1,     6,     1,     3,     3,     3,     6,     6,     1,
       1,     2,     3,     3,     3,     3,     4,     4,     2,     3,
       3,     3,     3,     3,     3,     3,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     4,     4,     3,     3,
       2,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     2,     1,     3,     2,
       3,     5,     1,     2,     4,     2,     2,     1,     3,     3,
       3,     3,     3,     1,     3,     3,     1,     0,     3,     4,
       1,     3,     2,     1,     1,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       5,    16,    31,     5,     0,     0,    13,    12,     0,     0,
       0,   117,   117,    31,     0,    31,     8,     0,     4,     1,
     124,   125,   126,   127,   120,     0,     0,    14,   123,    12,
       0,     5,     0,     0,   129,   128,     0,     0,     0,   113,
     116,    33,   102,   107,    34,    31,     5,     9,     0,    74,
      75,     0,    32,     0,     0,   122,    11,    18,    19,     2,
     126,   127,    49,     0,     0,     0,    43,     0,    50,    24,
       0,     0,   106,   105,     0,     0,     0,     0,     0,   117,
       0,   103,    10,     3,     0,     0,    80,     0,    66,    67,
      68,    69,    70,    71,    72,    73,     0,     0,   121,     0,
       0,   117,    20,     0,     0,     0,     0,   120,     0,     0,
      51,    58,     0,     0,     0,     0,     0,     0,     0,    29,
       0,     0,     0,   112,   109,   108,   110,   111,   114,   115,
       0,     0,     0,     0,    81,    79,     0,     0,     0,     0,
      97,    95,     0,    78,    18,    16,     0,     0,     0,    54,
      52,    55,    46,    64,     0,     0,    63,   124,    61,    39,
     125,    62,    41,    59,    60,    44,    45,    31,    65,    53,
       0,   104,     0,   118,     0,     0,    96,    99,    76,    91,
      89,    87,    90,    88,    82,    83,    84,    85,    86,    92,
      93,    94,     0,     0,    77,    17,    15,    21,     0,     0,
      56,    57,     0,     0,    31,    30,    22,    25,   119,   100,
       0,    98,     0,     0,     0,     0,    31,     0,     0,    47,
      48,     0,     0,    23,     0,   101,     0,     0,    26,    40,
      42,     0,     0,     0,     0,     0,    36,    27,    35,     0,
       0,    31,     0,    38,    28,     0,    37
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     4,    14,    15,    30,     7,     8,   100,    57,   205,
      16,   122,   217,   231,   241,   167,    17,   234,   158,   161,
     162,    66,    67,    96,    51,    52,   193,   139,   140,    40,
      41,    81,    68,    28,    43
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -111
static const yytype_int16 yypact[] =
{
     294,    10,   149,   294,    22,   371,  -111,   -15,    25,    54,
      79,   301,   301,    43,    32,   421,  -111,   180,  -111,  -111,
    -111,  -111,  -111,  -111,  -111,    46,    78,  -111,  -111,   -15,
     371,   294,   337,   371,  -111,  -111,   190,   346,   346,  -111,
     147,    84,    42,  -111,    84,    66,   294,  -111,   180,  -111,
    -111,   346,  -111,   491,    81,  -111,  -111,    87,    88,  -111,
      85,    97,   -32,   357,   451,   519,  -111,    -8,    86,  -111,
     291,    92,  -111,  -111,   291,   291,   291,   291,   301,   301,
     384,   119,  -111,  -111,    94,   109,  -111,   257,  -111,  -111,
    -111,  -111,  -111,  -111,  -111,  -111,   257,   131,  -111,   371,
     135,   301,   119,   384,   384,   467,   166,    81,   471,    -6,
     156,   158,   384,   393,   404,   393,   404,   384,   384,  -111,
     337,   171,   153,  -111,    47,    47,  -111,  -111,  -111,   147,
      18,   431,    20,   384,  -111,  -111,   269,   310,   257,   203,
    -111,  -111,   223,  -111,    87,    10,    38,   113,   136,  -111,
    -111,  -111,  -111,  -111,   371,   183,   145,   182,  -111,   145,
     189,  -111,   145,  -111,  -111,  -111,  -111,   324,  -111,  -111,
     404,  -111,   -14,  -111,   116,   134,  -111,  -111,  -111,  -111,
    -111,  -111,  -111,  -111,  -111,  -111,  -111,  -111,  -111,  -111,
    -111,  -111,   257,   257,  -111,  -111,  -111,  -111,   440,   440,
    -111,  -111,   393,   404,    43,  -111,  -111,  -111,  -111,  -111,
     490,  -111,   193,   198,   194,   195,    75,   204,   257,  -111,
    -111,   393,   404,  -111,   337,  -111,   199,   210,   185,  -111,
    -111,   215,   476,   371,   221,   -13,  -111,  -111,  -111,   371,
     232,   324,   216,  -111,  -111,   234,  -111
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -111,    17,  -111,   -12,   235,  -111,   118,   137,   179,    39,
     -11,  -111,  -111,  -111,  -111,  -111,  -111,  -111,  -110,  -108,
       5,    16,   -60,  -111,   237,   244,  -111,   -94,     2,    70,
       4,   238,    -5,  -111,   -66
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -8
static const yytype_int16 yytable[] =
{
      27,    45,   142,   109,    47,   163,    42,    42,   164,   152,
     105,   106,    42,    39,    39,   119,    44,   153,     5,    53,
      18,   141,    19,   117,   118,    58,   238,    29,    69,    31,
     141,     6,    42,    42,    47,   239,   240,    65,    46,    72,
      73,   171,   175,    42,   177,   173,    42,    78,    59,   120,
      53,   120,     9,    85,    10,    11,    12,   117,   118,   111,
     168,   197,   207,    83,    79,    54,    80,    78,   108,    13,
     141,   176,   141,    42,    42,     9,    32,    10,    11,    12,
     128,    39,    42,   130,     9,   132,    10,    11,    12,    76,
      77,    42,   214,    82,    58,   215,    42,    55,   210,   211,
     150,    33,   223,    39,    98,   146,    71,   103,   147,   148,
     101,   226,    80,    78,   227,   123,    99,   156,   159,   104,
     159,   134,   165,   166,   225,    65,   141,   141,   121,    74,
      75,    42,    42,    42,    76,    77,   172,   135,   174,    73,
      71,   208,   198,   133,   124,   125,   126,   127,   129,   200,
     117,   118,   141,   117,   118,    -7,   206,   209,     9,   143,
      10,    11,    12,   145,   228,   199,   179,   180,   181,   182,
     183,   184,   185,   117,   118,    13,   186,   187,   188,   189,
     190,   191,   117,   118,    74,    75,   151,    42,    42,    76,
      77,   169,   216,   192,    20,    21,    22,    23,   154,    24,
     155,   170,    25,   201,   202,    47,    48,   159,    34,    54,
      35,   203,    70,    42,   212,   213,   219,    37,    38,    49,
      50,   220,   229,   221,   222,    26,   159,   235,   236,    65,
     206,   178,   224,   230,   242,   179,   180,   181,   182,   183,
     184,   185,   120,   232,   237,   186,   187,   188,   189,   190,
     191,   194,   243,   245,   246,   179,   180,   181,   182,   183,
     184,   185,   192,   196,    56,   186,   187,   188,   189,   190,
     191,    20,    21,    22,    23,    34,    24,    35,   144,   136,
     244,   195,   192,    20,    21,    22,    23,    34,   107,    35,
      97,   136,    84,     0,    37,   137,   102,     1,   138,     2,
       0,     3,    26,     0,     0,     0,    37,   137,     0,    34,
     138,    35,     0,    70,    26,    20,    21,    22,    23,    34,
      24,    35,     0,    36,    20,    21,    22,    23,    34,    24,
      35,     0,    25,     9,     0,    10,    11,    12,    37,    38,
       0,     0,     0,     0,     0,     0,    26,    37,    38,     0,
     204,    20,    21,    60,    61,    26,    24,    62,     0,    63,
      20,    21,    22,    23,     0,    24,     0,     0,    25,     0,
       0,    20,    21,    60,    61,    64,   107,    62,     0,    63,
       0,     0,    26,    37,    38,    20,    21,    22,    23,     0,
      24,    26,     0,    25,     0,    64,     0,     0,    20,    21,
      60,    61,    26,    24,    62,     0,   131,   157,    21,    60,
      61,     0,    24,    62,     0,   131,    26,     0,    20,   160,
      60,    61,    64,    24,    62,     0,   131,    -6,     0,    26,
       9,    64,    10,    11,    12,     0,     0,     0,    26,     0,
       0,     0,    64,     0,     0,    20,    21,    60,    61,    26,
     107,    62,     0,   131,    20,    21,    22,    23,     0,    24,
      62,     0,    25,     0,     0,    20,    21,    22,    23,    64,
      24,   110,     0,    25,     0,     0,    26,     0,    64,     0,
       0,    20,    21,    22,    23,    26,    24,   149,     0,    25,
      20,    21,    22,    23,   152,    24,    26,     0,    25,     0,
       0,     0,     0,   112,   113,   114,   115,   116,   117,   118,
       0,     0,    26,     0,     0,   233,     0,     0,     0,    86,
       0,    26,   179,   180,   181,   182,   183,   184,   185,     0,
      49,    50,   186,   187,   188,   189,   190,   191,     0,    87,
      88,    89,    90,    91,    92,    93,    94,    95,     0,   192,
     218,   112,   113,   114,   115,   116,   117,   118
};

static const yytype_int16 yycheck[] =
{
       5,    13,    96,    63,    15,   115,    11,    12,   116,    23,
      42,    43,    17,    11,    12,    23,    12,    23,     8,    17,
       3,    87,     0,    37,    38,    30,    39,    42,    33,     4,
      96,    21,    37,    38,    45,    48,    49,    32,     6,    37,
      38,    23,   136,    48,   138,    25,    51,    29,    31,    57,
      48,    57,     9,    51,    11,    12,    13,    37,    38,    64,
     120,    23,   170,    46,    22,    19,    24,    29,    63,    26,
     136,   137,   138,    78,    79,     9,    22,    11,    12,    13,
      78,    79,    87,    79,     9,    80,    11,    12,    13,    42,
      43,    96,   202,    27,    99,   203,   101,    19,   192,   193,
     105,    22,    27,   101,    23,   101,    36,    22,   103,   104,
      22,   221,    24,    29,   222,    23,    29,   112,   113,    22,
     115,    27,   117,   118,   218,   120,   192,   193,    42,    37,
      38,   136,   137,   138,    42,    43,   131,    28,   133,   137,
      70,    25,    29,    24,    74,    75,    76,    77,    78,   154,
      37,    38,   218,    37,    38,     6,   167,    23,     9,    28,
      11,    12,    13,    28,   224,    29,    32,    33,    34,    35,
      36,    37,    38,    37,    38,    26,    42,    43,    44,    45,
      46,    47,    37,    38,    37,    38,    20,   192,   193,    42,
      43,    20,   204,    59,    14,    15,    16,    17,    42,    19,
      42,    48,    22,    20,    22,   216,    26,   202,    18,    19,
      20,    22,    22,   218,   198,   199,    23,    37,    38,    39,
      40,    23,    23,    29,    29,    45,   221,   232,   233,   224,
     241,    28,    28,    23,   239,    32,    33,    34,    35,    36,
      37,    38,    57,    28,    23,    42,    43,    44,    45,    46,
      47,    28,    20,    37,    20,    32,    33,    34,    35,    36,
      37,    38,    59,   145,    29,    42,    43,    44,    45,    46,
      47,    14,    15,    16,    17,    18,    19,    20,    99,    22,
     241,   144,    59,    14,    15,    16,    17,    18,    19,    20,
      53,    22,    48,    -1,    37,    38,    58,     3,    41,     5,
      -1,     7,    45,    -1,    -1,    -1,    37,    38,    -1,    18,
      41,    20,    -1,    22,    45,    14,    15,    16,    17,    18,
      19,    20,    -1,    22,    14,    15,    16,    17,    18,    19,
      20,    -1,    22,     9,    -1,    11,    12,    13,    37,    38,
      -1,    -1,    -1,    -1,    -1,    -1,    45,    37,    38,    -1,
      26,    14,    15,    16,    17,    45,    19,    20,    -1,    22,
      14,    15,    16,    17,    -1,    19,    -1,    -1,    22,    -1,
      -1,    14,    15,    16,    17,    38,    19,    20,    -1,    22,
      -1,    -1,    45,    37,    38,    14,    15,    16,    17,    -1,
      19,    45,    -1,    22,    -1,    38,    -1,    -1,    14,    15,
      16,    17,    45,    19,    20,    -1,    22,    14,    15,    16,
      17,    -1,    19,    20,    -1,    22,    45,    -1,    14,    15,
      16,    17,    38,    19,    20,    -1,    22,     6,    -1,    45,
       9,    38,    11,    12,    13,    -1,    -1,    -1,    45,    -1,
      -1,    -1,    38,    -1,    -1,    14,    15,    16,    17,    45,
      19,    20,    -1,    22,    14,    15,    16,    17,    -1,    19,
      20,    -1,    22,    -1,    -1,    14,    15,    16,    17,    38,
      19,    20,    -1,    22,    -1,    -1,    45,    -1,    38,    -1,
      -1,    14,    15,    16,    17,    45,    19,    20,    -1,    22,
      14,    15,    16,    17,    23,    19,    45,    -1,    22,    -1,
      -1,    -1,    -1,    32,    33,    34,    35,    36,    37,    38,
      -1,    -1,    45,    -1,    -1,    39,    -1,    -1,    -1,    28,
      -1,    45,    32,    33,    34,    35,    36,    37,    38,    -1,
      39,    40,    42,    43,    44,    45,    46,    47,    -1,    48,
      49,    50,    51,    52,    53,    54,    55,    56,    -1,    59,
      60,    32,    33,    34,    35,    36,    37,    38
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,     3,     5,     7,    63,     8,    21,    67,    68,     9,
      11,    12,    13,    26,    64,    65,    72,    78,    63,     0,
      14,    15,    16,    17,    19,    22,    45,    94,    95,    42,
      66,     4,    22,    22,    18,    20,    22,    37,    38,    90,
      91,    92,    94,    96,    92,    65,     6,    72,    26,    39,
      40,    86,    87,    90,    19,    19,    66,    70,    94,    63,
      16,    17,    20,    22,    38,    82,    83,    84,    94,    94,
      22,    91,    90,    90,    37,    38,    42,    43,    29,    22,
      24,    93,    27,    63,    87,    90,    28,    48,    49,    50,
      51,    52,    53,    54,    55,    56,    85,    86,    23,    29,
      69,    22,    93,    22,    22,    42,    43,    19,    82,    84,
      20,    94,    32,    33,    34,    35,    36,    37,    38,    23,
      57,    42,    73,    23,    91,    91,    91,    91,    90,    91,
      92,    22,    82,    24,    27,    28,    22,    38,    41,    89,
      90,    96,    89,    28,    70,    28,    92,    82,    82,    20,
      94,    20,    23,    23,    42,    42,    82,    14,    80,    82,
      15,    81,    82,    80,    81,    82,    82,    77,    84,    20,
      48,    23,    82,    25,    82,    89,    96,    89,    28,    32,
      33,    34,    35,    36,    37,    38,    42,    43,    44,    45,
      46,    47,    59,    88,    28,    69,    68,    23,    29,    29,
      94,    20,    22,    22,    26,    71,    72,    81,    25,    23,
      89,    89,    83,    83,    80,    81,    65,    74,    60,    23,
      23,    29,    29,    27,    28,    89,    80,    81,    84,    23,
      23,    75,    28,    39,    79,    94,    94,    23,    39,    48,
      49,    76,    94,    20,    71,    37,    20
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
#line 164 "parser.y"
    {}
    break;

  case 3:

/* Line 1455 of yacc.c  */
#line 166 "parser.y"
    {}
    break;

  case 4:

/* Line 1455 of yacc.c  */
#line 168 "parser.y"
    {}
    break;

  case 5:

/* Line 1455 of yacc.c  */
#line 170 "parser.y"
    {}
    break;

  case 6:

/* Line 1455 of yacc.c  */
#line 177 "parser.y"
    {
	/* The full program was parsed. Allocate and fill the final
	   .scop structures. */
	int nb_parameters, nb_arrays;

        parser_scop->parameters = clan_symbol_id_array(parser_clan_symbol,
                                                       SCOPLIB_TYPE_PARAMETER,
                                                       &nb_parameters);
        parser_scop->nb_parameters = nb_parameters;
        parser_scop->arrays = clan_symbol_id_array(parser_clan_symbol,
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
#line 220 "parser.y"
    {}
    break;

  case 12:

/* Line 1455 of yacc.c  */
#line 222 "parser.y"
    {}
    break;

  case 13:

/* Line 1455 of yacc.c  */
#line 228 "parser.y"
    {  parser_symbol_datatype = strdup((yyvsp[(1) - (1)].symbol)); }
    break;

  case 14:

/* Line 1455 of yacc.c  */
#line 230 "parser.y"
    { parser_symbol_datatype = strcat(strdup("struct "),strdup((yyvsp[(2) - (2)].symbol)));}
    break;

  case 15:

/* Line 1455 of yacc.c  */
#line 236 "parser.y"
    {}
    break;

  case 16:

/* Line 1455 of yacc.c  */
#line 238 "parser.y"
    {}
    break;

  case 17:

/* Line 1455 of yacc.c  */
#line 244 "parser.y"
    {}
    break;

  case 18:

/* Line 1455 of yacc.c  */
#line 246 "parser.y"
    {}
    break;

  case 19:

/* Line 1455 of yacc.c  */
#line 261 "parser.y"
    {     
    
    (yyval.symbol_table) = scoplib_symbol_add(&parser_scop_symbol,(yyvsp[(1) - (1)].symbol),SCOPLIB_TYPE_ARRAY,NULL,0); 
    parser_scop_symbol->data_type = parser_symbol_datatype;
    
    }
    break;

  case 20:

/* Line 1455 of yacc.c  */
#line 272 "parser.y"
    { 
       
      scoplib_matrix_p dimensions_bounds = (yyvsp[(2) - (2)].setex);     
      
      (yyval.symbol_table) = scoplib_symbol_add(&parser_scop_symbol,(yyvsp[(1) - (2)].symbol),SCOPLIB_TYPE_ARRAY,
                              dimensions_bounds,(yyvsp[(2) - (2)].setex)->NbRows); 
      parser_scop_symbol->data_type = parser_symbol_datatype;                              
                            
    }
    break;

  case 21:

/* Line 1455 of yacc.c  */
#line 286 "parser.y"
    { }
    break;

  case 24:

/* Line 1455 of yacc.c  */
#line 319 "parser.y"
    {
        clan_symbol_p symbol;
        symbol = clan_symbol_add(&parser_clan_symbol,(yyvsp[(3) - (3)].symbol),
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

  case 25:

/* Line 1455 of yacc.c  */
#line 341 "parser.y"
    {
        scoplib_vector_p parser_i_term = clan_vector_term(parser_clan_symbol,1,(yyvsp[(3) - (6)].symbol));
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

  case 26:

/* Line 1455 of yacc.c  */
#line 361 "parser.y"
    {
	scoplib_matrix_insert_matrix(parser_domain,(yyvsp[(9) - (9)].setex),parser_nb_cons);
        parser_nb_cons += (yyvsp[(9) - (9)].setex)->NbRows;
        parser_consperdim[parser_depth] += (yyvsp[(9) - (9)].setex)->NbRows;
      }
    break;

  case 27:

/* Line 1455 of yacc.c  */
#line 369 "parser.y"
    {
        parser_depth++;
        parser_scheduling[parser_depth] = 0;
      }
    break;

  case 28:

/* Line 1455 of yacc.c  */
#line 374 "parser.y"
    {
        parser_depth--;
        parser_scheduling[parser_depth]++;
        parser_nb_cons -= parser_consperdim[parser_depth];
        parser_consperdim[parser_depth] = 0;
	clan_symbol_remove(&parser_clan_symbol, parser_iterators[parser_depth]);
      }
    break;

  case 29:

/* Line 1455 of yacc.c  */
#line 386 "parser.y"
    {
	/* Insert the condition constraint in the current parser domain. */
	scoplib_matrix_insert_matrix(parser_domain,(yyvsp[(3) - (4)].setex),parser_nb_cons);
        parser_nb_cons += (yyvsp[(3) - (4)].setex)->NbRows;
      }
    break;

  case 30:

/* Line 1455 of yacc.c  */
#line 392 "parser.y"
    {
        parser_nb_cons -= (yyvsp[(3) - (6)].setex)->NbRows;
	/* Remove the condition constraint from the current parser domain. */
	int i, j;
	for (i = parser_nb_cons; i < parser_domain->NbRows - 1; ++i)
	  for (j = 0; j < parser_domain->NbColumns; ++j)
	    SCOPVAL_assign(parser_domain->p[i][j],parser_domain->p[i+1][j]);
      }
    break;

  case 31:

/* Line 1455 of yacc.c  */
#line 404 "parser.y"
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

  case 32:

/* Line 1455 of yacc.c  */
#line 415 "parser.y"
    {
	/* Deal with statements without surrounding loop by adding a
	   fake iterator */
	int old_parser_depth = parser_depth;
	if (parser_depth == 0)
	  {
	    char* fakeiter = strdup("fakeiter");
	    clan_symbol_p symbol = clan_symbol_lookup(parser_clan_symbol, fakeiter);
	    if (symbol)
	      free(fakeiter);
	    else
	      symbol = clan_symbol_add(&parser_clan_symbol,fakeiter,
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

  case 33:

/* Line 1455 of yacc.c  */
#line 473 "parser.y"
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

  case 34:

/* Line 1455 of yacc.c  */
#line 497 "parser.y"
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

  case 35:

/* Line 1455 of yacc.c  */
#line 527 "parser.y"
    {
        free((yyvsp[(1) - (2)].symbol));
      }
    break;

  case 36:

/* Line 1455 of yacc.c  */
#line 531 "parser.y"
    {
        free((yyvsp[(2) - (2)].symbol));
      }
    break;

  case 37:

/* Line 1455 of yacc.c  */
#line 535 "parser.y"
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

  case 38:

/* Line 1455 of yacc.c  */
#line 545 "parser.y"
    {
       if ((yyvsp[(3) - (3)].value) != 1)
	 {
	   yyerror("[Clan] Error: loop increment is not 1\n");
	   return 0;
	 }
       free ((yyvsp[(1) - (3)].symbol));
     }
    break;

  case 39:

/* Line 1455 of yacc.c  */
#line 563 "parser.y"
    {
	(yyval.setex) = scoplib_matrix_from_vector((yyvsp[(1) - (1)].affex));
        scoplib_vector_free((yyvsp[(1) - (1)].affex));
      }
    break;

  case 40:

/* Line 1455 of yacc.c  */
#line 569 "parser.y"
    {
       (yyval.setex) = scoplib_matrix_concat((yyvsp[(3) - (6)].setex), (yyvsp[(5) - (6)].setex));
     }
    break;

  case 41:

/* Line 1455 of yacc.c  */
#line 582 "parser.y"
    {
	(yyval.setex) = scoplib_matrix_from_vector((yyvsp[(1) - (1)].affex));
        scoplib_vector_free((yyvsp[(1) - (1)].affex));
      }
    break;

  case 42:

/* Line 1455 of yacc.c  */
#line 588 "parser.y"
    {
       (yyval.setex) = scoplib_matrix_concat((yyvsp[(3) - (6)].setex), (yyvsp[(5) - (6)].setex));
     }
    break;

  case 43:

/* Line 1455 of yacc.c  */
#line 601 "parser.y"
    {
        (yyval.affex) = (yyvsp[(1) - (1)].affex);
      }
    break;

  case 44:

/* Line 1455 of yacc.c  */
#line 605 "parser.y"
    {
        (yyval.affex) = scoplib_vector_add((yyvsp[(1) - (3)].affex),(yyvsp[(3) - (3)].affex));
        scoplib_vector_free((yyvsp[(1) - (3)].affex));
        scoplib_vector_free((yyvsp[(3) - (3)].affex));
      }
    break;

  case 45:

/* Line 1455 of yacc.c  */
#line 611 "parser.y"
    {
        (yyval.affex) = scoplib_vector_sub((yyvsp[(1) - (3)].affex),(yyvsp[(3) - (3)].affex));
	scoplib_vector_free((yyvsp[(1) - (3)].affex));
        scoplib_vector_free((yyvsp[(3) - (3)].affex));
      }
    break;

  case 46:

/* Line 1455 of yacc.c  */
#line 617 "parser.y"
    {
        (yyval.affex) = (yyvsp[(2) - (3)].affex);
      }
    break;

  case 47:

/* Line 1455 of yacc.c  */
#line 621 "parser.y"
    {
	SCOPVAL_assign((yyvsp[(3) - (6)].affex)->p[0], (yyvsp[(5) - (6)].affex)->p[(yyvsp[(5) - (6)].affex)->Size - 1]);
	(yyval.affex) = (yyvsp[(3) - (6)].affex);
      }
    break;

  case 48:

/* Line 1455 of yacc.c  */
#line 626 "parser.y"
    {
	SCOPVAL_assign((yyvsp[(3) - (6)].affex)->p[0], (yyvsp[(5) - (6)].affex)->p[(yyvsp[(5) - (6)].affex)->Size - 1]);
	(yyval.affex) = (yyvsp[(3) - (6)].affex);
      }
    break;

  case 49:

/* Line 1455 of yacc.c  */
#line 643 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_clan_symbol,(yyvsp[(1) - (1)].value),NULL);
      }
    break;

  case 50:

/* Line 1455 of yacc.c  */
#line 650 "parser.y"
    {
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(1) - (1)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_clan_symbol,1,(yyvsp[(1) - (1)].symbol));
        free((yyvsp[(1) - (1)].symbol));
      }
    break;

  case 51:

/* Line 1455 of yacc.c  */
#line 659 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_clan_symbol,-((yyvsp[(2) - (2)].value)),NULL);
      }
    break;

  case 52:

/* Line 1455 of yacc.c  */
#line 666 "parser.y"
    {
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(3) - (3)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_clan_symbol,(yyvsp[(1) - (3)].value),(yyvsp[(3) - (3)].symbol));
        free((yyvsp[(3) - (3)].symbol));
      }
    break;

  case 53:

/* Line 1455 of yacc.c  */
#line 675 "parser.y"
    {
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(1) - (3)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_clan_symbol,(yyvsp[(3) - (3)].value),(yyvsp[(1) - (3)].symbol));
        free((yyvsp[(1) - (3)].symbol));
      }
    break;

  case 54:

/* Line 1455 of yacc.c  */
#line 684 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_clan_symbol, ((yyvsp[(1) - (3)].value)) * ((yyvsp[(3) - (3)].value)), NULL);
      }
    break;

  case 55:

/* Line 1455 of yacc.c  */
#line 691 "parser.y"
    {
        (yyval.affex) = clan_vector_term(parser_clan_symbol, ((yyvsp[(1) - (3)].value)) / ((yyvsp[(3) - (3)].value)), NULL);
      }
    break;

  case 56:

/* Line 1455 of yacc.c  */
#line 698 "parser.y"
    {
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(4) - (4)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_clan_symbol,-((yyvsp[(2) - (4)].value)),(yyvsp[(4) - (4)].symbol));
        free((yyvsp[(4) - (4)].symbol));
      }
    break;

  case 57:

/* Line 1455 of yacc.c  */
#line 707 "parser.y"
    {
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(2) - (4)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_clan_symbol,-((yyvsp[(4) - (4)].value)),(yyvsp[(2) - (4)].symbol));
        free((yyvsp[(2) - (4)].symbol));
      }
    break;

  case 58:

/* Line 1455 of yacc.c  */
#line 716 "parser.y"
    {
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(2) - (2)].symbol),SCOPLIB_TYPE_UNKNOWN,parser_depth);
        (yyval.affex) = clan_vector_term(parser_clan_symbol,-1,(yyvsp[(2) - (2)].symbol));
        free((yyvsp[(2) - (2)].symbol));
      }
    break;

  case 59:

/* Line 1455 of yacc.c  */
#line 736 "parser.y"
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

  case 60:

/* Line 1455 of yacc.c  */
#line 771 "parser.y"
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

  case 61:

/* Line 1455 of yacc.c  */
#line 808 "parser.y"
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

  case 62:

/* Line 1455 of yacc.c  */
#line 841 "parser.y"
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

  case 63:

/* Line 1455 of yacc.c  */
#line 876 "parser.y"
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

  case 64:

/* Line 1455 of yacc.c  */
#line 893 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (3)].setex);
      }
    break;

  case 65:

/* Line 1455 of yacc.c  */
#line 900 "parser.y"
    {
       (yyval.setex) = scoplib_matrix_concat((yyvsp[(1) - (3)].setex),(yyvsp[(3) - (3)].setex));
       scoplib_matrix_free((yyvsp[(1) - (3)].setex));
       scoplib_matrix_free((yyvsp[(3) - (3)].setex));
     }
    break;

  case 76:

/* Line 1455 of yacc.c  */
#line 944 "parser.y"
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

  case 77:

/* Line 1455 of yacc.c  */
#line 957 "parser.y"
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

  case 78:

/* Line 1455 of yacc.c  */
#line 971 "parser.y"
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

  case 79:

/* Line 1455 of yacc.c  */
#line 984 "parser.y"
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

  case 80:

/* Line 1455 of yacc.c  */
#line 997 "parser.y"
    {
       (yyval.rw)[0] = (yyvsp[(1) - (2)].setex);
       (yyval.rw)[1] = NULL;
     }
    break;

  case 81:

/* Line 1455 of yacc.c  */
#line 1005 "parser.y"
    {
       (yyval.rw)[0] = (yyvsp[(2) - (3)].rw)[0];
       (yyval.rw)[1] = (yyvsp[(2) - (3)].rw)[1];
     }
    break;

  case 95:

/* Line 1455 of yacc.c  */
#line 1041 "parser.y"
    {
        (yyval.setex) = NULL;
      }
    break;

  case 96:

/* Line 1455 of yacc.c  */
#line 1048 "parser.y"
    {
        (yyval.setex) = NULL;
      }
    break;

  case 97:

/* Line 1455 of yacc.c  */
#line 1055 "parser.y"
    {
        (yyval.setex) = (yyvsp[(1) - (1)].setex);
      }
    break;

  case 98:

/* Line 1455 of yacc.c  */
#line 1063 "parser.y"
    {
        (yyval.setex) = scoplib_matrix_concat((yyvsp[(1) - (3)].setex),(yyvsp[(3) - (3)].setex));
	scoplib_matrix_free((yyvsp[(1) - (3)].setex));
        scoplib_matrix_free((yyvsp[(3) - (3)].setex));
      }
    break;

  case 99:

/* Line 1455 of yacc.c  */
#line 1072 "parser.y"
    {
        (yyval.setex) = (yyvsp[(2) - (2)].setex);
      }
    break;

  case 100:

/* Line 1455 of yacc.c  */
#line 1079 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (3)].setex);
      }
    break;

  case 101:

/* Line 1455 of yacc.c  */
#line 1086 "parser.y"
    {
	scoplib_matrix_p tmp = scoplib_matrix_concat((yyvsp[(1) - (5)].setex),(yyvsp[(3) - (5)].setex));
        (yyval.setex) = scoplib_matrix_concat(tmp,(yyvsp[(5) - (5)].setex));
	scoplib_matrix_free(tmp);
	scoplib_matrix_free((yyvsp[(1) - (5)].setex));
	scoplib_matrix_free((yyvsp[(3) - (5)].setex));
	scoplib_matrix_free((yyvsp[(5) - (5)].setex));
      }
    break;

  case 102:

/* Line 1455 of yacc.c  */
#line 1108 "parser.y"
    {
        int rank;
        scoplib_matrix_p matrix;
	char* s = (char*) (yyvsp[(1) - (1)].symbol);
	clan_symbol_p symbol = clan_symbol_lookup(parser_clan_symbol, s);
	// If the variable is an iterator or a parameter, discard it
	// from the read/write clause.
	if ((symbol && symbol->type == SCOPLIB_TYPE_ITERATOR) ||
	     (symbol && symbol->type == SCOPLIB_TYPE_PARAMETER))
	  (yyval.setex) = NULL;
	else
	  {
	    clan_symbol_add(&parser_clan_symbol, s, SCOPLIB_TYPE_ARRAY,parser_depth);
	    rank = clan_symbol_get_rank(parser_clan_symbol, s);
	    matrix = scoplib_matrix_malloc
	      (1, CLAN_MAX_DEPTH + CLAN_MAX_PARAMETERS + 2);
	    clan_matrix_tag_array(matrix, rank);
	    (yyval.setex) = matrix;
	  }
        free((yyvsp[(1) - (1)].symbol));
      }
    break;

  case 103:

/* Line 1455 of yacc.c  */
#line 1134 "parser.y"
    {
        int rank;
        clan_symbol_add(&parser_clan_symbol,(yyvsp[(1) - (2)].symbol),SCOPLIB_TYPE_ARRAY,parser_depth);
        rank = clan_symbol_get_rank(parser_clan_symbol,(yyvsp[(1) - (2)].symbol));
        clan_matrix_tag_array((yyvsp[(2) - (2)].setex),rank);
        (yyval.setex) = (yyvsp[(2) - (2)].setex);
        free((yyvsp[(1) - (2)].symbol));
      }
    break;

  case 104:

/* Line 1455 of yacc.c  */
#line 1147 "parser.y"
    {
	(yyval.setex) = (yyvsp[(3) - (4)].setex);
	free((yyvsp[(1) - (4)].symbol));
      }
    break;

  case 105:

/* Line 1455 of yacc.c  */
#line 1155 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (2)].setex);
      }
    break;

  case 106:

/* Line 1455 of yacc.c  */
#line 1162 "parser.y"
    {
	(yyval.setex) = (yyvsp[(2) - (2)].setex);
      }
    break;

  case 113:

/* Line 1455 of yacc.c  */
#line 1190 "parser.y"
    {
	(yyval.setex) = (yyvsp[(1) - (1)].setex);
      }
    break;

  case 114:

/* Line 1455 of yacc.c  */
#line 1197 "parser.y"
    {
	(yyval.setex) = scoplib_matrix_concat((yyvsp[(1) - (3)].setex),(yyvsp[(3) - (3)].setex));
      }
    break;

  case 115:

/* Line 1455 of yacc.c  */
#line 1204 "parser.y"
    {
	(yyval.setex) = (yyvsp[(1) - (3)].setex);
      }
    break;

  case 116:

/* Line 1455 of yacc.c  */
#line 1211 "parser.y"
    {
	(yyval.setex) = NULL;
      }
    break;

  case 117:

/* Line 1455 of yacc.c  */
#line 1218 "parser.y"
    {
	(yyval.setex) = NULL;
      }
    break;

  case 118:

/* Line 1455 of yacc.c  */
#line 1234 "parser.y"
    {
        (yyval.setex) = scoplib_matrix_from_vector((yyvsp[(2) - (3)].affex));
        scoplib_vector_free((yyvsp[(2) - (3)].affex));
      }
    break;

  case 119:

/* Line 1455 of yacc.c  */
#line 1242 "parser.y"
    {
	if ((yyvsp[(1) - (4)].setex) != NULL)
	  scoplib_matrix_insert_vector((yyvsp[(1) - (4)].setex),(yyvsp[(3) - (4)].affex),(yyvsp[(1) - (4)].setex)->NbRows);
        scoplib_vector_free((yyvsp[(3) - (4)].affex));
        (yyval.setex) = (yyvsp[(1) - (4)].setex);
      }
    break;

  case 120:

/* Line 1455 of yacc.c  */
#line 1263 "parser.y"
    {
       (yyval.symbol) = (yyvsp[(1) - (1)].symbol);
     }
    break;

  case 121:

/* Line 1455 of yacc.c  */
#line 1270 "parser.y"
    {
       (yyval.symbol) = (yyvsp[(2) - (3)].symbol);
     }
    break;

  case 122:

/* Line 1455 of yacc.c  */
#line 1277 "parser.y"
    {
       (yyval.symbol) = (yyvsp[(2) - (2)].symbol);
     }
    break;

  case 123:

/* Line 1455 of yacc.c  */
#line 1284 "parser.y"
    {
       (yyval.symbol) = NULL;
     }
    break;



/* Line 1455 of yacc.c  */
#line 2942 "parser.c"
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
#line 1296 "parser.y"




void
yyerror(char *s)
{
  fprintf(stderr, "%s\n", s);
  clan_parse_error = 1;
}

/**
  * clan_parser_get_symbols_from_clan_symbol_table function:
  * It tries to get the extra symbol information other than in pragma declaration
  * scop for smoothining of symbol table 
  */

void
clan_parser_get_symbols_from_clan_symbol_table(scoplib_symbol_p* parser_scop_symbol,
                                               clan_symbol_p parser_clan_symbol) {
                                               
  clan_symbol_p clan_symbol = parser_clan_symbol;
  while(clan_symbol != NULL ) {
    char* symbol_name = strdup(clan_symbol->identifier);
    if(scoplib_symbol_lookup(*parser_scop_symbol,symbol_name) == NULL ) {
      scoplib_symbol_add(parser_scop_symbol,symbol_name,clan_symbol->type,NULL,0);
    } 
    clan_symbol = clan_symbol->next;
  }
  
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
  parser_clan_symbol = NULL;
  parser_scop_symbol = NULL;

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
  clan_symbol_free(parser_clan_symbol);
//  scoplib_symbol_free(parser_scop_symbol);
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
  
  // printf("The parsing is successful.........\n");
  parser_scop->symbol_table = parser_scop_symbol;       
  
  if (! clan_parse_error)
    {
      if (parser_variables_localvars[0] != -1 ||
	  parser_variables_liveout[0] != -1)
	clan_scop_fill_options(parser_scop, parser_variables_localvars,
			       parser_variables_liveout);
      clan_scop_compact(parser_scop,CLAN_MAX_DEPTH);
    }
  else
    parser_scop = NULL;
  
 
  clan_parser_get_symbols_from_clan_symbol_table(&parser_scop_symbol,parser_clan_symbol);  
  parser_scop->symbol_table = parser_scop_symbol;      
  clan_parser_free_state();
  return parser_scop;
}

