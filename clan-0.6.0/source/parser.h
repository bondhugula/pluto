/* A Bison parser, made by GNU Bison 2.4.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006,
   2009, 2010 Free Software Foundation, Inc.
   
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

/* Line 1685 of yacc.c  */
#line 104 "parser.y"
 int value;                     /**< An integer value for integers */
         char * symbol;                 /**< A string for identifiers */
         scoplib_vector_p affex;        /**< An affine expression */
         scoplib_matrix_p setex;        /**< A set of affine expressions */
         scoplib_matrix_p rw[2];        /**< Read and write array accesses */
         scoplib_symbol_p symbol_table;         
       


/* Line 1685 of yacc.c  */
#line 183 "parser.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


