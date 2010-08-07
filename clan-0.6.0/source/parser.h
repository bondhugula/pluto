
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
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

/* Line 1676 of yacc.c  */
#line 90 "parser.y"
 int value;                     /**< An integer value for integers */
         char * symbol;                 /**< A string for identifiers */
         scoplib_vector_p affex;        /**< An affine expression */
         scoplib_matrix_p setex;        /**< A set of affine expressions */
         scoplib_matrix_p rw[2];        /**< Read and write array accesses */
       


/* Line 1676 of yacc.c  */
#line 171 "parser.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;


