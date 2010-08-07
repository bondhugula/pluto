
   /*+------- <| --------------------------------------------------------**
    **         A                     Clan                                **
    **---     /.\   -----------------------------------------------------**
    **   <|  [""M#                 scanner.l                             **
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
 *::8:::::::::SUNDOGa8a::. Software Foundation, either version 3 of the       *
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


%{
   #include <stdlib.h>
   #include <string.h>
   #include <clan/macros.h>
   #include <clan/vector.h>
   #include <clan/matrix.h>
   #include <parser.h>

   void yyerror(char *);
   int  clan_scanner(int, char *, char *);

   extern int           parser_recording; /**< Do we record what we read ? */
   extern char *        parser_record;    /**< What we record */
   extern clan_symbol_p parser_symbol;    /**< Top of the symbol table */

   int  scanner_parsing = CLAN_FALSE;     /**< Do we parse or not ? */
   int  scanner_comment = 0;              /**< Comment nesting depth */
   char scanner_latest_text[SCOPLIB_MAX_STRING]; /**< Latest text read */
%}

%x LINECOMMENT
%x FULLCOMMENT

%%

<FULLCOMMENT>"*/"      {
                         scanner_comment --;
                         if (scanner_comment == 0)
                           BEGIN INITIAL;   /* Quit any mode */
                       }
<FULLCOMMENT>"/*"      { scanner_comment ++; }
<FULLCOMMENT>.         ; /* Do nothing */

<LINECOMMENT>\n        { BEGIN INITIAL;     /* Quit any mode */ }
<LINECOMMENT>.         ; /* Do nothing */

[ \t\n]*               ; /* Skip whitespaces */
"//"                   { BEGIN LINECOMMENT; /* Enter LINECOMMENT mode */ }
"/*"                   {
                         BEGIN FULLCOMMENT; /* Enter FULLCOMMENT mode */
                         scanner_comment ++;
                       }
"#pragma"[ ]+"scop"    { scanner_parsing = CLAN_TRUE;  /* Start parsing */ }
"#pragma"[ ]+"endscop" { scanner_parsing = CLAN_FALSE; /* Stop  parsing */ }
"#pragma"[ ]+"local-vars" { return clan_scanner(PRAGMALOCALVARS, yytext,"PRAGMALOCALVARS");   }
"#pragma"[ ]+"live-out" { return clan_scanner(PRAGMALIVEOUT, yytext,"PRAGMALIVEOUT");   }
"if"                   { return clan_scanner(IF,       yytext,"IF");   }
"else"                 { return clan_scanner(ELSE,     yytext,"ELSE"); }
"for"                  { return clan_scanner(FOR,      yytext,"FOR");  }
"min"                  { return clan_scanner(MIN,      yytext,"min");  }
"max"                  { return clan_scanner(MAX,      yytext,"max");  }
"ceild"                { return clan_scanner(CEILD,  yytext,"ceild");  }
"floord"               { return clan_scanner(FLOORD, yytext,"floord");  }

[_A-Za-z]"."[A-Za-z0-9_]* |
[_A-Za-z]"->"[A-Za-z0-9_]* |
[_A-Za-z][A-Za-z0-9_]* {
                         char * word;
                         word = (char *)malloc((strlen(yytext)+1)*sizeof(char));
			 strcpy(word,yytext);
			 yylval.symbol = word;
			 return clan_scanner(ID,yytext,"ID");
                       }
[0-9]+"."[0-9]*        { return clan_scanner(REAL,yytext,"REAL"); }
[0-9]+                 {
		         yylval.value = atoi(yytext) ;
                         return clan_scanner(INTEGER,yytext,"INTEGER");
		       }

"("                    { return clan_scanner(syRPARENTHESIS,   yytext,"(");  }
")"                    { return clan_scanner(syLPARENTHESIS,   yytext,")");  }
"["                    { return clan_scanner(syRBRACKET,       yytext,"[");  }
"]"                    { return clan_scanner(syLBRACKET,       yytext,"]");  }
"{"                    { return clan_scanner(syRBRACE,         yytext,"{");  }
"}"                    { return clan_scanner(syLBRACE,         yytext,"}");  }
","                    { return clan_scanner(syCOMMA,          yytext,",");  }
";"                    { return clan_scanner(sySEMICOLON,      yytext,";");  }
"."                    { return clan_scanner(syPOINT,          yytext,".");  }
"->"                   { return clan_scanner(syARROW,          yytext,"->"); }
":"                    { return clan_scanner(opCOLON,          yytext,":");  }
"!"                    { return clan_scanner(opNOT,            yytext,"!");  }
"?"                    { return clan_scanner(opQMARK,          yytext,"?");  }
"++"                   { return clan_scanner(opINCREMENTATION, yytext,"++"); }
"--"                   { return clan_scanner(opDECREMENTATION, yytext,"--"); }
"+="                   { return clan_scanner(opPLUSEQUAL,      yytext,"+="); }
"-="                   { return clan_scanner(opMINUSEQUAL,     yytext,"-="); }
"*="                   { return clan_scanner(opMULTIPLYEQUAL,  yytext,"*="); }
"/="                   { return clan_scanner(opDIVIDEEQUAL,    yytext,"/="); }
"|="                   { return clan_scanner(opOREQUAL,        yytext,"|="); }
"^="                   { return clan_scanner(opCOMPEQUAL,      yytext,"^="); }
"&="                   { return clan_scanner(opANDEQUAL,       yytext,"&="); }
"%="                   { return clan_scanner(opMODEQUAL,       yytext,"%="); }
"=="                   { return clan_scanner(opEQUAL,          yytext,"=="); }
"<="                   { return clan_scanner(opLEQ,            yytext,"<="); }
">="                   { return clan_scanner(opGEQ,            yytext,">="); }
"<"                    { return clan_scanner(opLOWER,          yytext,"<");  }
">"                    { return clan_scanner(opGREATER,        yytext,">");  }
"+"                    { return clan_scanner(opPLUS,           yytext,"+");  }
"-"                    { return clan_scanner(opMINUS,          yytext,"-");  }
"*"                    { return clan_scanner(opMULTIPLY,       yytext,"*");  }
"/"                    { return clan_scanner(opDIVIDE,         yytext,"/");  }
"%"                    { return clan_scanner(opMOD,            yytext,"%");  }
"&"                    { return clan_scanner(opAND,            yytext,"&");  }
"&&"                   { return clan_scanner(opLAND,           yytext,"&&"); }
"|"                    { return clan_scanner(opOR,             yytext,"|");  }
"||"                   { return clan_scanner(opLOR,            yytext,"||"); }
"^"                    { return clan_scanner(opCOMP,           yytext,"^");  }
"="                    { return clan_scanner(opASSIGNMENT,     yytext,"=");  }

.                      { if (scanner_parsing == CLAN_TRUE)
                         {
			    //if (CLAN_DEBUG)
			    //   printf("Lex: Unknown character (%s)\n",yytext);
			   //yyerror("Lex: Unknown character");
                         }
		       }

%%

int yywrap(void)
{
  return 1;
}


/**
 * clan_scanner function:
 * This function achieves some basic work when Lex recognize something: it
 * prints a log information if necessary for debugging, it updates the string
 * scanner_latest_text with the string that have been read and returns the
 * token if we are scanning a SCoP or the special token IGNORE otherwise.
 * \param token   The token code to send to Yacc for the Lex item.
 * \param text    The textual Lex item.
 * \param message A string to be printed for debugging.
 **
 * - 30/04/2008: first version.
 */

int clan_scanner(int token, char * text, char * message)
{
  //if (CLAN_DEBUG)
  //  printf("Lex: %s (%s)\n",message,text);

  sprintf(scanner_latest_text,"%s",text);
  if (parser_recording)
    strcat(parser_record,text);

  if (scanner_parsing == CLAN_TRUE)
    return token;
  else
    return IGNORE;
}
