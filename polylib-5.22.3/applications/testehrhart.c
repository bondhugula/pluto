/***********************************************************************/
/*                Ehrhart V4.20                                        */
/*                copyright 1997, Doran Wilde                          */
/*                copyright 1997-2000, Vincent Loechner                */
/*       Permission is granted to copy, use, and distribute            */
/*       for any commercial or noncommercial purpose under the terms   */
/*       of the GNU General Public license, version 2, June 1991       */
/*       (see file : LICENSING).                                       */
/***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include <polylib/polylib.h>
#include <polylib/homogenization.h>
#include "config.h"

#ifndef HAVE_GETOPT_H
#define getopt_long(a,b,c,d,e) getopt(a,b,c)
#else
#include <getopt.h>
struct option options[] = {
    { "homogenized",   no_argument,  0,  'h' },
    { 0, 0, 0, 0 }
};
#endif

#define WS 0

/** 
    
define this to print all constraints on the validity domains if not
defined, only new constraints (not in validity domain given by the
user) are printed

*/
#define EPRINT_ALL_VALIDITY_CONSTRAINTS

/** 

The following are mainly for debug purposes. You shouldn't need to
change anything for daily usage...

*/

/** you may define each macro independently 
<ol>
<li> #define EDEBUG minimal debug 
<li> #define EDEBUG1 prints enumeration points
<li> #define EDEBUG11 prints number of points
<li> #define EDEBUG2 prints domains
<li> #define EDEBUG21 prints more domains
<li> #define EDEBUG3 prints systems of equations that are solved
<li> #define EDEBUG4 prints message for degree reduction
<li> #define EDEBUG5 prints result before simplification 
<li> #define EDEBUG6 prints domains in Preprocess 
<li> #define EDEBUG61 prints even more in Preprocess
<li> #define EDEBUG62 prints domains in Preprocess2
</ol>
*/

/* #define EDEBUG 	*/		/* minimal debug */
/* #define EDEBUG1	*/		/* prints enumeration points */
/* #define EDEBUG11	*/		/* prints number of points */
/* #define EDEBUG2	*/		/* prints domains */
/* #define EDEBUG21	*/		/* prints more domains */
/* #define EDEBUG3	*/		/* prints systems of equations that are solved */
/* #define EDEBUG4	*/		/* prints message for degree reduction */
/* #define EDEBUG5	*/		/* prints result before simplification */
/* #define EDEBUG6	*/		/* prints domains in Preprocess */
/* #define EDEBUG61	*/		/* prints even more in Preprocess */
/* #define EDEBUG62	*/		/* prints domains in Preprocess2 */


/**

 Reduce the degree of resulting polynomials

*/
#define REDUCE_DEGREE

/** 

define this to print one warning message per domain overflow these
overflows should no longer happen since version 4.20

*/
#define ALL_OVERFLOW_WARNINGS

/**

EPRINT : print results while computing the ehrhart polynomial.  this
is done by default if you build the executable ehrhart.  (If EMAIN is
defined).  Don't define EMAIN here, it is defined when necessary in
the makefile.  

<p>

Notice: you may however define EPRINT without defining EMAIN, but in
this case, you have to initialize the global variable param_name by
calling Read_ParamNames before any call to ehrhart.  This is NOT
recommanded, unless you know what you do.  EPRINT causes more debug
messages to be printed.

*/
/* #define EPRINT */

int main(int argc, char **argv)
{
    int i;
    char str[1024];
    Matrix *C1, *P1;
    Polyhedron *C, *P;
    Enumeration *en;
    char **param_name;
    int c, ind = 0;
    int hom = 0;
  
#ifdef EP_EVALUATION
    Value *p, *tmp;
    int k;
#endif

    while ((c = getopt_long(argc, argv, "h", options, &ind)) != -1) {
	switch (c) {
	case 'h':
	    hom = 1;
	    break;
	}
    }

    P1 = Matrix_Read();
    C1 = Matrix_Read();
    if(C1->NbColumns < 2) {
        fprintf( stderr, "Not enough parameters !\n" );
        exit(0);
    }
    if (hom) {
	Matrix *C2, *P2;
	P2 = AddANullColumn(P1);
	Matrix_Free(P1);
	P1 = P2;
	C2 = AddANullColumn(C1);
	Matrix_Free(C1);
	C1 = C2;
    }
    P = Constraints2Polyhedron(P1,WS);
    C = Constraints2Polyhedron(C1,WS);
    Matrix_Free(P1);
    Matrix_Free(C1);
  
    /* Read the name of the parameters */
    param_name = Read_ParamNames(stdin,C->Dimension - hom);
    if (hom) {
	char **param_name2 = (char**)malloc(sizeof(char*) * (C->Dimension));
	for (i = 0; i < C->Dimension - 1; i++)
	    param_name2[i] = param_name[i];
	param_name2[C->Dimension-1] = "_H";
	free(param_name);
	param_name=param_name2;
    }

    en = Polyhedron_Enumerate(P,C,WS,param_name);

    if (hom) {
	Enumeration *en2;

	printf("inhomogeneous form:\n");
      
	dehomogenize_enumeration(en, C->Dimension, WS);
	for (en2 = en; en2; en2 = en2->next) {
	    Print_Domain(stdout, en2->ValidityDomain, param_name);
	    print_evalue(stdout, &en2->EP, param_name);
	}
    }

#ifdef EP_EVALUATION
    if( isatty(0) && C->Dimension != 0)
        {  /* no tty input or no polyhedron -> no evaluation. */
            printf("Evaluation of the Ehrhart polynomial :\n");
            p = (Value *)malloc(sizeof(Value) * (C->Dimension));
            for(i=0;i<C->Dimension;i++) 
                value_init(p[i]);
            FOREVER {
                fflush(stdin);
                printf("Enter %d parameters : ",C->Dimension);
                for(k=0;k<C->Dimension;++k) {
                    scanf("%s",str);
                    value_read(p[k],str);
                }
                fprintf(stdout,"EP( ");
                value_print(stdout,VALUE_FMT,p[0]);
                for(k=1;k<C->Dimension;++k) {
                    fprintf(stdout,",");
                    value_print(stdout,VALUE_FMT,p[k]);
                }  
                fprintf(stdout," ) = ");
                value_print(stdout,VALUE_FMT,*(tmp=compute_poly(en,p)));
                free(tmp);
                fprintf(stdout,"\n");  
            }
        }
#endif /* EP_EVALUATION */
  
    Enumeration_Free(en);
    for( i=0 ; i < C->Dimension - hom; i++ )
        free( param_name[i] );
    free(param_name);
    Polyhedron_Free( P );
    Polyhedron_Free( C );

    return 0;
}

