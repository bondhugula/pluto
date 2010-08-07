/* types-polylib.h
     COPYRIGHT
          Both this software and its documentation are

              Copyright 1993, IRISA /Universite de Rennes I - France
              Copyright 1996,1997,1998, Doran Wilde and Vincent Loechner
              All rights reserved.

          Permission is granted to copy, use, and distribute
          for any commercial or noncommercial purpose under the terms
          of the GNU General Public license, version 2, June 1991
          (see file : LICENSING).
*/

#ifndef _types_polylib_h_
#define _types_polylib_h_

#ifdef GNUMP
#include<gmp.h>
#endif 

#include <limits.h>

/*********************** USER DEFINES ******************************/

/* first parameter name char.  */
#define FIRST_PARAMETER_NAME 'P'

/******************* END OF USER DEFINES ***************************/


#define PCHAR (FIRST_PARAMETER_NAME-1)
#define MAXNOOFRAYS 200 

#if defined(LINEAR_VALUE_IS_LONGLONG)
#define P_VALUE_FMT "%4lld "
#elif defined(LINEAR_VALUE_IS_LONG)
#define P_VALUE_FMT "%4ld "
#elif defined(LINEAR_VALUE_IS_CHARS)
#define P_VALUE_FMT "%s "
#elif defined(LINEAR_VALUE_IS_INT) 
#define P_VALUE_FMT "%4d "
#else  /* GNUMP */
#define P_VALUE_FMT "%4s "
#endif

/* Used in lower_upper_bounds */
#define LB_INFINITY 1
#define UB_INFINITY 2

/* MSB, TOP, and NEXT are defined over integer type, not on value type */
/* Put a one in the most significant bit of an int (portable) */
#define MSB ((unsigned)(((unsigned)1)<<(sizeof(int)*8-1)))

/* Largest representable positive number */
#define TOP ((int)(MSB-1))

/* Right shift the one bit in b and increment j if the last bit in b is one */
#define NEXT(j,b) { if (!((b)>>=1)) { (b)=MSB; (j)++; } }

/* Status of last Polyhedron operation */
extern int Pol_status;

#define POL_HIGH_BIT	(UINT_MAX - (UINT_MAX >> 1))
#define POL_NO_DUAL	(POL_HIGH_BIT | 0x0001)
#define POL_INTEGER	(POL_HIGH_BIT | 0x0002)
#define POL_ISSET(flags,f)  ((flags & f) == f)

typedef struct  {
  unsigned Size;
  Value *p;
} Vector;

typedef struct matrix {
  unsigned NbRows, NbColumns;
  Value **p;
  Value *p_Init;
  int p_Init_size;	/* needed to free the memory allocated by mpz_init */
} Matrix;

/* Macros to init/set/clear/test flags. */
#define FL_INIT(l, f)   (l) = (f)               /* Specific flags location. */
#define FL_SET(l, f)    ((l) |= (f))
#define FL_CLR(l, f)    ((l) &= ~(f))
#define FL_ISSET(l, f)  ((l) & (f))

#define F_INIT(p, f)    FL_INIT((p)->flags, f)  /* Structure element flags. */
#define F_SET(p, f)     FL_SET((p)->flags, f)
#define F_CLR(p, f)     FL_CLR((p)->flags, f)
#define F_ISSET(p, f)   FL_ISSET((p)->flags, f)

typedef struct polyhedron { 
  unsigned Dimension, NbConstraints, NbRays, NbEq, NbBid;
  Value **Constraint;
  Value **Ray;
  Value *p_Init;
  int p_Init_size;
  struct polyhedron *next;
#define    POL_INEQUALITIES	0x00000001
#define    POL_FACETS		0x00000002
#define    POL_POINTS		0x00000004
#define    POL_VERTICES		0x00000008
/* The flags field contains "valid" information,
 * i.e., the structure was created by PolyLib.
 */
#define	   POL_VALID		0x00000010
  unsigned flags;
} Polyhedron;

typedef struct interval {
  Value MaxN, MaxD;
  Value MinN, MinD; 
  int MaxI, MinI;
} Interval;

/* Test whether P is an empty polyhedron */
#define emptyQ(P) (P->NbRays==0)

/* Test whether P is a universe polyheron */
#define universeQ(P) (P->Dimension==P->NbBid)

typedef struct _Param_Vertex {  	
  Matrix *Vertex; /* Each row is a coordinate of the vertex. The first  */
	          /* "m" values of each row are the coefficients of the */
	          /* parameters. The (m+1)th value is the constant, the */
	          /* The (m+2)th value is the common denominator.       */
  Matrix *Domain; /* Constraints on parameters (in Polyhedral format)   */
  struct _Param_Vertex *next;          /* Pointer to the next structure */
} Param_Vertices;

typedef struct _Param_Domain {
  unsigned *F;         /* Bit array of faces */
  Polyhedron *Domain;  /* Pointer to Domain (constraints on parameters) */
  struct _Param_Domain *next; /* Pointer to the next structure  */
} Param_Domain;

typedef struct _Param_Polyhedron {
	int nbV;	    /* Number of parameterized vertices            */
	Param_Vertices *V;  /* Pointer to the list of parameteric vertices */
	Param_Domain *D;    /* Pointer to the list of validity domains     */
} Param_Polyhedron;

#define FORALL_PVertex_in_ParamPolyhedron(_V, _D, _P)   \
{     int _i, _ix;                                   \
      unsigned _bx;                                  \
      for( _i=0, _ix=0, _bx=MSB, _V=_P->V ;            \
           _V && (_i<_P->nbV) ; _i++, _V=_V->next )      \
      {       if (_D->F[_ix] & _bx)                   \
              {

#define END_FORALL_PVertex_in_ParamPolyhedron  \
              }                                \
              NEXT(_ix, _bx);                  \
      }                                        \
}

/* Data structures for pseudo-polynomial */

typedef enum { polynomial, periodic, evector } enode_type;

#ifdef CLN
#define POLY_UNION_OR_STRUCT struct
#else
#define POLY_UNION_OR_STRUCT union
#endif

typedef struct _evalue {
  Value d;              /* denominator */
  POLY_UNION_OR_STRUCT {
    Value n;            /* numerator (if denominator != 0) */
    struct _enode *p;	/* pointer   (if denominator == 0) */
  } x;
} evalue;

typedef struct _enode {
  enode_type type;      /* polynomial or periodic or evector */
  int size;             /* number of attached pointers */
  int pos;	        /* parameter position */
  evalue arr[1];        /* array of rational/pointer */
} enode;

typedef struct _enumeration {
  
  Polyhedron *ValidityDomain;    /* contraints on the parameters     */
  evalue EP;                     /* dimension = combined space       */
  struct _enumeration *next;     /* Ehrhart Polynomial, corresponding
	                            to parameter values inside the
                                    domain ValidityDomain below      */
} Enumeration;

/*-----------------------------Example Usage------------------------------*/
/* enode *e                                                               */
/*     e->type = polynomial     e->type = periodic   e->type = evector    */
/*     e->size = degree+1       e->size = period     e->size = length     */
/*     e->pos  = [1..nb_param]                                            */
/*     e->arr[i].d = denominator (Value)                                  */
/*     e->arr[i].x.p = pointer to another enode (if denominator is zero)  */
/*     e->arr[i].x.n = numerator (Value) (if denominator is non-zero)     */
/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/
/* This representation has the following advantages:                      */
/*   -- its dynamic, it can grow/shrink easily                            */
/*   -- it is easy to evaluate for a given context (values of parameters) */
/*   -- it allows pseudo-polynomial to be reduced with rules              */
/*   -- it can be constructed recursively                                 */
/*------------------------------------------------------------------------*/

/* *********************** |Represnting Z-Polyhedron| ******************* */


typedef enum {False = 0, True = 1} Bool;
typedef Matrix Lattice;
typedef struct LatticeUnion {
  Lattice *M;
  struct LatticeUnion *next;
} LatticeUnion;

typedef struct ZPolyhedron {
  Lattice *Lat ;
  Polyhedron *P;
  struct ZPolyhedron *next;
} ZPolyhedron;

#ifndef FOREVER
#define FOREVER for(;;)
#endif

#endif /* _types_polylib_h_ */












