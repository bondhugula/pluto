/*************************************************/
/*     verif_ehrhart.c                           */
/* program to compare effective number of points */
/* in a polytope with the corresponding          */
/* evaluation of the Ehrhart polynomial.         */
/* Parameters vary in range -RANGE to RANGE      */
/* (define below) by default.                    */
/* Can be overridden by specifying               */
/* -r<RANGE>, or -m<min> and -M<max>             */
/*                                               */
/* written by Vincent Loechner (c) 2000.         */
/*  loechner@icps.u-strasbg.fr                   */
/*************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <polylib/polylib.h>
#define MAXRAYS 1024

/* define this to print all the results */
/* else, only a progress bar is printed */
/* #define PRINT_ALL_RESULTS	 */
/* define this to continue the test after first error found */
/* #define DONT_BREAK_ON_ERROR */

/* RANGE : normal range for evalutations (-RANGE -> RANGE) */
#define RANGE 50

/* SRANGE : small range for evalutations */
#define SRANGE 15

/* if dimension >= BIDDIM, use SRANGE */
#define BIGDIM 5

/* VSRANGE : very small range for evalutations */
#define VSRANGE 5

/* if dimension >= VBIDDIM, use VSRANGE */
#define VBIGDIM 8

Value min, max;

#ifdef DONT_BREAK_ON_ERROR
#define PRINT_ALL_RESULTS
#endif

#ifndef PRINT_ALL_RESULTS
int st;
#endif

/****************************************************/
/* function check_poly :                            */
/* scans the parameter space from min to max (all   */
/* directions). Computes the number of points in    */
/* the polytope using both methods, and compare them*/
/* returns 1 on success                             */
/****************************************************/

int check_poly(Polyhedron *S,Polyhedron *C,Enumeration *en,
	       int nparam,int pos,Value *z) {
  
  int cc,k;
  Value c,tmp,*ctmp;
  
  value_init(c); value_init(tmp);
  
  if(pos == nparam) {
    
    /* Computes the ehrhart polynomial */
    value_assign(c,*(ctmp=compute_poly(en,&z[S->Dimension-nparam+1])));
    free(ctmp);
    /* if c=0 we may be out of context. */
    /* scanning is useless in this case*/
    if(!in_domain(C,&z[S->Dimension-nparam+1])) {
   
      /* ok */ ;
    }
    else {
      
#ifdef PRINT_ALL_RESULTS
      printf("EP( ");
      value_print(stdout,VALUE_FMT,z[S->Dimension-nparam+1]);
      for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
	printf(", ");
	value_print(stdout,VALUE_FMT,z[k]);
      }
      printf(" ) = ");
      value_print(stdout,VALUE_FMT,c);
      printf(" ");
#endif

      /* Count manually the number of points */
      count_points(1,S,z,&tmp);
#ifdef PRINT_ALL_RESULTS
	printf(", count = ");
	value_print(stdout, P_VALUE_FMT, tmp);
	printf(". ");
#endif

      if(value_ne(tmp,c)) {
	printf("\n"); 
	fflush(stdout);
	fprintf(stderr,"Error !\n");
	fprintf(stderr,"EP( ");
	value_print(stderr,VALUE_FMT,z[S->Dimension-nparam+1]);
	for(k=S->Dimension-nparam+2;k<=S->Dimension;++k) {
	  fprintf(stderr,", ");
	  value_print(stderr,VALUE_FMT,z[k]);
	}
	fprintf(stderr," ) should be ");
	value_print(stderr,VALUE_FMT,tmp);
	fprintf(stderr,", while EP eval gives ");
	value_print(stderr,VALUE_FMT,c);
	fprintf(stderr,".\n");
#ifndef DONT_BREAK_ON_ERROR
	value_clear(c); value_clear(tmp);
	return(0);
#endif
      }

#ifdef PRINT_ALL_RESULTS
      else
	printf("OK.\n");
#endif
    }
  }
  else
    for(value_assign(tmp,min); value_le(tmp,max); value_increment(tmp,tmp)) {

#ifndef PRINT_ALL_RESULTS
      k = VALUE_TO_INT(tmp);
      if(!pos && !(k%st)) {
	printf("o");
	fflush(stdout);
      }
#endif
      
      value_assign(z[pos+S->Dimension-nparam+1],tmp);
      if(!check_poly(S,C,en,nparam,pos+1,z)) {
	value_clear(c); value_clear(tmp);
	return(0);
      }
    }
  value_clear(c); value_clear(tmp);
  return(1);
} /* check_poly */

int main(int argc,char *argv[]) {
	
  Matrix *C1, *P1;
  Polyhedron *C, *P, *S;
  Polyhedron *CC, *PP;
  Enumeration *en;
  Value *p, tmp;
  int i,j;
  int m,M;
  
/******* Read the input *********/
  P1 = Matrix_Read();
  C1 = Matrix_Read();

  if(C1->NbColumns < 2) {
    fprintf(stderr,"Not enough parameters !\n");
    exit(0);
  }
  
  P = Constraints2Polyhedron(P1,MAXRAYS);
  C = Constraints2Polyhedron(C1,MAXRAYS);
  Matrix_Free(C1);
  Matrix_Free(P1);

  /******* Read the options: initialize min and max ********/
  if(P->Dimension >= VBIGDIM)
    M = VSRANGE;
  else if(P->Dimension >= BIGDIM)
    M = SRANGE;
  else
    M = RANGE;
  m = -M;
  if(argc != 1 ) {
    for(i=1;i<argc;i++) {
      if(!strncmp(argv[i],"-m",2)) {
	
	/* min specified */
	m = atoi(&argv[i][2]);
      }
      else if(!strncmp(argv[i],"-M",2)) {
	
	/* max specified */
	M = atoi(&argv[i][2]);
      }
      else if(!strncmp(argv[i], "-r", 2)) {
	
	/* range specified */
	M = atoi(&argv[i][2]);
	m = -M;
      }
      else {
	fprintf(stderr,"Unknown option: %s\n",argv[i]);
	fprintf(stderr,"Usage: %s [-m<>][-M<>][-r<>]\n",argv[0]);
	return(-1);
      }
    }
  }
  if(m > M) {
    fprintf(stderr,"Nothing to do: min > max !\n");
    return(0);
  }
  value_init(min);
  value_init(max);
  value_set_si(min,m);
  value_set_si(max,M);
  value_init(tmp);

  /******* Compute true context *******/
  CC = align_context(C,P->Dimension,MAXRAYS);
  PP = DomainIntersection(P,CC,MAXRAYS);
  Domain_Free(CC);
  C1 = Matrix_Alloc(C->Dimension+1,P->Dimension+1);

  for(i=0;i<C1->NbRows;i++)
    for(j=0;j<C1->NbColumns;j++)
      if(i==j-P->Dimension+C->Dimension)
	value_set_si(C1->p[i][j],1);
      else
	value_set_si(C1->p[i][j],0);
  CC = Polyhedron_Image(PP,C1,MAXRAYS);
  Domain_Free(C);
  C = CC;

  /******* Compute EP *********/
  en = Polyhedron_Enumerate(P,C,MAXRAYS,NULL);
  
  /******* Initializations for check *********/
  p = (Value *)malloc(sizeof(Value) * (P->Dimension+2));
  for(i=0;i<=P->Dimension;i++) {
    value_init(p[i]);
    value_set_si(p[i],0);
  }
  value_init(p[i]);
  value_set_si(p[i],1);

  /* S = scanning list of polyhedra */
  S = Polyhedron_Scan(P,C,MAXRAYS);

#ifndef PRINT_ALL_RESULTS
  if(C->Dimension > 0) {
    value_subtract(tmp,max,min);
    if (VALUE_TO_INT(tmp) > 80)
      st = 1+(VALUE_TO_INT(tmp))/80;
    else
      st=1;
    for(i=VALUE_TO_INT(min);i<=VALUE_TO_INT(max);i+=st)
      printf(".");
    printf( "\r" );
    fflush(stdout);
  }
#endif

  /******* CHECK NOW *********/
  if(S && !check_poly(S,C,en,C->Dimension,0,p)) {
    fprintf(stderr,"Check failed !\n");
    for(i=0;i<=(P->Dimension+1);i++) 
      value_clear(p[i]);
    value_clear(tmp);  
    return(-1);
  }
    
#ifndef PRINT_ALL_RESULTS
  printf( "\n" );
#endif
  
  for(i=0;i<=(P->Dimension+1);i++) 
    value_clear(p[i]);
  value_clear(tmp);
  return(0);
} /* main */




