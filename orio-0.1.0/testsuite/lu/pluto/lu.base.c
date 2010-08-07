
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

double L[N][N];
double U[N][N];
double A[N][N+13];

void print_array() 
{ 
  int i, j; 
 
  for (i=0; i<N; i++) { 
    for (j=0; j<N; j++) { 
      fprintf(stderr, "%lf ", round(A[i][j])); 
      if (j%80 == 79) fprintf(stderr, "\n"); 
    } 
    fprintf(stderr, "\n"); 
  } 
  fprintf(stderr, "\n"); 
} 

void init_arrays()
{
    int i, j, k;

    /* have to initialize this matrix properly to prevent
    * division by zero
    */
    for (i=0; i<N; i++) {
     for (j=0; j<N; j++) {
      L[i][j] = 0.0;
      U[i][j] = 0.0;
     }
    }
    
    for (i=0; i<N; i++) {
     for (j=0; j<=i; j++) {
      L[i][j] = i+j+1;
      U[j][i] = i+j+1;
     }
    }
    
    for (i=0; i<N; i++) {
     for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
       A[i][j] += L[i][k]*U[k][j]; 
      }
     }
    }
}

double rtclock()
{
  struct timezone tzp;
  struct timeval tp;
  int stat;
  gettimeofday (&tp, &tzp);
  return (tp.tv_sec + tp.tv_usec*1.0e-6);
}

int main()
{
  init_arrays();

  double annot_t_start=0, annot_t_end=0, annot_t_total=0;
  int annot_i;

  for (annot_i=0; annot_i<REPS; annot_i++)
  {
    annot_t_start = rtclock();
   
    register int i,j,k;

    for (k=0; k<=N-1; k++) 
      {
	for (j=k+1; j<=N-1; j++)
	  A[k][j] = A[k][j]/A[k][k];
	for(i=k+1; i<=N-1; i++)   
	  for (j=k+1; j<=N-1; j++)   
	    A[i][j] = A[i][j] - A[i][k]*A[k][j];
      }
    
    annot_t_end = rtclock();
    annot_t_total += annot_t_end - annot_t_start;
  }

  //print_array();
  
  annot_t_total = annot_t_total / REPS;
  printf("%f\n", annot_t_total);
  
  return ((int) A[0][0]); 

}
                                    
