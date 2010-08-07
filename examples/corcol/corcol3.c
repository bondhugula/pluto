#include <stdio.h>
#include <unistd.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#define my_sqrt_array(x,j) sqrt(x[j])

#define N 2000
#define M 2000

#define FLOAT_TYPE float

FLOAT_TYPE float_n = (FLOAT_TYPE) N;
FLOAT_TYPE eps = 0.005;
FLOAT_TYPE stddev[M + 1];
FLOAT_TYPE data[M + 1][N + 1];
FLOAT_TYPE mean[M + 1];
FLOAT_TYPE symmat[M + 1][M + 1];


#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif

void init_array()
{
    int i, j;

    for (i=0; i<M; i++) {
        for (j=0; j<N; j++) {
	  data[i][j] = ((FLOAT_TYPE) i*j)/N;
        }
    }
}


void print_array(char** argv)
{
    int i, j;
    if (! strcmp(argv[0], ""))
      {
	for (i=0; i<M; i++) {
	  for (i=0; i<M; i++) {
	  fprintf(stderr, "%0.2lf ", symmat[i][j]);
	  if (i%80 == 20) fprintf(stderr, "\n");
	  }
	}
	fprintf(stderr, "\n");
      }
}


double rtclock()
{
    struct timezone Tzp;
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, &Tzp);
    if (stat != 0) printf("Error return from gettimeofday: %d",stat);
    return(Tp.tv_sec + Tp.tv_usec*1.0e-6);
}

#ifdef TIME
#define IF_TIME(foo) foo;
#else
#define IF_TIME(foo)
#endif


int main(int argc, char** argv)
{
  double t_start, t_end;
  int i, j, j1, j2;
  int n = N;
  int m = M;
  init_array();

  IF_TIME(t_start = rtclock());



  /* Determine mean of column vectors of input data matrix */
#pragma scop
  for (j = 1; j <= m; j++)
    {
      mean[j] = 0.0;
      for (i = 1; i <= n; i++)
	mean[j] += data[i][j];
      mean[j] /= float_n;
    }

/* Determine standard deviations of column vectors of data matrix. */
  for (j = 1; j <= m; j++)
    {
      stddev[j] = 0.0;
      for (i = 1; i <= n; i++)
	stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] /= float_n;
      stddev[j] = my_sqrt_array(stddev, j);
      /* The following in an inelegant but usual way to handle
	 near-zero std. dev. values, which below would cause a zero-
	 divide. */
      stddev[j] = stddev[j] <= eps ? 1.0 : stddev[j];
    }
#pragma endscop

 /* Center and reduce the column vectors. */
  for (i = 1; i <= n; i++)
    for (j = 1; j <= m; j++)
      {
	data[i][j] -= mean[j];
	data[i][j] /= sqrt(float_n) * stddev[j];
      }

/* Calculate the m * m correlation matrix. */
  for (j1 = 1; j1 <= m-1; j1++)
    {
      symmat[j1][j1] = 1.0;
      for (j2 = j1+1; j2 <= m; j2++)
	{
	  symmat[j1][j2] = 0.0;
	  for (i = 1; i <= n; i++)
	    symmat[j1][j2] += ( data[i][j1] * data[i][j2]);
	  symmat[j2][j1] = symmat[j1][j2];
	}
    }
 symmat[m][m] = 1.0;

 IF_TIME(t_end = rtclock());
 IF_TIME(fprintf(stderr, "%0.6lfs\n", t_end - t_start));


 print_array(argv);


 return 0;
}
