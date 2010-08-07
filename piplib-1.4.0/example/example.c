/* This is a very simple example of how to use the PipLib inside your programs.
 * You should compile it by typing 'make' (after edition of the makefile), then
 * test it for instance by typing 'more FILE.pol | ./example'. Finally you can
 * compare results given by PIP by typing 'pip32 FILE.dat'
 */

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <piplib/piplib.h>

static PipOptions *options_read(FILE *f)
{
  char s[1024];
  PipOptions *options = pip_options_init();
  while (fgets(s, 1024, f)) {
    if (strncasecmp(s, "Maximize", 8) == 0)
      options->Maximize = 1;
    if (strncasecmp(s, "Urs_parms", 9) == 0)
      options->Urs_parms = 1;
    if (strncasecmp(s, "Urs_unknowns", 12) == 0)
      options->Urs_unknowns = 1;
    if (strncasecmp(s, "Rational", 8) == 0)
      options->Nq = 0;
    if (strncasecmp(s, "Dual", 4) == 0)
      options->Compute_dual = 1;
  }
  return options;
}

int main(int argc, const char **argv)
{ int bignum ;
  PipMatrix  * domain, * context  ;
  PipQuast   * solution ;
  PipOptions * options ;
  int verbose = 0;

  while (argc > 1) {
    if (strncmp(argv[1], "-v", 2) == 0) {
      const char *v = argv[1]+2;
      ++verbose;
      while (*v++ == 'v')
	++verbose;
    } else
      break;
    ++argv;
    --argc;
  }
  
  printf("[PIP2-like future input] Please enter:\n- the context matrix,\n") ;
  context = pip_matrix_read(stdin) ;
  pip_matrix_print(stdout,context) ;

  printf("- the bignum column (start at 0, -1 if no bignum),\n") ;
  fscanf(stdin," %d",&bignum) ;
  printf("%d\n",bignum) ;

  printf("- the constraint matrix.\n") ;
  domain = pip_matrix_read(stdin) ;
  pip_matrix_print(stdout,domain) ;
  printf("\n") ;
  
  if (isatty(0))
    printf("- options (EOF to stop).\n") ;
  options = options_read(stdin);
  options->Verbose = verbose;
  if (isatty(0))
    pip_options_print(stdout, options);

  /* The bignum in PIP1 is fixed on the constraint matrix, here is
   * the translation.
   */
  if (bignum > 0)
  bignum += domain->NbColumns - context->NbColumns ;
  
  solution = pip_solve(domain,context,bignum,options) ;

  pip_options_free(options) ;
  pip_matrix_free(domain) ;
  pip_matrix_free(context) ;

  pip_quast_print(stdout,solution,0) ;

  pip_quast_free(solution) ;
  return 0 ;
}
