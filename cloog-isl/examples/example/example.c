/* This is a very simple example of how to use the CLooGLib inside your
 * programs. You should compile it by typing 'make' (after edition of the
 * makefile), then test it for instance by typing
 * 'more FILE.cloog | ./example' (or example.exe under Cygwin).
 */

# include <stdio.h>
# include <cloog/cloog.h>

int main()
{
  CloogState *state;
  CloogProgram *program;
  CloogOptions * options ;
  
  state = cloog_state_malloc();
  options = cloog_options_malloc(state);
  program = cloog_program_read(stdin,options) ;
  
  program = cloog_program_generate(program,options) ;
  cloog_program_pprint(stdout,program,options) ;

  cloog_options_free(options) ;
  cloog_program_free(program) ;
  cloog_state_free(state);
  
  return 0 ;
}
