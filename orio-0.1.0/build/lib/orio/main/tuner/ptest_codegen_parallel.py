#
# The code generator for performance-testing code
#

import re, sys
import ptest_codegen

#-----------------------------------------------------

CODE_TEMPLATE = r'''

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

typedef struct {
  int testid;
  char coord[1024];
  double tm;
} TimingInfo;

#ifndef REPS
#define REPS 1000
#endif

#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

/*@ decls @*/

void init_arrays()
{
/*@ inits @*/
}

int main(int argc, char *argv[])
{
  int numprocs, myid, _i;
  TimingInfo mytimeinfo;
  TimingInfo *timevec;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  /* Construct the MPI type for the timing info (what a pain!) */
  MPI_Datatype TimingInfoMPIType; 
  MPI_Datatype type[3] = {MPI_INT, MPI_CHAR, MPI_DOUBLE}; 
  int blocklen[3] = {1,1024,1};
  MPI_Aint disp[3];
  int base;
  MPI_Address( &mytimeinfo.testid, disp);
  MPI_Address( &mytimeinfo.coord, disp+1);
  MPI_Address( &mytimeinfo.tm, disp+2);
  base = disp[0];
  for (_i=0; _i <3; _i++) disp[_i] -= base;
  MPI_Type_struct( 3, blocklen, disp, type, &TimingInfoMPIType);
  MPI_Type_commit( &TimingInfoMPIType);
  /* end of MPI type construction */
  
  if (myid == 0) timevec = (TimingInfo*) malloc(numprocs * sizeof(TimingInfo));
  init_arrays();

  double annot_t_start=0, annot_t_end=0, annot_t_total=0;
  int annot_i;

  /*@ code @*/

  MPI_Gather(&mytimeinfo, 1, TimingInfoMPIType, timevec, 1, TimingInfoMPIType, 0, MPI_COMM_WORLD);
  
  if (myid == 0) {
    printf("{'%s' : %g", mytimeinfo.coord, mytimeinfo.tm);
    for (_i = 1; _i < numprocs; _i++) {
      printf(", '%s' : %g", timevec[_i].coord, timevec[_i].tm);
    }
    printf("}\n");
  }
      
  MPI_Finalize();
/*@ return @*/
}
                                                                   
'''

#-----------------------------------------------------

class PerfTestCodeGenParallel(ptest_codegen.PerfTestCodeGen):
    '''The code generator used to produce a performance-testing code'''

    # regular expressions
    __DECLS_TAG = r'/\*@\s*decls\s*@\*/'
    __INITS_TAG = r'/\*@\s*inits\s*@\*/'
    __CODE_TAG = r'/\*@\s*code\s*@\*/'
    __RETURN_TAG = r'/\*@\s*return\s*@\*/'

    #-----------------------------------------------------
    
    def __init__(self, input_params, input_decls):
        '''To instantiate the testing code generator'''
        ptest_codegen.PerfTestCodeGen.__init__(self, input_params, input_decls)

    #-----------------------------------------------------

    def generate(self, codedict):
        '''Generate the testing code used to get performance cost'''

        # generate the performance-testing code
        ptest_code = CODE_TEMPLATE
        ptest_code = re.sub(self.__DECLS_TAG, self.decl_code, ptest_code)
        ptest_code = re.sub(self.__INITS_TAG, self.init_code, ptest_code)
        ptest_code = re.sub(self.__RETURN_TAG, self.return_code, ptest_code)
        myid = 0
        code = ''
        for coord_key in codedict.keys():
            if myid == 0:
                code += 'if (myid == 0) {\n'
            else:
                code += '  } else if (myid == %s) {\n' % myid
            code += '    strcpy(mytimeinfo.coord,"%s");\n' % coord_key 
            code += '    mytimeinfo.testid = myid;\n'
            code += '    for (annot_i=0; annot_i<REPS; annot_i++) {\n' \
                    + '      annot_t_start = MPI_Wtime();\n\n'
                
            # The transformed code
            code += codedict[coord_key]
            code += '    }\n'
            code += '    annot_t_end = MPI_Wtime();\n' \
                + '    annot_t_total += annot_t_end - annot_t_start;\n' \
                + '    mytimeinfo.tm = annot_t_total / REPS;\n' \
                + ' \n\n' 
            myid += 1
        code += '}\n'
  
        ptest_code = re.sub(self.__CODE_TAG, code, ptest_code)

        # return the performance-testing code
        return ptest_code

    
