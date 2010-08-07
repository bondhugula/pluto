#
# The skeleton code used for performance testing
#

import os, re, sys

#-----------------------------------------------------

DEFAULT = r'''

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

/*@ global @*/

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
  /*@ prologue @*/

  double orio_t_start=0, orio_t_end=0, orio_t_total=0;
  int orio_i;

  for (orio_i=0; orio_i<REPS; orio_i++)
  {
    orio_t_start = rtclock();
    
    /*@ tested code @*/

    orio_t_end = rtclock();
    orio_t_total += orio_t_end - orio_t_start;
  }
  
  orio_t_total = orio_t_total / REPS;
  printf("%f\n", orio_t_total);
  
  /*@ epilogue @*/

  return 0;
}
                                    
'''

#-----------------------------------------------------

class PerfTestSkeletonCode:
    '''The skeleton code used in the performance testing'''

    # tags
    __GLOBAL_TAG = r'/\*@\s*global\s*@\*/'
    __PROLOGUE_TAG = r'/\*@\s*prologue\s*@\*/'
    __TCODE_TAG = r'/\*@\s*tested\s+code\s*@\*/'
    __EPILOGUE_TAG = r'/\*@\s*epilogue\s*@\*/'

    #-----------------------------------------------------
    
    def __init__(self, code):
        '''To instantiate the skeleton code for the performance testing'''

        if code == None:
            code = DEFAULT

        self.__checkSkeletonCode(code)
        self.code = code

    #-----------------------------------------------------

    def __checkSkeletonCode(self, code):
        '''To check the validity of the skeleton code'''

        match_obj = re.search(self.__GLOBAL_TAG, code)
        if not match_obj:
            print 'error: missing "global" tag in the skeleton code'
            sys.exit(1)

        match_obj = re.search(self.__PROLOGUE_TAG, code)
        if not match_obj:
            print 'error: missing "prologue" tag in the skeleton code'
            sys.exit(1)

        match_obj = re.search(self.__TCODE_TAG, code)
        if not match_obj:
            print 'error: missing "tested code" tag in the skeleton code'
            sys.exit(1)

        match_obj = re.search(self.__EPILOGUE_TAG, code)
        if not match_obj:
            print 'error: missing "epilogue" tag in the skeleton code'
            sys.exit(1)

    #-----------------------------------------------------

    def insertCode(self, global_code, prologue_code, tested_code, epilogue_code):
        '''To insert code fragments into the skeleton code'''

        # initialize the performance-testing code
        code = self.code

        # insert code in the skeleton code
        code = re.sub(self.__GLOBAL_TAG, global_code, code)
        code = re.sub(self.__PROLOGUE_TAG, prologue_code, code)
        code = re.sub(self.__TCODE_TAG, tested_code, code)        
        code = re.sub(self.__EPILOGUE_TAG, epilogue_code, code)

        # return the performance-testing code
        return code
        
