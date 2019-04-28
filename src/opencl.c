/*
 * =====================================================================================
 *
 *       Filename:  opencl.c
 *
 *    Description:  contains functions to generate opencl code
 *
 *        Version:  1.0
 *        Created:  Saturday 28 July 2012 09:11:25  IST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Thejas C R
 *   Organization:
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <stdio.h>

/*
 * we have to generate the following functions.
 *
    opencl_codegen_generic_code(); // this function will generate the code that
 remains unchanged across programs. (like CreateContext etc)
    opencl_codegen_macros(); // generate the macros that conveys the size of our
 variables
    opencl_codegen_variables(); // generate the variables to be used by the CPU
 version of the opencl code
    opencl_codegen_headerfile(); // generates the header file and its entries
 specific to this program
    opencl_codegen_createbuffer(); // generate the createBuffers function
    opencl_codegen_runcomppackkernel(); // generate the runComputePackKernel
 code
    opencl_codegen_initdatadistrib(); // the InitialDataDistribution code
    opencl_codegen_setinitkernelargs();
    opencl_codegen_copyflowoutdata(); // the code for copying out data from
 their originating devices
    opencl_codegen_copyflowindata(); // the code for copying in the data to
 their receiving devices
    opencl_codegen_computeglobaloffsets(); // computing the kernel offsets to be
 used in each device
    opencl_codegen_enqueuecomputekernels(); // the GPU version of the
 runComputePackKernel
    opencl_codegen_cleanupresources(); // CleanupResources. Does the opposite of
 what CreateBuffers does
    opencl_codegen_destructor(); // the destructor
    opencl_codegen_constructor(); // the destructor
    opencl_codegen_printarray(); // the code to print the output array
    opencl_codegen_computekernel(); // the pluto generated tiled and paralleled
 code for CPU computations
*/

void opencl_codegen_generic_code();
void opencl_codegen_macros();
void opencl_codegen_variables();
void opencl_codegen_headerfile();
void opencl_codegen_generic_headerfile(char *filename);
void opencl_codegen_createbuffer();
void opencl_codegen_runcomppackkernel();
void opencl_codegen_initdatadistrib();
void opencl_codegen_setinitkernelargs();
void opencl_codegen_copyflowoutdata();
void opencl_codegen_copyflowindata();
void opencl_codegen_computeglobaloffsets();
void opencl_codegen_enqueuecomputekernels();
void opencl_codegen_cleanupresources();
void opencl_codegen_destructor();
void opencl_codegen_printarray();
void opencl_codegen_computekernel();

void pluto_opencl_codegen() {
  // we have to generate the following functions. for now we assume the opencl
  // kernel is given to us.

  opencl_codegen_headerfile("floyd");
  opencl_codegen_generic_code();
  opencl_codegen_macros();
  opencl_codegen_variables();
  opencl_codegen_createbuffer();
  opencl_codegen_runcomppackkernel();
  opencl_codegen_initdatadistrib();
  opencl_codegen_setinitkernelargs();
  opencl_codegen_copyflowoutdata();
  opencl_codegen_copyflowindata();
  opencl_codegen_computeglobaloffsets();
  opencl_codegen_enqueuecomputekernels();
  opencl_codegen_cleanupresources();
  opencl_codegen_destructor();
  opencl_codegen_printarray();
  opencl_codegen_computekernel();

  return;
}

void opencl_codegen_generic_code() {}

void opencl_codegen_macros() {}

void opencl_codegen_variables() {}

void opencl_codegen_generic_headerfile(char *filename) {
  char cmd[512] = { 0 };

  sprintf(cmd, "echo generic.hpp.templ >> %s.hpp", filename);

  system(cmd);

  return;
}

void opencl_codegen_headerfile(char *filename) {
  opencl_codegen_generic_headerfile(filename);

  return;
}

void opencl_codegen_createbuffer() {}

void opencl_codegen_runcomppackkernel() {}

void opencl_codegen_initdatadistrib() {}

void opencl_codegen_setinitkernelargs() {}

void opencl_codegen_copyflowoutdata() {}

void opencl_codegen_copyflowindata() {}

void opencl_codegen_computeglobaloffsets() {}

void opencl_codegen_enqueuecomputekernels() {}

void opencl_codegen_cleanupresources() {}

void opencl_codegen_destructor() {}

void opencl_codegen_printarray() {}

void opencl_codegen_computekernel() {}
