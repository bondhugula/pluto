#
# The basic code generator for performance-testing code
#

import re, sys
import ptest_codegen, skeleton_code

#-----------------------------------------------------

class PerfTestCodeGenBasic(ptest_codegen.PerfTestCodeGen):
    '''The code generator used to produce a performance-testing code'''

    # function names
    malloc_func_name = 'malloc_arrays'
    init_func_name = 'init_input_vars'

    #-----------------------------------------------------

    def __init__(self, input_params, input_decls, decl_file, init_file, skeleton_code_file):
        '''To instantiate the testing code generator'''
        
        ptest_codegen.PerfTestCodeGen.__init__(self, input_params, input_decls)

        self.decl_file = decl_file
        self.init_file = init_file
        self.skeleton_code_file = skeleton_code_file
        
        self.__checkDeclFile()
        self.__checkInitFile()
        scode = self.__checkSkeletonCodeFile()

        self.ptest_skeleton_code = skeleton_code.PerfTestSkeletonCode(scode)
        
    #-----------------------------------------------------

    def __checkDeclFile(self):
        '''To check the declaration file'''

        # do not perform checking if no declaration file is specified
        if not self.decl_file:
            return

        # check if the file can be opened
        try:
            f = open(self.decl_file)
            f.close()
        except:
            print 'error: cannot read file: "%s"' % self.decl_file
            sys.exit(1)

    #-----------------------------------------------------

    def __checkInitFile(self):
        '''To check the initialization file'''

        # do not perform checking if no initialization file is specified
        if not self.init_file:
            return

        # read the content of the file
        try:
            f = open(self.init_file)
            init_code = f.read()
            f.close()
        except:
            print 'error: cannot read file: "%s"' % self.init_file
            sys.exit(1)

        # check if the file contains the initialization function
        init_func_re = r'void\s+%s\(\s*\)\s*\{' % self.init_func_name
        match_obj = re.search(init_func_re, init_code)
        if not match_obj:
            print (('error: no initialization function (named "%s") can be found in the ' +
                    'initialization file: "%s"') % (self.init_func_name, self.init_file))
            sys.exit(1)

    #-----------------------------------------------------

    def __checkSkeletonCodeFile(self):
        '''To check the skeleton-code file, and return the skeleton code'''

        # do not perform checking if no skeleton-code file is specified
        if not self.skeleton_code_file:
            return None

        # read the content of the file
        try:
            f = open(self.skeleton_code_file)
            skton_code = f.read()
            f.close()
        except:
            print 'error: cannot read file: "%s"' % self.decl_file
            sys.exit(1)

        # return the skeleton code
        return skton_code

    #-----------------------------------------------------

    def generate(self, code):
        '''Generate the testing code used to get performance cost'''

        # generate the macro definition codes for the input parameters
        iparam_code = self.iparam_code

        # generate the declaration code
        if self.decl_file:
            decl_code = '#include "%s"\n' % self.decl_file
        else:
            decl_code = self.decl_code + '\n'
            decl_code += 'void %s() {\n%s\n}\n' % (self.malloc_func_name, self.malloc_code)

        # generate the initialization code
        if self.init_file:
            init_code = '#include "%s"\n' % self.init_file
        else:
            init_code = 'void %s() {\n%s\n}\n' % (self.init_func_name, self.init_code)

        # create code for the global definition section
        global_code = ''
        global_code += iparam_code + '\n'
        global_code += decl_code + '\n'
        global_code += init_code + '\n'

        # create code for the prologue section
        prologue_code = ''
        if not self.decl_file:
            prologue_code += ('%s();' % self.malloc_func_name) + '\n'
        prologue_code += ('%s();' % self.init_func_name) + '\n'

        # create code for the epilogue section
        epilogue_code = ''

        # get the performance-testing code
        ptest_code = self.ptest_skeleton_code.insertCode(global_code, prologue_code, 
                                                         code, epilogue_code)

        # return the performance-testing code
        return ptest_code

    
