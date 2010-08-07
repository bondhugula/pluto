#
# Contain unparsing procedures.
#

import sys

#-------------------------------------------

class CodeGen:
    '''The code generator for the Blue Gene's memory alignment optimizer'''

    def __init__(self, vars, annot_body_code, indent):
        '''To instantiate a code generator instance'''
        
        self.vars = vars
        self.annot_body_code = annot_body_code
        self.indent = indent

    #------------------------------------------------------

    def __printAddress(self, var):
        '''To Return the starting address location of the given variable (in C/C++)'''

        vname = var.vname
        dims = var.dims[:]
        dims.remove(None)
        s = str(vname)
        if len(dims) > 0:
            s += '[' + ']['.join(map(str, dims)) + ']'
        return s
        
    #------------------------------------------------------

    def generate(self):
        '''To generate the memory-alignment checking code'''

        # initialize the original indentation and the extra indentation
        indent = self.indent
        extra_indent = '  ' 

        # generate the disjoint pragma
        s = '\n'
        s += indent + '#pragma disjoint ('
        s += ','.join(map(lambda v: '*' + self.__printAddress(v), self.vars))
        s += ') \n'

        # generate the alignment test
        s += indent + 'if ((('
        s += '|'.join(map(lambda v: '(int)(' + self.__printAddress(v) + ')', self.vars))
        s += ') & 0xF) == 0) {\n'

        # generate a sequence of alignment intrinsics
        for v in self.vars:
            s += indent + extra_indent + '__alignx(16,' + self.__printAddress(v) + '); \n'

        # append the annotation body code
        s += self.annot_body_code.replace('\n', '\n' + extra_indent) + '\n'

        # generate the unaligned version
        s += indent + '} else {\n'
        s += self.annot_body_code.replace('\n', '\n' + extra_indent) + '\n'
        s += indent + '} \n'
        s += indent

        # return the generated code
        return s
    
