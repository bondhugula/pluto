#
# Variable class that denotes the variable instances to be checked for alignment.
#

import sys

#-----------------------------------------------

class Variable:
    
    def __init__(self, var_name, dimensions, line_no):
        '''Instantiate the variable to be checked for its alignment'''
        self.var_name = var_name
        self.dimensions = dimensions
        self.line_no = line_no

    def __repr__(self):
        '''Return the string representation of this variable'''
        s = ''
        s += str(self.var_name)
        if len(self.dimensions) > 0:
            s += '['
        for i, d in enumerate(self.dimensions):
            if i > 0:
                s += ']['
            if d:
                s += str(d)
        if len(self.dimensions) > 0:
            s += ']'
        return s

    def __str__(self):
        '''Return the string representation of this variable'''
        return repr(self)

    def rowMajorCheck(self):
        '''Check if the starting address fits with C/C++ row-major array allocation'''

        if self.dimensions.count(None) != 1:
            print 'error:%s: there must be one empty bracket: "%s"' % (self.line_no, self)
            sys.exit(1)
        if self.dimensions[len(self.dimensions) - 1] != None:
            print (('error:%s: the last dimension must be an empty bracket ' + 
                    '(C/C++ row-major array location): "%s"') % (self.line_no, self))
            sys.exit(1)

    def columnMajorCheck(self):
        '''Check if the starting address fits with Fortran column-major array allocation'''

        if self.dimensions.count(None) != 1:
            print 'error:%s: there must be one empty bracket: "%s"' % (self.line_no, self)
            sys.exit(1)
        if self.dimensions[0] != None:
            print (('error:%s: the first dimension must be an empty bracket ' + 
                    '(Fortran column-major array location): "%s"') % (self.line_no, self))
            sys.exit(1)

