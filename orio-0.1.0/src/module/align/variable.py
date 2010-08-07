#
# Variable class denoting the variable instances that need to be checked for alignment
#

import sys

#-----------------------------------------------

class Variable:
    '''The variable that needs to be check for alignment'''

    #------------------------------------------
    
    def __init__(self, vname, dims, line_no):
        '''To instantiate a variable instance'''
        
        self.vname = vname
        self.dims = dims
        self.line_no = line_no

    #------------------------------------------
    
    def __str__(self):
        '''Return the string representation of this variable'''
        return repr(self)

    def __repr__(self):
        '''Return the string representation of this variable'''
        s = ''
        s += str(self.vname)
        if len(self.dims) > 0:
            s += '[' + ']['.join(map(str, self.dims)) + ']'
        return s
    
    #------------------------------------------

    def semantCheck(self):
        '''Perform a semantic check'''

        # check if the variable has no dimensions at all
        if len(self.dims) == 0:
            print 'error:%s: at least one dimension must be defined: "%s"' % (self.line_no, self)
            sys.exit(1)
        
        # one of the dimensions must contain nothing
        if self.dims.count(None) != 1:
            print 'error:%s: there must be one empty bracket: "%s"' % (self.line_no, self)
            sys.exit(1)

        # the last the dimension must be empty
        if self.dims[-1] != None:
            print (('error:%s: the last dimension must be an empty bracket ' + 
                    '(due to the row-major array allocation): "%s"') % (self.line_no, self))
            sys.exit(1)

