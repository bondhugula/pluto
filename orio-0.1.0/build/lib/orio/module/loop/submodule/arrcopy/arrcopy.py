#
# Array copy transformation
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class ArrCopy(module.loop.submodule.submodule.SubModule):
    '''The array copy transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate an array copy transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

    #-----------------------------------------------------------------
    
    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        AREF = 'aref'
        SUFFIX = 'suffix'
        DTYPE = 'dtype'
        DIMS = 'dimsizes'
        
        # all expected transformation arguments
        aref = None
        suffix = None
        dtype = None
        dimsizes = None

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:
            
            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)

            # array reference
            if aname == AREF:
                aref = (rhs, line_no)

            # suffix for the array buffer name
            elif aname == SUFFIX:
                suffix = (rhs, line_no)

            # array elements' data type
            elif aname == DTYPE:
                dtype = (rhs, line_no)

            # array dimension sizes
            elif aname == DIMS:
                dimsizes = (rhs, line_no)

            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check for undefined transformation arguments
        if aref == None:
            print 'error:%s: missing array reference argument' % self.__class__.__name__
            sys.exit(1)
        if dimsizes == None:
            print 'error:%s: missing array dimension sizes argument' % self.__class__.__name__
            sys.exit(1)

        # check semantics of the transformation arguments
        (aref, suffix, dtype, dimsizes) = self.checkTransfArgs(aref, suffix, dtype, dimsizes)

        # return information about the transformation arguments
        return (aref, suffix, dtype, dimsizes)

    #-----------------------------------------------------------------

    def checkTransfArgs(self, aref, suffix, dtype, dimsizes):
        '''Check the semantics of the given transformation arguments'''

        # evaluate the array reference
        rhs, line_no = aref
        if not isinstance(rhs, str):
            print 'error:%s: array reference argument must be a string: %s' % (line_no, rhs)
            sys.exit(1)
        aref = rhs

        # evaluate the suffix of the array buffer name
        if suffix != None:
            rhs, line_no = suffix
            if rhs != None and not isinstance(rhs, str):
                print 'error:%s: suffix argument must be a string: %s' % (line_no, rhs)
                sys.exit(1)
            suffix = rhs

        # evaluate the data type of the array elements
        if dtype != None:
            rhs, line_no = dtype
            if rhs != None and not isinstance(rhs, str):
                print 'error:%s: data type argument must be a string: %s' % (line_no, rhs)
                sys.exit(1)
            dtype = rhs

        # evaluate the data type of the array elements
        rhs, line_no = dimsizes
        if ((not isinstance(rhs, list) and not isinstance(rhs, tuple)) or
            not reduce(lambda x,y: x and y, map(lambda x: isinstance(x, int), rhs), True) or 
            not reduce(lambda x,y: x and y, map(lambda x: x > 0, rhs), True)):
            print (('error:%s: array dimension sizes argument must be a list of positive ' +
                   'integers: %s') % (line_no, rhs))
            sys.exit(1)
        dimsizes = rhs
            
        # return information about the transformation arguments
        return (aref, suffix, dtype, dimsizes)

    #-----------------------------------------------------------------

    def optimizeArrayCopy(self, aref, suffix, dtype, dimsizes, stmt):
        '''To apply array copy optimization'''

        # perform the bound replacement transformation
        t = transformator.Transformator(aref, suffix, dtype, dimsizes, stmt)
        transformed_stmt = t.transform()
        
        # return the transformed statement
        return transformed_stmt
    
    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        aref, suffix, dtype, dimsizes = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the bound replacement transformation
        transformed_stmt = self.optimizeArrayCopy(aref, suffix, dtype, dimsizes, self.stmt)

        # return the transformed statement
        return transformed_stmt



    
