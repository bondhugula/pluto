#
# Scalar replacement transformation
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class ScalarReplace(module.loop.submodule.submodule.SubModule):
    '''The scalar replacement transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        'To instantiate a scalar replacement transformation submodule'
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)
        
    #-----------------------------------------------------------------
    
    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        DTYPE = 'dtype'
        PREFIX = 'prefix'

        # all expected transformation arguments
        dtype = None
        prefix = None
        
        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:

            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
                
            # data type
            if aname in DTYPE:
                dtype = (rhs, line_no)

            # prefix name for scalars
            if aname in PREFIX:
                prefix = (rhs, line_no)

            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check semantics of the transformation arguments
        dtype, prefix = self.checkTransfArgs(dtype, prefix) 

        # return information about the transformation arguments
        return (dtype, prefix)

    #-----------------------------------------------------------------

    def checkTransfArgs(self, dtype, prefix):
        '''Check the semantics of the given transformation arguments'''
                
        # evaluate the data type
        if dtype != None:
            rhs, line_no = dtype
            if dtype != None and not isinstance(rhs, str):
                print 'error:%s: data type argument must be a string: %s' % (line_no, rhs)
                sys.exit(1)
            dtype = rhs
        
        # evaluate the prefix name for scalars variables
        if prefix != None:
            rhs, line_no = prefix
            if rhs != None and not isinstance(rhs, str):
                print 'error:%s: the prefix name of scalars must be a string: %s' % (line_no, rhs)
                sys.exit(1)
            prefix = rhs
            
        # return information about the transformation arguments
        return (dtype, prefix)

    #-----------------------------------------------------------------

    def replaceScalars(self, dtype, prefix, stmt):
        '''To apply scalar replacement transformation'''
        
        # perform the scalar replacement transformation
        t = transformator.Transformator(dtype, prefix, stmt)
        transformed_stmt = t.transform()
        
        # return the transformed statement
        return transformed_stmt
    
    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        dtype, prefix = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the bound replacement transformation
        transformed_stmt = self.replaceScalars(dtype, prefix, self.stmt)

        # return the transformed statement
        return transformed_stmt



    
