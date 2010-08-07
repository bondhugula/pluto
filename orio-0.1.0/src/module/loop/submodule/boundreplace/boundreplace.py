#
# Bound replacement transformation
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class BoundReplace(module.loop.submodule.submodule.SubModule):
    '''The bound replacement transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate a bound replacement transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

    #-----------------------------------------------------------------
    
    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        LPREFIX = 'lprefix'
        UPREFIX = 'uprefix'
        
        # all expected transformation arguments
        lprefix = None
        uprefix = None

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:

            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
                
            # prefix name for lower bound
            if aname == LPREFIX:
                lprefix = (rhs, line_no)
                
            # prefix name for upper bound
            elif aname == UPREFIX:
                uprefix = (rhs, line_no)

            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check semantics of the transformation arguments
        lprefix, uprefix = self.checkTransfArgs(lprefix, uprefix)

        # return information about the transformation arguments
        return (lprefix, uprefix)
    
    #-----------------------------------------------------------------

    def checkTransfArgs(self, lprefix, uprefix):
        '''Check the semantics of the given transformation arguments'''

        # evaluate the prefix name for lower/upper bounds
        for i, prefix in enumerate([lprefix, uprefix]):
            if prefix != None:
                rhs, line_no = prefix
                if rhs != None and not isinstance(rhs, str):
                    print (('error:%s: the prefix name of the lower/upper bound must be ' +
                            'a string: %s') % (line_no, rhs))
                    sys.exit(1)
                if i == 0:
                    lprefix = rhs
                elif i == 1:
                    uprefix = rhs
            
        # return information about the transformation arguments
        return (lprefix, uprefix)

    #-----------------------------------------------------------------

    def replaceBounds(self, lprefix, uprefix, stmt):
        '''To apply bound replacement transformation'''

        # perform the bound replacement transformation
        t = transformator.Transformator(lprefix, uprefix, stmt)
        transformed_stmt = t.transform()
        
        # return the transformed statement
        return transformed_stmt
    
    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        lprefix, uprefix = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the bound replacement transformation
        transformed_stmt = self.replaceBounds(lprefix, uprefix, self.stmt)

        # return the transformed statement
        return transformed_stmt



    
