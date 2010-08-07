#
# Loop transformation submodule that enables pragma directive insertions.
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class Pragma(module.loop.submodule.submodule.SubModule):
    '''The pragma directive insertion submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate a pragma insertion submodule'''

        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

    #-----------------------------------------------------------------

    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        PRAGMAS = 'pragmas'

        # all expected transformation arguments
        pragmas = []

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:

            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)

            # pragma directives
            if aname == PRAGMAS:
                pragmas = (rhs, line_no)

            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check semantics of the transformation arguments
        pragmas, = self.checkTransfArgs(pragmas)
        
        # return information about the transformation arguments
        return (pragmas, )

    #-----------------------------------------------------------------

    def checkTransfArgs(self, pragmas):
        '''Check the semantics of the given transformation arguments'''

        # evaluate the pragma directives
        rhs, line_no = pragmas
        if isinstance(rhs, str):
            pragmas = [rhs]
        else:
            if ((not isinstance(rhs, list) and not isinstance(rhs, tuple)) or
                not reduce(lambda x,y: x and y, map(lambda x: isinstance(x, str), rhs), True)):
                print ('error:%s: pragma directives must be a list/tuple of strings: %s' %
                       (line_no, rhs))
                sys.exit(1)
            pragmas = rhs
        
        # return information about the transformation arguments
        return (pragmas, )

    #-----------------------------------------------------------------

    def insertPragmas(self, pragmas, stmt):
        '''To apply pragma directive insertion'''

        # perform the pragma directive insertion
        t = transformator.Transformator(pragmas, stmt)
        transformed_stmt = t.transform()

        # return the transformed statement
        return transformed_stmt

    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        pragmas, = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the pragma directive insertion 
        transformed_stmt = self.insertPragmas(pragmas, self.stmt)

        # return the transformed statement
        return transformed_stmt



    
