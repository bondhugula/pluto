#
# Loop tiling transformation
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class Tile(module.loop.submodule.submodule.SubModule):
    '''The loop tiling transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate a loop tiling transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)
        
    #-----------------------------------------------------------------
    
    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        TSIZE = 'tsize'
        TINDEX = 'tindex'
        
        # all expected transformation arguments
        tsize = None
        tindex = None

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:

            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)

            # tile size
            if aname == TSIZE:
                tsize = (rhs, line_no)
                
            # tile loop index name
            elif aname == TINDEX:
                tindex = (rhs, line_no)

            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check for undefined transformation arguments
        if tsize == None:
            print 'error:%s: missing tile size argument' % self.__class__.__name__
            sys.exit(1)
        if tindex == None:
            print 'error:%s: missing tile loop index name argument' % self.__class__.__name__
            sys.exit(1)

        # check semantics of the transformation arguments
        tsize, tindex = self.checkTransfArgs(tsize, tindex)
        
        # return information about the transformation arguments
        return (tsize, tindex)

    #-----------------------------------------------------------------

    def checkTransfArgs(self, tsize, tindex):
        '''Check the semantics of the given transformation arguments'''
    
        # evaluate the tile size
        rhs, line_no = tsize
        if not isinstance(rhs, int) or rhs <= 0:
            print 'error:%s: tile size must be a positive integer: %s' % (line_no, rhs)
            sys.exit(1)
        tsize = rhs
            
        # evaluate the tile loop index name
        rhs, line_no = tindex
        if not isinstance(rhs, str):
            print 'error:%s: tile loop index name must be a string: %s' % (line_no, rhs)
            sys.exit(1)
        tindex = rhs

        # return information about the transformation arguments
        return (tsize, tindex)
                
    #-----------------------------------------------------------------

    def tile(self, tsize, tindex, stmt):
        '''To apply loop tiling transformation'''
        
        # perform the loop tiling transformation
        t = transformator.Transformator(tsize, tindex, stmt)
        transformed_stmt = t.transform()
        
        # return the transformed statement
        return transformed_stmt
    
    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        tsize, tindex = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the loop tiling transformation
        transformed_stmt = self.tile(tsize, tindex, self.stmt)

        # return the transformed statement
        return transformed_stmt



    
