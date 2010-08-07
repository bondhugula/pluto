#
# Loop transformation submodule that implements a general approach of register tiling (also
# called loop unrolling-and-jamming), which can handle loop transformation in non-rectangular
# iteration spaces.
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class RegTile(module.loop.submodule.submodule.SubModule):
    '''The register tiling transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate a register tiling transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

    #-----------------------------------------------------------------
    
    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        LOOPS = 'loops'
        UFACTORS = 'ufactors'

        # all expected transformation arguments
        loops = None
        ufactors = None

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:

            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
                
            # unroll factors
            if aname == LOOPS:
                loops = (rhs, line_no)
    
            # unroll factors
            elif aname == UFACTORS:
                ufactors = (rhs, line_no)
    
            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check for undefined transformation arguments
        if loops == None:
            print 'error:%s: missing loops argument' % self.__class__.__name__
            sys.exit(1)
        if ufactors == None:
            print 'error:%s: missing unroll factors argument' % self.__class__.__name__
            sys.exit(1)

        # check semantics of the transformation arguments
        loops, ufactors = self.checkTransfArgs(loops, ufactors)
                
        # return information about the transformation arguments
        return (loops, ufactors)

    #-----------------------------------------------------------------

    def checkTransfArgs(self, loops, ufactors):
        '''Check the semantics of the given transformation arguments'''
        
        # evaluate the unroll factors
        rhs, line_no = loops
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: loops value must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        for e in rhs:
            if not isinstance(e, str):
                print 'error:%s: loops element must be a string, found: %s' % (line_no, e)
                sys.exit(1)
        for e in rhs:
            if rhs.count(e) > 1:
                print 'error:%s: loops value contains duplication: "%s"' % (line_no, e)
                sys.exit(1)
        loops = rhs

        # evaluate the unroll factors
        rhs, line_no = ufactors
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: unroll factors value must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        for e in rhs:
            if not isinstance(e, int) or e <= 0:
                print 'error:%s: unroll factor must be a positive integer, found: %s' % (line_no, e)
                sys.exit(1)
        ufactors = rhs

        # compare the loops and unroll factors
        if len(loops) != len(ufactors):
            print 'error:%s: mismatch on the number of loops and unroll factors' % line_no
            sys.exit(1)

        # return information about the transformation arguments
        return (loops, ufactors)

    #-----------------------------------------------------------------

    def tileForRegs(self, loops, ufactors, stmt):
        '''To apply register tiling transformation'''

        # perform the register tiling transformation
        t = transformator.Transformator(loops, ufactors, stmt)
        transformed_stmt = t.transform()

        # return the transformed statement
        return transformed_stmt

    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        loops, ufactors = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the register tiling transformation
        transformed_stmt = self.tileForRegs(loops, ufactors, self.stmt)
        
        # return the transformed statement
        return transformed_stmt



    
