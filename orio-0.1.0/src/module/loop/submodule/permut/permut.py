#
# Loop transformation submodule that implements loop permutation/interchange
#

import sys
import module.loop.submodule.submodule, transformator

#---------------------------------------------------------------------

class Permut(module.loop.submodule.submodule.SubModule):
    '''The loop permutation transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate a loop permutation transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

    #-----------------------------------------------------------------

    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        SEQ = 'seq'

        # all expected transformation arguments
        seq = None

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:

            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
                
            # permutation sequence
            if aname == SEQ:
                seq = (rhs, line_no)
                
            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check for undefined transformation arguments
        if seq == None:
            print 'error:%s: missing permutation sequence argument' % self.__class__.__name__
            sys.exit(1)

        # check semantics of the transformation arguments
        seq, = self.checkTransfArgs(seq)
        
        # return information about the transformation arguments
        return (seq, )

    #-----------------------------------------------------------------

    def checkTransfArgs(self, seq):
        '''Check the semantics of the given transformation arguments'''

        # evaluate the permutation sequence
        rhs, line_no = seq
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print ('error:%s: permutation sequence must be a list/tuple of loop index names: %s' %
                   (line_no, rhs))
            sys.exit(1)
        inames = {}
        for i in rhs:
            if isinstance(i, str):
                pass
            elif isinstance(i, list) and len(i) == 1 and isinstance(i[0], str):
                i = i[0]
            else:
                print ('error:%s: invalid element of the permutation sequence: %s' %
                       (line_no, i))
                sys.exit(1)
            if i in inames:
                print ('error:%s: permutation sequence contains repeated loop index: %s' %
                       (line_no, i))
                sys.exit(1)
            inames[i] = None
        seq = rhs
        
        # return information about the transformation arguments
        return (seq, )

    #-----------------------------------------------------------------

    def permute(self, seq, stmt):
        '''To apply loop permutation/interchange transformation'''

        # perform the loop permutation transformation
        t = transformator.Transformator(seq, stmt)
        transformed_stmt = t.transform()

        # return the transformed statement
        return transformed_stmt

    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        seq, = self.readTransfArgs(self.perf_params, self.transf_args)

        # perform the loop permutation transformation
        transformed_stmt = self.permute(seq, self.stmt)
        
        # return the transformed statement
        return transformed_stmt



    
