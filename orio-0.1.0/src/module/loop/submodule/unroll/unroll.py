#
# Loop transformation submodule that implements a pure loop unrolling
#

import sys
import module.loop.submodule.submodule, module.loop.submodule.unrolljam.unrolljam

#---------------------------------------------------------------------

class Unroll(module.loop.submodule.submodule.SubModule):
    '''The unrolling transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate an unrolling transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

        self.ujam_smod = module.loop.submodule.unrolljam.unrolljam.UnrollJam()
        
    #-----------------------------------------------------------------

    def readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''
        return self.ujam_smod.readTransfArgs(perf_params, transf_args)

    #-----------------------------------------------------------------

    def unroll(self, ufactor, stmt):
        '''To apply unroll-and-jam transformation'''
        return self.ujam_smod.unrollAndJam(ufactor, False, stmt)
    
    #-----------------------------------------------------------------

    def transform(self):
        '''To perform code transformations'''

        # read all transformation arguments
        ufactor, = self.readTransfArgs(self.perf_params, self.transf_args)
        
        # perform the unroll-and-jam transformation
        transformed_stmt = self.unroll(ufactor, self.stmt)
        
        # return the transformed statement
        return transformed_stmt
                                                      
