#
# The abstract class for transformation submodule
#

class SubModule:
    '''Transformation submodule'''

    def __init__(self, perf_params, transf_args, stmt):
        '''
        To instantiate a transformation submodule used to transform the annotated code.
        
        The class variables consist of the following:
           perf_params        a table/mapping that maps each performance parameter to its value
           transf_args        a list of transformation arguments
           stmt               the statement AST to be transformed
        '''

        self.perf_params = perf_params
        self.transf_args = transf_args
        self.stmt = stmt

    #----------------------------------------------------------------
    
    def transform(self):
        '''
        The transformation procedure that is used to transform the transformation statement.
        The returned value is the transformed statement AST.
        '''

        raise NotImplementedError('%s: unimplemented abstract function "transform"' %
                                  (self.__class__.__name__))
