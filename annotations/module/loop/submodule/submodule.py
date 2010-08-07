#
# The abstract class for annotation submodule.
#

import sys

class SubModule:
    def __init__(self):
        '''To instantiate an annotation submodule'''
        pass

    def transform(self, trans_stmt):
        '''
        The transformation procedure that is used to transform the "transform" statement.

        The input parameters:
          trans_stmt           the "transform" statement to be transformed

        The returned value: the transformed statement
        '''

        raise NotImplementedError('%s: unimplemented abstract function "transform"' %
                                  (self.__class__.__name__))
