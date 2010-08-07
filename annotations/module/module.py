#
# The abstract class for annotation module.
#

class Module:
    def __init__(self):
        '''To instantiate an annotation module'''
        pass
    
    def transform(self, leader_info, annot_body_code, trailer_code, lang):
        '''
        The transformation procedure that is used to transform the annotated code
        region.

        The input parameters:
          leader_info          a list that contains information about the leader annotation
             = (code, indent, line_no, module_name, module_body)
             where:
                code           the code of the leader annotation
                indent         the indentation preceeding the leader annotation
                line_no        the line number of the leader annotation is located in
                               the source file (for debugging purposes)
                module_name    the name of the annotation module
                module_body    the code of the leader annotation
          annot_body_code      the code of the annotation body
          trailer_code         the code of the trailer annotation
          lang                 the language of the source code (see "src/main")

        The returned value: the transformed code (in string)
        '''
        
        raise NotImplementedError('%s: unimplemented abstract function "transform"' %
                                  (self.__class__.__name__))
