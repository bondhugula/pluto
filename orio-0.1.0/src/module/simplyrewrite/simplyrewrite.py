#
# This module does not perform any significant code transformations at all.
# The purpose of this module is merely to provide a simple example to application developers
# for how to extend Orio with a new program transformation module.
#

import module.module

#-----------------------------------------

class SimplyRewrite(module.module.Module):
    '''A simple rewriting module'''

    def __init__(self, perf_params, module_body_code, annot_body_code, cmd_line_opts, line_no, indent_size):
        '''To instantiate a simple rewriting module'''

        module.module.Module.__init__(self, perf_params, module_body_code, annot_body_code,
                                      cmd_line_opts, line_no, indent_size)

    #---------------------------------------------------------------------
    
    def transform(self):
        '''To simply rewrite the annotated code'''

        # to create a comment containing information about the class attributes
        comment = '''
        /*
         perf_params = %s
         module_body_code = "%s"
         annot_body_code = "%s"
         line_no = %s
         indent_size = %s
        */
        ''' % (self.perf_params, self.module_body_code, self.annot_body_code, self.line_no, self.indent_size)

        # to rewrite the annotated code, with the class-attribute comment being prepended
        output_code = comment + self.annot_body_code

        # return the output code
        return output_code

