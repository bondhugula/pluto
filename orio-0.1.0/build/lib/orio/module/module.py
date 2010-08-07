#
# File: src/module/module.py
#

class Module:
    '''The abstract class of Orio's code transformation module'''
    
    def __init__(self, perf_params, module_body_code, annot_body_code,
                 cmd_line_opts, line_no, indent_size):
        '''
        The class constructor used to instantiate a program transformation module.
        
        The following are the class attributes:
          perf_params         a table that maps each performance parameter to its value
          module_body_code    the code inside the module body block
          annot_body_code     the code contained in the annotation body block
          cmd_line_opts       information about the command line options
                              (see src/main/cmd_line_opts.py for more details)
          line_no             the starting line position of the module code in the source code
          indent_size         an integer representing the number of whitespace characters that
                              preceed the leader annotation
        '''

        self.perf_params = perf_params
        self.module_body_code = module_body_code
        self.annot_body_code = annot_body_code
        self.cmd_line_opts = cmd_line_opts
        self.line_no = line_no
        self.indent_size = indent_size
        
        # a boolean value to indicate if the results of the running transformation need to be shown
        self.verbose = self.cmd_line_opts.verbose

    #--------------------------------------------------------

    def transform(self):
        '''
        The main code transformation procedure. The returned value is a string value that
        represents the transformed/optimized code.
        '''

        raise NotImplementedError('%s: unimplemented abstract function "transform"' %
                                  (self.__class__.__name__)) 

