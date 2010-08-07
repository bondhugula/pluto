#
# The class for various loop optimization module
#

import sys
import module.module, parser, transformator

#-----------------------------------------

class Loop(module.module.Module):
    
    def __init__(self):
        '''To instantiate a loop optimization module'''
        module.module.Module.__init__(self)

    def transform(self, leader_annot_info, annot_body_code, trailer_annot_code, lang):
        '''To transform the annotated code region'''

        # extract leader annotation information
        leader_annot_code, indent, line_no, module_name, module_body = leader_annot_info

        # parse the module body
        p = parser.getParser(line_no)
        stmts = p.parse(module_body)

        # apply transformation procedure on each statement
        transformed_stmts = transformator.transform(stmts)
        
        # unparse
        transformed_code = ''
        extra_indent = '  '
        include_orig_loop = False
        if include_orig_loop:
            transformed_code += '\n' + indent + '#if ORIG_LOOP' + '\n'
            transformed_code += annot_body_code.replace('\n', '\n' + extra_indent)
            transformed_code += '\n' + indent + '#else ' + '\n\n'
            for s in transformed_stmts:
                transformed_code += s.unparseToC(indent + extra_indent, extra_indent)
            transformed_code += '\n' + indent + '#endif ' + '\n' + indent
        else:
            for s in transformed_stmts:
                transformed_code += s.unparseToC(indent, extra_indent)

        return leader_annot_code + transformed_code + trailer_annot_code

