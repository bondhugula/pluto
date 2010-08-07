#
# The class for alignment annotation module for Blue Gene/L architecture.
#

import sys
import module.module, src.main
import parser, unparser

#-----------------------------------------

class Align(module.module.Module):

    def __init__(self):
        '''To instantiate an alignment annotation module'''
        module.module.Module.__init__(self)

    def transform(self, leader_annot_info, annot_body_code, trailer_annot_code, lang):
        '''To transform the annotated code region'''

        # extract leader annotation information
        leader_annot_code, indent, line_no, module_name, module_body = leader_annot_info

        # parse the module body
        p = parser.getParser(line_no)
        vars = p.parse(module_body)

        # semantic check
        for v in vars:
            if lang == src.main.C_CPP:
                v.rowMajorCheck()
            elif lang == src.main.FORTRAN:
                v.columnMajorCheck()        
            else:
                print 'internal error: unknown source language'
                sys.exit(1)
        
        # unparse
        if lang == src.main.C_CPP:
            extra_indent = ' ' * 2
            transformed_code = unparser.unparseToC(vars, annot_body_code,
                                                   indent, extra_indent)
        elif lang == src.main.FORTRAN:
            indent = ' ' * 6
            extra_indent = ' ' * 2
            transformed_code = unparser.unparseToFortran(vars, annot_body_code,
                                                         indent, extra_indent)
        else:
            print 'internal error: unknown source language'
            sys.exit(1)

        return leader_annot_code + transformed_code + trailer_annot_code

