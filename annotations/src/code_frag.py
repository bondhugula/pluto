# 
# Below is the class structure of code fragment.
#
# CodeFragment
#  |
#  +-- Annotation
#  |    |
#  |    +-- LeaderAnnotation
#  |    +-- TrailerAnnotation
#  |
#  +-- NonAnnotation
#  |
#  +-- AnnotatedCodeRegion
#  |
#  +-- TransformedCode
#

#-----------------------------------------

class CodeFragment:

    def __init__(self):
        '''Create a code fragment'''
        pass

#-----------------------------------------

class Annotation(CodeFragment):

    def __init__(self, code, indent, line_no):
        '''Create an annotation code fragment'''
        CodeFragment.__init__(self)
        self.code = code
        self.indent = indent
        self.line_no = line_no

#-----------------------------------------

class LeaderAnnotation(Annotation):

    def __init__(self, code, indent, line_no, module_name, module_body):
        '''Create a leader annotation code fragment'''
        Annotation.__init__(self, code, indent, line_no)
        self.module_name = module_name
        self.module_body = module_body

#-----------------------------------------

class TrailerAnnotation(Annotation):

    def __init__(self, code, indent, line_no):
        '''Create a trailer annotation code fragment'''
        Annotation.__init__(self, code, indent, line_no)

#-----------------------------------------

class NonAnnotation(CodeFragment):

    def __init__(self, code, indent, line_no):
        '''Create a non-annotation code fragment'''
        CodeFragment.__init__(self)
        self.code = code
        self.indent = indent
        self.line_no = line_no

#-----------------------------------------

class AnnotatedCodeRegion(CodeFragment):

    def __init__(self, leader_annot, annot_body, trailer_annot):
        '''Create an annotated code region'''
        CodeFragment.__init__(self)
        self.leader_annot = leader_annot
        self.annot_body = annot_body
        self.trailer_annot = trailer_annot

#-----------------------------------------

class TransformedCode(CodeFragment):

    def __init__(self, code):
        '''Create a transformed code'''
        CodeFragment.__init__(self)
        self.code = code



