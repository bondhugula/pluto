# 
# CodeFragment
#  |
#  +-- Ann
#  |    |
#  |    +-- LeaderAnn
#  |    +-- TrailerAnn
#  |
#  +-- NonAnn
#  |
#  +-- AnnCodeRegion
#

#-----------------------------------------

class CodeFragment:

    def __init__(self):
        '''To instantiate a code fragment'''
        pass

#-----------------------------------------

class Ann(CodeFragment):

    def __init__(self, code, line_no, indent_size):
        '''To instantiate an annotation'''

        CodeFragment.__init__(self)
        self.code = code
        self.line_no = line_no
        self.indent_size = indent_size

#-----------------------------------------

class LeaderAnn(Ann):

    def __init__(self, code, line_no, indent_size, mod_name, mod_name_line_no,
                 mod_code, mod_code_line_no):

        Ann.__init__(self, code, line_no, indent_size)
        self.mod_name = mod_name
        self.mod_name_line_no = mod_name_line_no
        self.mod_code = mod_code
        self.mod_code_line_no = mod_code_line_no

#-----------------------------------------

class TrailerAnn(Ann):

    def __init__(self, code, line_no, indent_size):
        '''To instantiate a trailer annotation'''

        Ann.__init__(self, code, line_no, indent_size)

#-----------------------------------------

class NonAnn(CodeFragment):

    def __init__(self, code, line_no, indent_size):
        '''To instantiate a non-annotation'''

        CodeFragment.__init__(self)
        self.code = code
        self.line_no = line_no
        self.indent_size = indent_size

#-----------------------------------------

class AnnCodeRegion(CodeFragment):

    def __init__(self, leader_ann, cfrags, trailer_ann):
        '''To instantiate an annotation code region'''

        CodeFragment.__init__(self)
        self.leader_ann = leader_ann
        self.cfrags = cfrags
        self.trailer_ann = trailer_ann

    
        

