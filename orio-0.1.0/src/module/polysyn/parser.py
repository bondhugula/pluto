# 
# The parser for the PolySyn transformation module
#

import re, sys
from tool.ZestyParser import *

#--------------------------------------------------------------------------------

# a callback function for a single token
def f_token(p, r, c):
    line_no, col_no = p.coord(c)
    line_no += p.start_line_no - 1
    return (r, line_no)

#--------------------------------------------------------------------------------

# tokens
ID          = Token('[A-Za-z_][A-Za-z0-9_]*', group=0) >> f_token
EQUALS      = Token('=')
SEMI        = Token(';')

# Python expression
PYEXP       = Token('[^;]+', group=0) >> f_token

# ignored strings (i.e. whitespaces)
SPACE       = Token('\s+')

#--------------------------------------------------------------------------------

# assignment
def f_assign(p, r, c):
    line_no, col_no = p.coord(c)
    line_no += p.start_line_no - 1
    return (line_no, r[0], r[1])

p_assign = ((Skip(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
             (Omit(EQUALS) ^ 'expected equal sign') + Skip(SPACE) + PYEXP + Skip(SPACE) +
             (Omit(SEMI) ^ 'expected semicolon'))
            >> f_assign)

#--------------------------------------------------------------------------------

# program
def f_program(p, r, c):
    return r[0]

p_program = ((TokenSeries(p_assign, skip=SPACE,
                          until=(Skip(SPACE) + EOF, 'not a valid assignment')) +
              Skip(SPACE) + (EOF ^ 'not a valid assignment')) 
             >> f_program)

#--------------------------------------------------------------------------------

class Parser:
    '''The parser of the PolySyn module'''

    def __init__(self):
        '''To instantiate a parser instance'''
        pass
        
    #----------------------------------------------------------------------------

    def parse(self, code, line_no):
        '''To parse the given code'''

        # remove all comments
        code = code + '\n'
        code = re.sub(r'#.*?\n', '\n', code)

        # if a blank code
        if code.strip() == '':
            return []

        # create the parser
        p = ZestyParser(code)

        # update the starting line number of the code
        p.start_line_no = line_no

        # parse the tuning specifications
        try:
            assigns = p.scan(p_program)
        except ParseError, e:
            print 'error: %s' % e
            sys.exit(1)

        # return the assignment sequence
        return assigns


