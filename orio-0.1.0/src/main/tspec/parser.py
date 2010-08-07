# 
# The parser for the TSpec (Tuning Specifier)
#

import re, sys
from tool.ZestyParser import *

#--------------------------------------------------------------------------------

# a callback function for a single token
def f_token(p, r, c):
    line_no, col_no = p.coord(c)
    return (r, line_no)

#--------------------------------------------------------------------------------

# keywords
ARG         = Token('arg')
CONSTRAINT  = Token('constraint')
DECL        = Token('decl')
DEF         = Token('def')
LET         = Token('let')
PARAM       = Token('param')
SPEC        = Token('spec')

# tokens
ID          = Token('[A-Za-z_][A-Za-z0-9_]*', group=0) >> f_token
EQUALS      = Token('=')
LBRACKET    = Token('\[')
RBRACKET    = Token('\]')
LBRACE      = Token('\{')
RBRACE      = Token('\}')
SEMI        = Token(';')

# Python expression
PYEXP_SM    = Token('[^;]+', group=0) >> f_token
PYEXP_BR    = Token('[^\]]+', group=0) >> f_token

# ignored strings (i.e. comments and whitespaces)
SPACE       = Token('\s+')

#--------------------------------------------------------------------------------

# let statement
def f_let(p, r, c):
    line_no, col_no = p.coord(c)
    return ('let', line_no, r[0], r[1])

p_let = ((Omit(LET) + Omit(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
          (Omit(EQUALS) ^ 'expected equal sign') + Skip(SPACE) + PYEXP_SM + Skip(SPACE) +
          (Omit(SEMI) ^ 'expected semicolon'))
         >> f_let)

#--------------------------------------------------------------------------------

# argument statement
def f_arg(p, r, c):
    line_no, col_no = p.coord(c)
    return ('arg', line_no, r[0], r[1])

p_arg = ((Omit(ARG) + Omit(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
          (Omit(EQUALS) ^ 'expected equal sign') + Skip(SPACE) + PYEXP_SM + Skip(SPACE) +
          (Omit(SEMI) ^ 'expected semicolon'))
         >> f_arg)

#--------------------------------------------------------------------------------

# parameter statement
def f_param(p, r, c):
    is_range = len(r[1]) > 0
    line_no, col_no = p.coord(c)
    return ('param', line_no, r[0], is_range, r[2])    

p_param = ((Omit(PARAM) + Omit(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
            ((LBRACKET + Skip(SPACE) + (RBRACKET ^ 'expected closing bracket')) | EmptyToken) +
            Skip(SPACE) + (Omit(EQUALS) ^ 'expected equal sign') + Skip(SPACE) +
            PYEXP_SM + Skip(SPACE) + (Omit(SEMI) ^ 'expected semicolon'))
           >> f_param)

#--------------------------------------------------------------------------------

# constraint statement
def f_constraint(p, r, c):
    line_no, col_no = p.coord(c)
    return ('constraint', line_no, r[0], r[1])

p_constraint = ((Omit(CONSTRAINT) + Omit(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
                 (Omit(EQUALS) ^ 'expected equal sign') + Skip(SPACE) + PYEXP_SM + Skip(SPACE) +
                 (Omit(SEMI) ^ 'expected semicolon'))
                >> f_constraint)

#--------------------------------------------------------------------------------

# declaration statement
def f_decl(p, r, c):
    id_name = r[0][-1]
    types = r[0][:-1]
    line_no, col_no = p.coord(c)
    return ('decl', line_no, id_name, types, r[1], r[2])

p_decl = ((Omit(DECL) + Omit(SPACE) +
           TokenSeries(ID, skip=SPACE, min=1) +
           TokenSeries(Omit(LBRACKET) + Skip(SPACE) + Only(PYEXP_BR) + Skip(SPACE) +
                       (Omit(RBRACKET) ^ 'expected closing bracket'), skip=SPACE) +
           ((SEMI >> (lambda r: (None, None))) |
            ((Omit(EQUALS) + Skip(SPACE) + PYEXP_SM + Skip(SPACE) + Omit(SEMI)) >> (lambda r: r[0]))))
          >> f_decl)

def f_decl2(p, r, c):
    id_name = r[0][-1]
    types = r[0][:-1]
    line_no, col_no = p.coord(c)
    return ('decl', line_no, id_name, types, r[1], [])

p_decl2 = ((Omit(DECL) + Omit(SPACE) +
            TokenSeries(ID, skip=SPACE, min=1) +
            TokenSeries(Omit(LBRACKET) + Skip(SPACE) + Only(PYEXP_BR) + Skip(SPACE) +
                        (Omit(RBRACKET) ^ 'expected closing bracket'), skip=SPACE) +
            Skip(SPACE) + Omit(SEMI))
           >> f_decl2)

#--------------------------------------------------------------------------------

# definition statement
def f_def(p, r, c):
    line_no, col_no = p.coord(c)
    return ('def', line_no, r[0], r[1])

p_def = ((Omit(DEF) + Omit(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
          (Omit(LBRACE) ^ 'expected opening brace') +
          TokenSeries((p_let | p_arg | p_param | p_constraint | p_decl | p_decl2), skip=SPACE) +
          (Omit(RBRACE) ^ 'expected closing brace'))
         >> f_def)

#--------------------------------------------------------------------------------

# specification statement
def f_spec(p, r, c):
    line_no, col_no = p.coord(c)
    return ('spec', line_no, r[0], r[1])

p_spec = ((Omit(SPEC) + Omit(SPACE) + (ID ^ 'expected identifier') + Skip(SPACE) +
           (Omit(LBRACE) ^ 'expected opening brace') +
           TokenSeries((p_let | p_def), skip=SPACE) +
           (Omit(RBRACE) ^ 'expected closing brace'))
          >> f_spec)

#--------------------------------------------------------------------------------

# specification body
def f_spec_body(p, r, c):
    return r[0]

p_spec_body = (TokenSeries((p_let | p_def), skip=SPACE) +
               (EOF ^ 'unrecognized statement (must be a let/def statement)')
               >> f_spec_body)

#--------------------------------------------------------------------------------

# program body
def f_program_body(p, r, c):
    return r[0]

p_program_body = (TokenSeries((p_let | p_spec), skip=SPACE) +
                  (EOF ^ 'unrecognized top-level statement (must be a let/spec statement)')
                  >> f_program_body)

#--------------------------------------------------------------------------------

class TSpecParser:
    '''The parser of the TSpec language'''

    def __init__(self):
        '''To instantiate a TSpec parser'''
        pass
        
    #----------------------------------------------------------------------------

    def __parse(self, code, line_no, token):
        '''To parse the given code and return a sequence of statements'''

        # append multiple newlines to imitate the actual line number
        code = ('\n' * (line_no-1)) + code

        # append a newline on the given code
        code += '\n'

        # remove all comments
        code = re.sub(r'#.*?\n', '\n', code)

        # create the parser
        p = ZestyParser(code)
        
        # parse the tuning specifications
        try:
            stmt_seq = p.scan(token)
        except ParseError, e:
            print 'error: %s' % e
            sys.exit(1)

        # return the statement sequence
        return stmt_seq


    #----------------------------------------------------------------------------

    def parseProgram(self, code, line_no = 1):
        '''To parse the given program body and return a sequence of statements'''
        return self.__parse(code, line_no, p_program_body)

    #----------------------------------------------------------------------------

    def parseSpec(self, code, line_no = 1):
        '''To parse the given specification body and return a sequence of statements'''
        return self.__parse(code, line_no, p_spec_body)



