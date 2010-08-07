#
# Contain syntax and grammar specifications
#

import sys
import variable, tool.ply.lex, tool.ply.yacc 

#------------------------------------------------

__start_line_no = 1

#------------------------------------------------

# reserved words
reserved = []

tokens = reserved + [
    # literals (identifier, integer constant)
    'ID', 'ICONST',

    # delimeters [ ] , 
    'LBRACKET', 'RBRACKET',
    'COMMA'
    ]

# tokens
t_LBRACKET         = r'\['
t_RBRACKET         = r'\]'
t_COMMA            = r','

# ignored characters
t_ignore = ' \t'

# newlines
def t_NEWLINE(t):
    r'\n+'
    t.lineno += t.value.count('\n')

# identifiers and reserved words
reserved_map = {}
for r in reserved:
    reserved_map[r.lower()] = r

def t_ID(t):
    r'[A-Za-z_]\w*'
    t.type = reserved_map.get(t.value,'ID')
    return t

# integer literal
t_ICONST = r'\d+'

# syntactical error
def t_error(t):
    print 'error:%s: syntactical error: "%s"' % ((t.lineno + __start_line_no - 1),
                                                 t.value[0])
    sys.exit(1)
    
#------------------------------------------------

# grammar
def p_annote(p):
    'annote : var_list_opt'
    p[0] = p[1]

def p_var_list_opt_1(p):
    'var_list_opt :'
    p[0] = []

def p_var_list_opt_2(p):
    'var_list_opt : var_list'
    p[0] = p[1]

def p_var_list_1(p):
    'var_list : var'
    p[0] = [p[1]]

def p_var_list_2(p):
    'var_list : var_list COMMA var'
    p[1].append(p[3])
    p[0] = p[1]

def p_var(p):
    'var : ID dimensions'
    p[0] = variable.Variable(p[1], p[2], p.lineno(1) + __start_line_no - 1)

def p_dimensions_1(p):
    'dimensions : LBRACKET RBRACKET'
    p[0] = [None]

def p_dimensions_2(p):
    'dimensions : LBRACKET index RBRACKET'
    p[0] = [p[2]]

def p_dimensions_3(p):
    'dimensions : dimensions LBRACKET RBRACKET'
    p[1].append(None)
    p[0] = p[1]

def p_dimensions_4(p):
    'dimensions : dimensions LBRACKET index RBRACKET'
    p[1].append(p[3])
    p[0] = p[1]

def p_index_1(p):
    'index : ID'
    p[0] = p[1]

def p_index_2(p):
    'index : ICONST'
    p[0] = p[1]

# grammatical error
def p_error(p):
    print 'error:%s: grammatical error: "%s"' % ((p.lineno + __start_line_no - 1), p.value)
    sys.exit(1)

#------------------------------------------------

def getParser(start_line_no):
    '''Create the parser'''

    # set the starting line number
    global __start_line_no 
    __start_line_no = start_line_no

    # create the lexer and parser
    lexer = tool.ply.lex.lex()
    parser = tool.ply.yacc.yacc(method='LALR', debug=0)

    # return the parser
    return parser

