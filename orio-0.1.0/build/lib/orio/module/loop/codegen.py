#
# The code generator (i.e. unparser) for the AST classes
#

import sys
import ast

#-------------------------------------------------

class CodeGen:
    '''The code generator for the AST classes'''

    def __init__(self):
        '''To instantiate a code generator'''
        pass

    #----------------------------------------------

    def generate(self, tnode, indent = '  ', extra_indent = '  '):
        '''To generate code that corresponds to the given AST'''

        s = ''

        if isinstance(tnode, ast.NumLitExp):
            s += str(tnode.val)

        elif isinstance(tnode, ast.StringLitExp):
            s += str(tnode.val)

        elif isinstance(tnode, ast.IdentExp):
            s += str(tnode.name)

        elif isinstance(tnode, ast.ArrayRefExp):
            s += self.generate(tnode.exp, indent, extra_indent)
            s += '[' + self.generate(tnode.sub_exp, indent, extra_indent) + ']'

        elif isinstance(tnode, ast.FunCallExp):
            s += self.generate(tnode.exp, indent, extra_indent) + '('
            s += ','.join(map(lambda x: self.generate(x, indent, extra_indent), tnode.args))
            s += ')'

        elif isinstance(tnode, ast.UnaryExp):
            s = self.generate(tnode.exp, indent, extra_indent)
            if tnode.op_type == tnode.PLUS:
                s = '+' + s
            elif tnode.op_type == tnode.MINUS:
                s = '-' + s
            elif tnode.op_type == tnode.LNOT:
                s = '!' + s
            elif tnode.op_type == tnode.PRE_INC:
                s = ' ++' + s
            elif tnode.op_type == tnode.PRE_DEC:
                s = ' --' + s
            elif tnode.op_type == tnode.POST_INC:
                s = s + '++ '
            elif tnode.op_type == tnode.POST_DEC:
                s = s + '-- '
            else:
                print 'internal error: unknown unary operator type: %s' % tnode.op_type
                sys.exit(1)

        elif isinstance(tnode, ast.BinOpExp):
            s += self.generate(tnode.lhs, indent, extra_indent)
            if tnode.op_type == tnode.MUL:
                s += '*'
            elif tnode.op_type == tnode.DIV:
                s += '/'
            elif tnode.op_type == tnode.MOD:
                s += '%'
            elif tnode.op_type == tnode.ADD:
                s += '+'
            elif tnode.op_type == tnode.SUB:
                s += '-'
            elif tnode.op_type == tnode.LT:
                s += '<'
            elif tnode.op_type == tnode.GT:
                s += '>'
            elif tnode.op_type == tnode.LE:
                s += '<='
            elif tnode.op_type == tnode.GE:
                s += '>='
            elif tnode.op_type == tnode.EQ:
                s += '=='
            elif tnode.op_type == tnode.NE:
                s += '!='
            elif tnode.op_type == tnode.LOR:
                s += '||'
            elif tnode.op_type == tnode.LAND:
                s += '&&'
            elif tnode.op_type == tnode.COMMA:
                s += ','
            elif tnode.op_type == tnode.EQ_ASGN:
                s += '='
            else:
                print 'internal error: unknown binary operator type: %s' % tnode.op_type
                sys.exit(1)
            s += self.generate(tnode.rhs, indent, extra_indent)

        elif isinstance(tnode, ast.ParenthExp):
            s += '(' + self.generate(tnode.exp, indent, extra_indent) + ')'

        elif isinstance(tnode, ast.ExpStmt):
            s += indent
            if tnode.exp:
                s += self.generate(tnode.exp, indent, extra_indent)
            s += ';\n'

        elif isinstance(tnode, ast.CompStmt):
            s += indent + '{\n'
            for stmt in tnode.stmts:
                s += self.generate(stmt, indent + extra_indent, extra_indent)
            s += indent + '}\n'

        elif isinstance(tnode, ast.IfStmt):
            s += indent + 'if (' + self.generate(tnode.test, indent, extra_indent) + ') '
            if isinstance(tnode.true_stmt, ast.CompStmt):
                tstmt_s = self.generate(tnode.true_stmt, indent, extra_indent)
                s += tstmt_s[tstmt_s.index('{'):]
                if tnode.false_stmt:
                    s = s[:-1] + ' else '
            else:
                s += '\n'
                s += self.generate(tnode.true_stmt, indent + extra_indent, extra_indent)
                if tnode.false_stmt:
                    s += indent + 'else '
            if tnode.false_stmt:
                if isinstance(tnode.false_stmt, ast.CompStmt):
                    tstmt_s = self.generate(tnode.false_stmt, indent, extra_indent)
                    s += tstmt_s[tstmt_s.index('{'):]
                else:
                    s += '\n'
                    s += self.generate(tnode.false_stmt, indent + extra_indent, extra_indent)

        elif isinstance(tnode, ast.ForStmt):
            s += indent + 'for ('
            if tnode.init:
                s += self.generate(tnode.init, indent, extra_indent)
            s += '; '
            if tnode.test:
                s += self.generate(tnode.test, indent, extra_indent)
            s += '; '
            if tnode.iter:
                s += self.generate(tnode.iter, indent, extra_indent)
            s += ') '
            if isinstance(tnode.stmt, ast.CompStmt): 
                stmt_s = self.generate(tnode.stmt, indent, extra_indent)
                s += stmt_s[stmt_s.index('{'):]
            else:
                s += '\n'
                s += self.generate(tnode.stmt, indent + extra_indent, extra_indent)

        elif isinstance(tnode, ast.TransformStmt):
            print 'internal error: a transformation statement is never generated as an output'
            sys.exit()

        elif isinstance(tnode, ast.VarDecl):
            s += indent + str(tnode.type_name) + ' '
            s += ', '.join(tnode.var_names)
            s += ';\n'

        elif isinstance(tnode, ast.Pragma):
            s += '#pragma ' + str(tnode.pstring) + '\n'

        elif isinstance(tnode, ast.Container):
            s += self.generate(tnode.ast, indent, extra_indent)

        else:
            print 'internal error: unrecognized type of AST: %s' % tnode.__class__.__name__
            sys.exit(1)

        return s


