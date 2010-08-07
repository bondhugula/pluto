#
# The classes of the abstract syntax tree 
#
#  AST 
#   |
#   +-- Exp 
#   |    |
#   |    +-- NumLitExp
#   |    +-- StringLitExp
#   |    +-- IdentExp
#   |    +-- ArrayRefExp 
#   |    +-- FunCallExp 
#   |    +-- UnaryExp 
#   |    +-- BinOpExp 
#   |    +-- ParenthExp
#   |
#   +-- Stmt 
#   |    |
#   |    +-- ExpStmt 
#   |    +-- CompStmt 
#   |    +-- IfStmt 
#   |    +-- ForStmt 
#   |    +-- TransformStmt 
#   |
#   +-- NewAST 
#        |
#        +-- VarDecl 
#
# - The NewAST is used to denote ASTs that are used only in the output code generation.
# - UnaryExp.AND never appears in the input code. It is used in the output code generation only.
#

import sys

#-----------------------------------------------
# AST - the base class
#-----------------------------------------------

class AST:
    def __init__(self, line_no = ''):
        '''Create an abstract syntax tree node'''
        self.line_no = line_no        # may be null

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        raise NotImplementedError('%s: abstract function "replicate" not implemented' %
                                  self.__class__.__name__)
        
    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        raise NotImplementedError('%s: abstract function "__repr__" not implemented' %
                                  self.__class__.__name__)

    def __str__(self):
        '''Return a string representation of this abstract syntax tree node'''
        return repr(self)

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        raise NotImplementedError('%s: abstract function "unparseToC" not implemented' %
                                  self.__class__.__name__)

    def unparseToFortran(self, indent, extra_indent):
        '''Generate Fortran code from this abstract syntax tree node'''
        raise NotImplementedError('%s: abstract function "unparseToFortran" not implemented'
                                  % self.__class__.__name__)

#-----------------------------------------------
# Expression
#-----------------------------------------------

class Exp(AST):
    def __init__(self, line_no = ''):
        '''Create an expression'''
        AST.__init__(self, line_no)

#-----------------------------------------------
# Number Literal
#-----------------------------------------------

class NumLitExp(Exp):
    INT = 1
    FLOAT = 2
    
    def __init__(self, val, lit_type, line_no = ''):
        '''Create a numeric literal'''
        Exp.__init__(self, line_no)
        self.val = val
        self.lit_type = lit_type

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return NumLitExp(self.val, self.lit_type, self.line_no)
        
    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        return str(self.val)

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        return str(self.val)

#-----------------------------------------------
# String Literal
#-----------------------------------------------

class StringLitExp(Exp):
    def __init__(self, val, line_no = ''):
        '''Create a string literal'''
        Exp.__init__(self, line_no)
        self.val = val

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return StringLitExp(self.val, self.line_no)
        
    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        return str(self.val)
    
    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        return str(self.val)

#-----------------------------------------------
# Identifier
#-----------------------------------------------

class IdentExp(Exp):
    def __init__(self, name, line_no = ''):
        '''Create an identifier'''
        Exp.__init__(self, line_no)
        self.name = name

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return IdentExp(self.name, self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        return str(self.name)

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        return str(self.name)

#-----------------------------------------------
# Array Reference
#-----------------------------------------------

class ArrayRefExp(Exp):
    def __init__(self, exp, sub_exp, line_no = ''):
        '''Create an array reference'''
        Exp.__init__(self, line_no)
        self.exp = exp
        self.sub_exp = sub_exp

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return ArrayRefExp(self.exp.replicate(), self.sub_exp.replicate(), self.line_no)
        
    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        return str(self.exp) + '[' + str(self.sub_exp) + ']'
        
    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = ''
        s += self.exp.unparseToC(indent, extra_indent) + '['
        s += self.sub_exp.unparseToC(indent, extra_indent) + ']'
        return s

#-----------------------------------------------
# Function Call
#-----------------------------------------------

class FunCallExp(Exp):
    def __init__(self, exp, args, line_no = ''):
        '''Create a function call'''
        Exp.__init__(self, line_no)
        self.exp = exp
        self.args = args

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return FunCallExp(self.exp.replicate(), [a.replicate() for a in self.args],
                          self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = str(self.exp) + '('
        s += ', '.join(map(str, self.args))
        s += ')'
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = self.exp.unparseToC(indent, extra_indent) + '('
        s += ', '.join(map(lambda x: x.unparseToC(indent, extra_indent), self.args))
        s += ')'
        return s
        
#-----------------------------------------------
# Unary Expression
#-----------------------------------------------

class UnaryExp(Exp):
    PLUS = 1
    MINUS = 2
    LNOT = 3
    PRE_INC = 4
    PRE_DEC = 5
    POST_INC = 6
    POST_DEC = 7
    AND = 8
    
    def __init__(self, exp, op_type, line_no = ''):
        '''Create a unary operation expression'''
        Exp.__init__(self, line_no)
        self.exp = exp
        self.op_type = op_type

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return UnaryExp(self.exp.replicate(), self.op_type, self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = str(self.exp)
        if self.op_type == UnaryExp.PLUS:
            s = '+' + s
        elif self.op_type == UnaryExp.MINUS:
            s = '-' + s
        elif self.op_type == UnaryExp.LNOT:
            s = '!' + s
        elif self.op_type == UnaryExp.PRE_INC:
            s = '++' + s
        elif self.op_type == UnaryExp.PRE_DEC:
            s = '--' + s
        elif self.op_type == UnaryExp.POST_INC:
            s = s + '++'
        elif self.op_type == UnaryExp.POST_DEC:
            s = s + '--'
        elif self.op_type == UnaryExp.AND:
            s = '&' + s
        else:
            print 'internal error: unknown unary operator type'
            sys.exit(1)            
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = self.exp.unparseToC(indent, extra_indent)
        if self.op_type == UnaryExp.PLUS:
            s = '+' + s
        elif self.op_type == UnaryExp.MINUS:
            s = '-' + s
        elif self.op_type == UnaryExp.LNOT:
            s = '!' + s
        elif self.op_type == UnaryExp.PRE_INC:
            s = '++' + s
        elif self.op_type == UnaryExp.PRE_DEC:
            s = '--' + s
        elif self.op_type == UnaryExp.POST_INC:
            s = s + '++'
        elif self.op_type == UnaryExp.POST_DEC:
            s = s + '--'
        elif self.op_type == UnaryExp.AND:
            s = '&' + s
        else:
            print 'internal error: unknown unary operator type'
            sys.exit(1)            
        return s

#-----------------------------------------------
# Binary Operation
#-----------------------------------------------

class BinOpExp(Exp):
    MUL = 1
    DIV = 2
    MOD = 3
    ADD = 4
    SUB = 5
    LT = 6
    GT = 7
    LE = 8
    GE = 9
    EQ = 10
    NE = 11
    LOR = 12
    LAND = 13
    COMMA = 14
    EQ_ASGN = 15

    def __init__(self, lhs, rhs, op_type, line_no = ''):
        '''Create a binary operation expression'''
        Exp.__init__(self, line_no)
        self.lhs = lhs
        self.rhs = rhs
        self.op_type = op_type

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return BinOpExp(self.lhs.replicate(), self.rhs.replicate(), self.op_type, self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = str(self.lhs)
        if (self.op_type == BinOpExp.MUL):
            s += ' * '
        elif (self.op_type == BinOpExp.DIV):
            s += ' / '
        elif (self.op_type == BinOpExp.MOD):
            s += ' % '
        elif (self.op_type == BinOpExp.ADD):
            s += ' + '
        elif (self.op_type == BinOpExp.SUB):
            s += ' - '
        elif (self.op_type == BinOpExp.LT):
            s += ' < '
        elif (self.op_type == BinOpExp.GT):
            s += ' > '
        elif (self.op_type == BinOpExp.LE):
            s += ' <= '
        elif (self.op_type == BinOpExp.GE):
            s += ' >= '
        elif (self.op_type == BinOpExp.EQ):
            s += ' == '
        elif (self.op_type == BinOpExp.NE):
            s += ' != '
        elif (self.op_type == BinOpExp.LOR):
            s += ' || '
        elif (self.op_type == BinOpExp.LAND):
            s += ' && '
        elif (self.op_type == BinOpExp.COMMA):
            s += ' , '
        elif (self.op_type == BinOpExp.EQ_ASGN):
            s += ' = '
        else:
            print 'internal error: unknown bin-op operator type'
            sys.exit(1)
        s += str(self.rhs)
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = self.lhs.unparseToC(indent, extra_indent)
        if (self.op_type == BinOpExp.MUL):
            s += ' * '
        elif (self.op_type == BinOpExp.DIV):
            s += ' / '
        elif (self.op_type == BinOpExp.MOD):
            s += ' % '
        elif (self.op_type == BinOpExp.ADD):
            s += ' + '
        elif (self.op_type == BinOpExp.SUB):
            s += ' - '
        elif (self.op_type == BinOpExp.LT):
            s += ' < '
        elif (self.op_type == BinOpExp.GT):
            s += ' > '
        elif (self.op_type == BinOpExp.LE):
            s += ' <= '
        elif (self.op_type == BinOpExp.GE):
            s += ' >= '
        elif (self.op_type == BinOpExp.EQ):
            s += ' == '
        elif (self.op_type == BinOpExp.NE):
            s += ' != '
        elif (self.op_type == BinOpExp.LOR):
            s += ' || '
        elif (self.op_type == BinOpExp.LAND):
            s += ' && '
        elif (self.op_type == BinOpExp.COMMA):
            s += ' , '
        elif (self.op_type == BinOpExp.EQ_ASGN):
            s += ' = '
        else:
            print 'internal error: unknown bin-op operator type'
            sys.exit(1)
        s += self.rhs.unparseToC(indent, extra_indent)
        return s
        
#-----------------------------------------------
# Parenthesized Expression
#-----------------------------------------------

class ParenthExp(Exp):
    def __init__(self, exp, line_no = ''):
        '''Create a parenthesized expression'''
        Exp.__init__(self, line_no)
        self.exp = exp

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return ParenthExp(self.exp.replicate(), self.line_no)
        
    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        return '(' + str(self.exp) + ')'

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        return '(' + self.exp.unparseToC(indent, extra_indent) + ')'        

#-----------------------------------------------
# Statement
#-----------------------------------------------

class Stmt(AST):
    def __init__(self, line_no = ''):
        '''Create a statement'''
        AST.__init__(self, line_no)

#-----------------------------------------------
# Expression Statement
#-----------------------------------------------

class ExpStmt(Stmt):
    def __init__(self, exp, line_no = ''):
        '''Create an expression statement'''
        Stmt.__init__(self, line_no)
        self.exp = exp      # may be null

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        r_e = self.exp
        if r_e != None:
            r_e = r_e.replicate()
        return ExpStmt(r_e, self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = ''
        if self.exp != None:
            s += str(self.exp)
        s += '; '
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = indent
        if self.exp != None:
            s += self.exp.unparseToC(indent, extra_indent)
        s += '; \n'
        return s
        
#-----------------------------------------------
# Compound Statement
#-----------------------------------------------

class CompStmt(Stmt):
    def __init__(self, stmts, line_no = ''):
        '''Create a compound statement'''
        Stmt.__init__(self, line_no)
        self.stmts = stmts

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return CompStmt([s.replicate() for s in self.stmts], self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = '{'
        for i, t in enumerate(self.stmts):
            if (i > 0):
                s += ' '
            s += str(t)
        s += '}'
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = indent + '{ \n'
        for stmt in self.stmts:
            s += stmt.unparseToC(indent + extra_indent, extra_indent)
        s += indent + '} \n'
        return s

#-----------------------------------------------
# If-Then-Else
#-----------------------------------------------

class IfStmt(Stmt):
    def __init__(self, test, true_stmt, false_stmt = None, line_no = ''):
        '''Create an if statement'''
        Stmt.__init__(self, line_no)
        self.test = test
        self.true_stmt = true_stmt
        self.false_stmt = false_stmt     # may be null

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        f_s = self.false_stmt
        if f_s:
            f_s = f_s.replicate()
        return IfStmt(self.test.replicate(), self.true_stmt.replicate(), f_s, self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = 'if (' + str(self.test) + ') ' + str(self.true_stmt)
        if self.false_stmt:
            s += ' else ' + str(self.false_stmt)
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = indent + 'if (' + self.test.unparseToC(indent, extra_indent) + ') \n'
        s += self.true_stmt.unparseToC(indent + extra_indent, extra_indent)
        if self.false_stmt:
            s += indent + 'else \n'
            s += self.false_stmt.unparseToC(indent + extra_indent, extra_indent)
        return s

#-----------------------------------------------
# For Loop
#-----------------------------------------------

class ForStmt(Stmt):
    def __init__(self, init, test, iter, stmt, line_no = ''):
        '''Create a for-loop statement'''
        Stmt.__init__(self, line_no)
        self.init = init    # may be null
        self.test = test    # may be null
        self.iter = iter    # may be null
        self.stmt = stmt

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        r_in = self.init
        r_t = self.test
        r_it = self.iter
        if r_in:
            r_in = r_in.replicate()
        if r_t:
            r_t = r_t.replicate()
        if r_it:
            r_it = r_it.replicate()
        return ForStmt(r_in, r_t, r_it, self.stmt.replicate(), self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = 'for ('
        if self.init:
            s += str(self.init)
        s += '; '
        if self.test:
            s += str(self.test)
        s += '; '
        if self.iter:
            s += str(self.iter)    
        s += ') ' + str(self.stmt)
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = '\n'
        s += indent + 'for ('
        if self.init:
            s += self.init.unparseToC(indent, extra_indent)
        s += '; '
        if self.test:
            s += self.test.unparseToC(indent, extra_indent)
        s += '; '
        if self.iter:
            s += self.iter.unparseToC(indent, extra_indent)
        s += ') ' + self.stmt.unparseToC(indent + extra_indent, extra_indent)
        return s
        
#-----------------------------------------------
# Transformation
#-----------------------------------------------

class TransformStmt(Stmt):
    def __init__(self, name, kw_args, stmt, line_no = ''):
        '''Create a transformation statement'''
        Stmt.__init__(self, line_no)
        self.name = name
        self.kw_args = kw_args
        self.stmt = stmt

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return TransformStmt(self.name, [k.replicate() for k in self.kw_args],
                             self.stmt.replicate(), self.line_no)

    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = 'transform ' + str(self.name) + ' ('
        for i, k in enumerate(self.kw_args):
            if i > 0:
                s += ', '
            s += str(k)
        s += ') ' + str(self.stmt)
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        print 'internal error: a transformation statement is never generated as output'
        sys.exit(1)

#-----------------------------------------------
# New AST
#-----------------------------------------------

class NewAST(AST):
    def __init__(self, line_no = ''):
        '''Create a newly-added statement'''
        AST.__init__(self, line_no)

#-----------------------------------------------
# Variable Declaration
#-----------------------------------------------

class VarDecl(NewAST):
    def __init__(self, type_name, var_names, line_no = ''):
        '''Create a variable declaration'''
        NewAST.__init__(self, line_no)
        self.type_name = type_name
        self.var_names = var_names

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return VarDecl(self.type_name, self.var_names[:], self.line_no)
    
    def __repr__(self):
        '''Return a string representation of this abstract syntax tree node'''
        s = ''
        s += self.type_name + ' '
        s += ', '.join(map(str, self.var_names))
        s += '; '
        return s

    def unparseToC(self, indent, extra_indent):
        '''Generate C/C++ code from this abstract syntax tree node'''
        s = ''
        s += indent + self.type_name + ' '
        s += ', '.join(map(lambda x: x.unparseToC(indent, extra_indent), self.var_names))
        s += '; \n'
        return s
    
