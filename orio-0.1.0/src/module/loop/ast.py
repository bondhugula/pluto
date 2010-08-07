#
# The classes for the abstract syntax tree (AST)
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
#        +-- Pragma 
#        +-- Container
#
# - The NewAST is an AST used only in the output code generation. Such separation is needed to
#   simplify the input language.
#

import sys
import codegen

#-----------------------------------------------
# AST - Abstract Syntax Tree
#-----------------------------------------------

class AST:

    def __init__(self, line_no = ''):
        '''Create an abstract syntax tree node'''
        self.line_no = line_no           # may be null (i.e. empty string)
        
    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        raise NotImplementedError('%s: abstract function "replicate" not implemented' %
                                  self.__class__.__name__)

    def __repr__(self):
        '''Return a string representation for this AST object'''
        return codegen.CodeGen().generate(self)

    def __str__(self):
        '''Return a string representation for this AST object'''
        return repr(self)
    
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
        return FunCallExp(self.exp.replicate(), [a.replicate() for a in self.args], self.line_no)

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

    def __init__(self, exp, op_type, line_no = ''):
        '''Create a unary operation expression'''
        Exp.__init__(self, line_no)
        self.exp = exp
        self.op_type = op_type

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return UnaryExp(self.exp.replicate(), self.op_type, self.line_no)

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
        self.exp = exp         # may be null

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        r_e = self.exp
        if r_e:
            r_e = r_e.replicate()
        return ExpStmt(r_e, self.line_no)

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
    
#-----------------------------------------------
# If-Then-Else
#-----------------------------------------------

class IfStmt(Stmt):

    def __init__(self, test, true_stmt, false_stmt = None, line_no = ''):
        '''Create an if statement'''
        Stmt.__init__(self, line_no)
        self.test = test
        self.true_stmt = true_stmt
        self.false_stmt = false_stmt           # may be null

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        f_s = self.false_stmt
        if f_s:
            f_s = f_s.replicate()
        return IfStmt(self.test.replicate(), self.true_stmt.replicate(), f_s, self.line_no)

#-----------------------------------------------
# For Loop
#-----------------------------------------------

class ForStmt(Stmt):

    def __init__(self, init, test, iter, stmt, line_no = ''):
        '''Create a for-loop statement'''
        Stmt.__init__(self, line_no)
        self.init = init      # may be null
        self.test = test      # may be null
        self.iter = iter      # may be null
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

#-----------------------------------------------
# Transformation
#-----------------------------------------------

class TransformStmt(Stmt):

    def __init__(self, name, args, stmt, line_no = ''):
        '''Create a transformation statement'''
        Stmt.__init__(self, line_no)
        self.name = name
        self.args = args
        self.stmt = stmt

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return TransformStmt(self.name, self.args[:], self.stmt.replicate(), self.line_no)

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

#-----------------------------------------------
# Pragma Directive
#-----------------------------------------------

class Pragma(NewAST):

    def __init__(self, pstring, line_no = ''):
        '''Create a pragma directive'''
        NewAST.__init__(self, line_no)
        self.pstring = pstring

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return Pragma(self.pstring, self.line_no)

#-----------------------------------------------
# Container
#-----------------------------------------------

class Container(NewAST):

    def __init__(self, ast, line_no = ''):
        '''Create a container AST (to protect the contained AST from any code transformations)'''
        NewAST.__init__(self, line_no)
        self.ast = ast

    def replicate(self):
        '''Replicate this abstract syntax tree node'''
        return Container(self.ast.replicate(), self.line_no)



