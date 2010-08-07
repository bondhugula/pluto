#
# Contain a class that provides a set of common library functions for AST processing
#

import sys
import module.loop.ast

#-----------------------------------------------------------
 
class CommonLib:
    '''A common library set for AST processing'''
    
    def __init__(self):
        '''To instantiate a common library object'''
        pass

    #-------------------------------------------------------

    def replaceIdent(self, tnode, iname_from, iname_to):
        '''Replace the names of all matching identifiers with the given name'''

        if isinstance(tnode, module.loop.ast.NumLitExp):
            return tnode
            
        elif isinstance(tnode, module.loop.ast.StringLitExp):
            return tnode
            
        elif isinstance(tnode, module.loop.ast.IdentExp):
            if tnode.name == iname_from:
                tnode.name = iname_to
            return tnode
            
        elif isinstance(tnode, module.loop.ast.ArrayRefExp):
            tnode.exp = self.replaceIdent(tnode.exp, iname_from, iname_to)
            tnode.sub_exp = self.replaceIdent(tnode.sub_exp, iname_from, iname_to)
            return tnode
            
        elif isinstance(tnode, module.loop.ast.FunCallExp):
            tnode.exp = self.replaceIdent(tnode.exp, iname_from, iname_to)
            tnode.args = [self.replaceIdent(a, iname_from, iname_to) for a in tnode.args]
            return tnode
            
        elif isinstance(tnode, module.loop.ast.UnaryExp):
            tnode.exp = self.replaceIdent(tnode.exp, iname_from, iname_to)
            return tnode
            
        elif isinstance(tnode, module.loop.ast.BinOpExp):
            tnode.lhs = self.replaceIdent(tnode.lhs, iname_from, iname_to)
            tnode.rhs = self.replaceIdent(tnode.rhs, iname_from, iname_to)
            return tnode
            
        elif isinstance(tnode, module.loop.ast.ParenthExp):
            tnode.exp = self.replaceIdent(tnode.exp, iname_from, iname_to)
            return tnode
            
        elif isinstance(tnode, module.loop.ast.ExpStmt):
            if tnode.exp:
                tnode.exp = self.replaceIdent(tnode.exp, iname_from, iname_to)
            return tnode
            
        elif isinstance(tnode, module.loop.ast.CompStmt):
            tnode.stmts = [self.replaceIdent(s, iname_from, iname_to) for s in tnode.stmts]
            return tnode
            
        elif isinstance(tnode, module.loop.ast.IfStmt):
            tnode.test = self.replaceIdent(tnode.test, iname_from, iname_to)
            tnode.true_stmt = self.replaceIdent(tnode.true_stmt, iname_from, iname_to)
            if tnode.false_stmt:
                tnode.false_stmt = self.replaceIdent(tnode.false_stmt, iname_from, iname_to)
            return tnode
            
        elif isinstance(tnode, module.loop.ast.ForStmt):
            if tnode.init:
                tnode.init = self.replaceIdent(tnode.init, iname_from, iname_to)
            if tnode.test:
                tnode.test = self.replaceIdent(tnode.test, iname_from, iname_to)
            if tnode.iter:
                tnode.iter = self.replaceIdent(tnode.iter, iname_from, iname_to)
            tnode.stmt = self.replaceIdent(tnode.stmt, iname_from, iname_to)
            return tnode

        elif isinstance(tnode, module.loop.ast.TransformStmt):
            print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
            sys.exit(1)
        
        elif isinstance(tnode, module.loop.ast.NewAST):
            return tnode
        
        else:
            print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
            sys.exit(1)
        
    #-------------------------------------------------------

    def containIdentName(self, exp, iname):
        '''
        Check if the given expression contains an identifier whose name matches to the given name
        '''

        if exp == None:
            return False
        
        if isinstance(exp, module.loop.ast.NumLitExp):
            return False
        
        elif isinstance(exp, module.loop.ast.StringLitExp):
            return False
        
        elif isinstance(exp, module.loop.ast.IdentExp):
            return exp.name == iname
        
        elif isinstance(exp, module.loop.ast.ArrayRefExp):
            return self.containIdentName(exp.exp, iname) or self.containIdentName(exp.sub_exp, iname)
        
        elif isinstance(exp, module.loop.ast.FunCallExp):
            has_match = reduce(lambda x,y: x or y,
                               [self.containIdentName(a, iname) for a in exp.args],
                               False)
            return self.containIdentName(exp.exp, iname) or has_match
        
        elif isinstance(exp, module.loop.ast.UnaryExp):
            return self.containIdentName(exp.exp, iname)
        
        elif isinstance(exp, module.loop.ast.BinOpExp):
            return self.containIdentName(exp.lhs, iname) or self.containIdentName(exp.rhs, iname)
        
        elif isinstance(exp, module.loop.ast.ParenthExp):
            return self.containIdentName(exp.exp, iname)
        
        elif isinstance(exp, module.loop.ast.NewAST):
            return False
        
        else:
            print 'internal error: unexpected AST type: "%s"' % exp.__class__.__name__
            sys.exit(1)
            
    #-------------------------------------------------------

    def isComplexExp(self, exp):
        '''
        To determine if the given expression is complex. Simple expressions contain only a variable
        or a number or a string.
        '''
        
        if isinstance(exp, module.loop.ast.NumLitExp):
            return False
        
        # a rare case
        elif isinstance(exp, module.loop.ast.StringLitExp):
            return False
        
        elif isinstance(exp, module.loop.ast.IdentExp):
            return False
        
        # a rare case
        elif isinstance(exp, module.loop.ast.ArrayRefExp):
            return True
        
        elif isinstance(exp, module.loop.ast.FunCallExp):
            return True
        
        elif isinstance(exp, module.loop.ast.UnaryExp):
            return self.isComplexExp(exp.exp)
        
        elif isinstance(exp, module.loop.ast.BinOpExp):
            return True
        
        elif isinstance(exp, module.loop.ast.ParenthExp):
            return self.isComplexExp(exp.exp)
        
        # a rare case
        elif isinstance(exp, module.loop.ast.NewAST):
            return True
        
        else:
            print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
            sys.exit(1)
            
            
            
            
