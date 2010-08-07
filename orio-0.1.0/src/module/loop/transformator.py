#
# The transformator that applies code transformation procedures
#

import sys
import ast, main.dyn_loader, module.loop.codegen

#-----------------------------------------

# the name of the transformation submodule
TSUBMOD_NAME = 'module.loop.submodule'

#-----------------------------------------

class Transformator:
    '''Code transformator'''

    def __init__(self, perf_params, verbose):
        '''To instantiate a code transformator object'''

        self.perf_params = perf_params
        self.verbose = verbose
        self.dloader = main.dyn_loader.DynLoader()
        
    #--------------------------------------

    def __transformStmt(self, stmt):
        '''Apply code transformation on the given statement'''
 
        if isinstance(stmt, ast.ExpStmt):
            return stmt

        elif isinstance(stmt, ast.CompStmt):
            stmt.stmts = [self.__transformStmt(s) for s in stmt.stmts]
            return stmt

        elif isinstance(stmt, ast.IfStmt):
            stmt.true_stmt = self.__transformStmt(stmt.true_stmt)
            if stmt.false_stmt:
                stmt.false_stmt = self.__transformStmt(stmt.false_stmt)
            return stmt

        elif isinstance(stmt, ast.ForStmt):
            stmt.stmt = self.__transformStmt(stmt.stmt)
            return stmt

        elif isinstance(stmt, ast.TransformStmt):

            # transform the nested statement
            stmt.stmt = self.__transformStmt(stmt.stmt)

            # check for repeated transformation argument names
            arg_names = {}
            for [aname, rhs, line_no] in stmt.args:
                if aname in arg_names:
                    print 'error:%s: repeated transformation argument: "%s"' % (line_no, aname)
                    sys.exit(1)
                arg_names[aname] = None

            # dynamically load the transformation submodule class
            class_name = stmt.name
            submod_name = '.'.join([TSUBMOD_NAME, class_name.lower(), class_name.lower()])
            submod_class = self.dloader.loadClass(submod_name, class_name)
            
            # apply code transformations
            try:
                t = submod_class(self.perf_params, stmt.args, stmt.stmt)
                transformed_stmt = t.transform()
            except Exception, e:
                print (('error:%s: encountered an error as optimizing the transformation ' +
                        'statement: "%s"') % (stmt.line_no, class_name))
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)

            # return the transformed statement
            return transformed_stmt

        else:
            print 'internal error: unknown statement type: %s' % stmt.__class__.__name__
            sys.exit(1)
   
    #--------------------------------------

    def transform(self, stmts):
        '''Apply code transformations on each statement in the given statement list'''

        return [self.__transformStmt(s) for s in stmts]

        

