#
# Contain the transformation procedure
#

import sys
import module.loop.ast, module.loop.ast_lib.common_lib, module.loop.ast_lib.forloop_lib

#-----------------------------------------

class Transformator:
    '''Code transformator'''

    def __init__(self, lprefix, uprefix, stmt):
        '''To instantiate a code transformator object'''

        self.lprefix = 'lb_'
        self.uprefix = 'ub_'
        if lprefix != None:
            self.lprefix = lprefix
        if uprefix != None:            
            self.uprefix = uprefix
        self.stmt = stmt
        
        self.dtype = 'register int'
        self.lcounter = 1
        self.ucounter = 1
        self.flib = module.loop.ast_lib.forloop_lib.ForLoopLib()
        self.clib = module.loop.ast_lib.common_lib.CommonLib()
        
    #----------------------------------------------------------

    def __replaceBounds(self, stmt, iter_name):
        '''
        To replace loop bounds with scalar intermediates and to hoist the loop bounds
        initializations out of the affecting loops
        '''
        
        if isinstance(stmt, module.loop.ast.ExpStmt):
            return (stmt, [])

        elif isinstance(stmt, module.loop.ast.CompStmt):
            stmts = []
            asgns = []
            for s in stmt.stmts:
                is_comp_before = isinstance(s, module.loop.ast.CompStmt)
                (ns, nasgns) = self.__replaceBounds(s, iter_name)
                is_comp_after = isinstance(ns, module.loop.ast.CompStmt)
                if not is_comp_before and is_comp_after:
                    stmts.extend(ns.stmts)
                else:
                    stmts.append(ns)
                asgns.extend(nasgns)
            stmt.stmts = stmts
            return (stmt, asgns)

        elif isinstance(stmt, module.loop.ast.IfStmt):
            asgns = []
            nstmt, nasgns = self.__replaceBounds(stmt.true_stmt, iter_name)
            stmt.true_stmt = nstmt
            asgns.extend(nasgns)
            if stmt.false_stmt:
                nstmt, nasgns = self.__replaceBounds(stmt.false_stmt, iter_name)
                stmt.false_stmt = nstmt
                asgns.extend(nasgns)
            return (stmt, asgns)
                                        
        elif isinstance(stmt, module.loop.ast.ForStmt):

            # extract the for-loop structure
            for_loop_info = self.flib.extractForLoopInfo(stmt)
            index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info 

            # get the iteration variable name
            niter_name = index_id.name

            # recursion on the loop body
            nstmt, nasgns = self.__replaceBounds(stmt.stmt, niter_name)
            stmt.stmt = nstmt

            # replace loop bounds if they are complex
            asgns = nasgns
            if lbound_exp and self.clib.isComplexExp(lbound_exp):
                intmd = module.loop.ast.IdentExp(self.lprefix + str(self.lcounter)) 
                self.lcounter += 1
                stmt.init.rhs = intmd.replicate()
                asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                                lbound_exp.replicate(),
                                                module.loop.ast.BinOpExp.EQ_ASGN) 
                asgn = module.loop.ast.ExpStmt(asgn)
                asgns.append(asgn)
            if ubound_exp and self.clib.isComplexExp(ubound_exp):
                intmd = module.loop.ast.IdentExp(self.uprefix + str(self.ucounter)) 
                self.ucounter += 1
                stmt.test.rhs = intmd.replicate()
                asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                                ubound_exp.replicate(),
                                                module.loop.ast.BinOpExp.EQ_ASGN) 
                asgn = module.loop.ast.ExpStmt(asgn)
                asgns.append(asgn)

            # generate the transformed loop
            cdecls = []
            casgns = []
            nasgns = []
            for asgn in asgns:
                bound_exp = asgn.exp.rhs
                if iter_name != None and self.clib.containIdentName(bound_exp, iter_name):
                    intmd_name = asgn.exp.lhs.name
                    if len(cdecls) == 0 or len(cdecls[-1].var_names) > 8: 
                        cdecls.append(module.loop.ast.VarDecl(self.dtype, [intmd_name]))
                    else:
                        cdecls[-1].var_names.append(intmd_name)
                        cdecls[-1].var_names.sort()
                    casgns.append(asgn)
                else:
                    nasgns.append(asgn)
            if len(casgns) > 0:
                nstmt = module.loop.ast.CompStmt(cdecls + casgns + [stmt])
            else:
                nstmt = stmt

            # return the transformed loop
            return (nstmt, nasgns)
            
        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)

        elif isinstance(stmt, module.loop.ast.NewAST):
            return (stmt, [])

        else:
            print 'internal error: unexpected AST type: "%s"' % exp.__class__.__name__
            sys.exit(1)
                                    
    #----------------------------------------------------------

    def transform(self):
        '''To replace loop bounds with scalars, and also perform hoisting on the loop bounds'''

        # reset counters
        self.lcounter = 1
        self.ucounter = 1

        # perform loop bounds replacement
        transformed_stmt, asgns = self.__replaceBounds(self.stmt, None)

        # create the loop bounds declarations
        decls = []
        for asgn in asgns:
            intmd_name = asgn.exp.lhs.name
            if len(decls) == 0 or len(decls[-1].var_names) > 8:
                decls.append(module.loop.ast.VarDecl(self.dtype, [intmd_name]))
            else:
                decls[-1].var_names.append(intmd_name)
                decls[-1].var_names.sort()

        if len(asgns) > 0:
            if isinstance(transformed_stmt, module.loop.ast.CompStmt):
                transformed_stmt.stmts = decls + asgns + transformed_stmt.stmts
            else:
                transformed_stmt = module.loop.ast.CompStmt(decls + asgns + [transformed_stmt])

        # return the transformed statement
        return transformed_stmt
