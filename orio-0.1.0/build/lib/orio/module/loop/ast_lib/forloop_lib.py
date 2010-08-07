#
# A library for for-loop statements
#

import sets, sys
import module.loop.ast

#-----------------------------------------

class ForLoopLib:
    '''A library tool used to provide a set of subroutines to process for-loop statements'''

    def __init__(self):
        '''To instantiate a for-loop library tool object'''
        pass

    #-------------------------------------------------

    def extractForLoopInfo(self, stmt):
        '''
        Given a for-loop statement, extract information about its loop structure
        Note that the for-loop must be in the following form:
          for (<id> = <exp>; <id> <= <exp>; <id> += <exp>)
            <stmt>
        Subtraction is not considered at the iteration expression for the sake of
        the implementation simplicity.
        '''

        # get rid of compound statement that contains only a single statement
        while isinstance(stmt, module.loop.ast.CompStmt) and len(stmt.stmts) == 1:
            stmt = stmt.stmts[0]

        # check if it is a for-loop statement
        if not isinstance(stmt, module.loop.ast.ForStmt):
            print 'error:%s: not a for-loop statement' % stmt.line_no
            sys.exit(1)

        # check initialization expression
        if stmt.init:
            while True:
                while isinstance(stmt.init, module.loop.ast.ParenthExp):
                    stmt.init = stmt.init.exp
                if (isinstance(stmt.init, module.loop.ast.BinOpExp) and
                    stmt.init.op_type == module.loop.ast.BinOpExp.EQ_ASGN):
                    while isinstance(stmt.init.lhs, module.loop.ast.ParenthExp):
                        stmt.init.lhs = stmt.init.lhs.exp
                    while isinstance(stmt.init.rhs, module.loop.ast.ParenthExp):
                        stmt.init.rhs = stmt.init.rhs.exp
                    if isinstance(stmt.init.lhs, module.loop.ast.IdentExp): 
                        break
                print ('error:%s: loop initialization expression not in "<id> = <exp>" form' %
                       stmt.init.line_no)
                sys.exit(1)
                
        # check test expression
        if stmt.test:
            while True:
                while isinstance(stmt.test, module.loop.ast.ParenthExp):
                    stmt.test = stmt.test.exp
                if (isinstance(stmt.test, module.loop.ast.BinOpExp) and
                    stmt.test.op_type == module.loop.ast.BinOpExp.LE):
                    while isinstance(stmt.test.lhs, module.loop.ast.ParenthExp):
                        stmt.test.lhs = stmt.test.lhs.exp
                    while isinstance(stmt.test.rhs, module.loop.ast.ParenthExp):
                        stmt.test.rhs = stmt.test.rhs.exp
                    if isinstance(stmt.test.lhs, module.loop.ast.IdentExp): 
                        break
                print ('error:%s: loop test expression not in "<id> <= <exp>" form' %
                       stmt.test.line_no)
                sys.exit(1)
            
        # check iteration expression
        if stmt.iter:
            while True:
                while isinstance(stmt.iter, module.loop.ast.ParenthExp):
                    stmt.iter = stmt.iter.exp
                if (isinstance(stmt.iter, module.loop.ast.BinOpExp) and
                    stmt.iter.op_type == module.loop.ast.BinOpExp.EQ_ASGN):
                    while isinstance(stmt.iter.lhs, module.loop.ast.ParenthExp):
                        stmt.iter.lhs = stmt.iter.lhs.exp
                    while isinstance(stmt.iter.rhs, module.loop.ast.ParenthExp):
                        stmt.iter.rhs = stmt.iter.rhs.exp
                    if isinstance(stmt.iter.lhs, module.loop.ast.IdentExp):
                        if (isinstance(stmt.iter.rhs, module.loop.ast.BinOpExp) and
                            stmt.iter.rhs.op_type in (module.loop.ast.BinOpExp.ADD,
                                                      module.loop.ast.BinOpExp.SUB)):
                            while isinstance(stmt.iter.rhs.lhs, module.loop.ast.ParenthExp):
                                stmt.iter.rhs.lhs = stmt.iter.rhs.lhs.exp
                            while isinstance(stmt.iter.rhs.rhs, module.loop.ast.ParenthExp):
                                stmt.iter.rhs.rhs = stmt.iter.rhs.rhs.exp
                            if (isinstance(stmt.iter.rhs.lhs, module.loop.ast.IdentExp) and
                                stmt.iter.lhs.name == stmt.iter.rhs.lhs.name):
                                break
                elif (isinstance(stmt.iter, module.loop.ast.UnaryExp) and
                      stmt.iter.op_type in (module.loop.ast.UnaryExp.POST_INC,
                                            module.loop.ast.UnaryExp.PRE_INC,
                                            module.loop.ast.UnaryExp.POST_DEC,
                                            module.loop.ast.UnaryExp.PRE_DEC)):
                    while isinstance(stmt.iter.exp, module.loop.ast.ParenthExp):
                        stmt.iter.exp = stmt.iter.exp.exp
                    if isinstance(stmt.iter.exp, module.loop.ast.IdentExp):
                        break
                print (('error:%s: loop iteration expression not in "<id>++" or "<id>--" or ' +
                        '"<id> += <exp>" or "<id> = <id> + <exp>" form') % stmt.iter.line_no)
                sys.exit(1)

        # check if the control expressions are all empty
        if not stmt.init and not stmt.test and not stmt.iter:
            print ('error:%s: a loop with an empty control expression cannot be handled' %
                   stmt.line_no)
            sys.exit(1)
    
        # check if the iterator names are all the same
        init_iname = None
        test_iname = None
        iter_iname = None
        if stmt.init:
            init_iname = stmt.init.lhs.name
        if stmt.test:
            test_iname = stmt.test.lhs.name
        if stmt.iter:
            if isinstance(stmt.iter, module.loop.ast.BinOpExp):
                iter_iname = stmt.iter.lhs.name
            else:
                assert(isinstance(stmt.iter, module.loop.ast.UnaryExp)), 'internal error: not unary'
                iter_iname = stmt.iter.exp.name
        inames = []
        if init_iname:
            inames.append(init_iname)
        if test_iname:
            inames.append(test_iname)
        if iter_iname:
            inames.append(iter_iname)
        if inames.count(inames[0]) != len(inames):
            print ('error:%s: iterator names across init, test, and iter exps must be the same'
                   % stmt.line_no)
            sys.exit(1)
        
        # extract for-loop structure information
        index_id = module.loop.ast.IdentExp(inames[0])
        lbound_exp = None
        ubound_exp = None
        stride_exp = None
        if stmt.init:
            lbound_exp = stmt.init.rhs.replicate()
        if stmt.test:
            ubound_exp = stmt.test.rhs.replicate()
        if stmt.iter:
            if isinstance(stmt.iter, module.loop.ast.BinOpExp):
                stride_exp = stmt.iter.rhs.rhs.replicate()
                if isinstance(stride_exp, module.loop.ast.BinOpExp):
                    stride_exp = module.loop.ast.ParenthExp(stride_exp)
                if stmt.iter.rhs.op_type == module.loop.ast.BinOpExp.SUB:
                    stride_exp = module.loop.ast.UnaryExp(stride_exp, module.loop.ast.UnaryExp.MINUS)
            elif isinstance(stmt.iter, module.loop.ast.UnaryExp):
                if stmt.iter.op_type in (module.loop.ast.UnaryExp.POST_INC,
                                         module.loop.ast.UnaryExp.PRE_INC):
                    stride_exp = module.loop.ast.NumLitExp(1, module.loop.ast.NumLitExp.INT)
                elif stmt.iter.op_type in (module.loop.ast.UnaryExp.POST_DEC,
                                           module.loop.ast.UnaryExp.PRE_DEC):
                    stride_exp = module.loop.ast.NumLitExp(-1, module.loop.ast.NumLitExp.INT)
                else:
                    print 'internal error: unexpected unary operation type'
                    sys.exit(1)
            else:
                print 'internal error: unexpected type of iteration expression'
                sys.exit(1)
        loop_body = stmt.stmt.replicate()
        for_loop_info = (index_id, lbound_exp, ubound_exp, stride_exp, loop_body)
        
        # return the for-loop structure information
        return for_loop_info

    #-------------------------------------------------
    
    def createForLoop(self, index_id, lbound_exp, ubound_exp, stride_exp, loop_body):
        '''
        Generate a for loop:
          for (index_id = lbound_exp; index_id <= ubound_exp; index_id = index_id + stride_exp)
            loop_body
        '''

        init_exp = None
        test_exp = None
        iter_exp = None
        if lbound_exp:
            init_exp = module.loop.ast.BinOpExp(index_id.replicate(),
                                                lbound_exp.replicate(),
                                                module.loop.ast.BinOpExp.EQ_ASGN)
        if ubound_exp:
            test_exp = module.loop.ast.BinOpExp(index_id.replicate(),
                                                ubound_exp.replicate(),
                                                module.loop.ast.BinOpExp.LE)
        if stride_exp:
            while isinstance(stride_exp, module.loop.ast.ParenthExp):
                stride_exp = stride_exp.exp
            it = module.loop.ast.BinOpExp(index_id.replicate(),
                                          stride_exp.replicate(),
                                          module.loop.ast.BinOpExp.ADD)
            iter_exp = module.loop.ast.BinOpExp(index_id.replicate(),
                                                it,
                                                module.loop.ast.BinOpExp.EQ_ASGN)
        return module.loop.ast.ForStmt(init_exp, test_exp, iter_exp, loop_body.replicate())
    
    #-------------------------------------------------

    def getLoopIndexNames(self, stmt):
        '''Return a list of all loop index names'''

        if isinstance(stmt, module.loop.ast.ExpStmt):
            return []

        elif isinstance(stmt, module.loop.ast.CompStmt):
            inames = []
            for s in stmt.stmts:
                inames.extend(self.getLoopIndexNames(s))
            return list(sets.Set(inames))

        elif isinstance(stmt, module.loop.ast.IfStmt):
            inames = []
            inames.extend(self.getLoopIndexNames(stmt.true_stmt))
            if stmt.false_stmt:
                inames.extend(self.getLoopIndexNames(stmt.false_stmt))
            return list(sets.Set(inames))

        elif isinstance(stmt, module.loop.ast.ForStmt):
            inames = []
            inames.extend(self.getLoopIndexNames(stmt.stmt))
            index_id, lbound_exp, ubound_exp, stride_exp, loop_body = self.extractForLoopInfo(stmt)
            if index_id.name not in inames:
                inames.append(index_id.name)
            return inames

        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)
                        
        elif isinstance(stmt, module.loop.ast.NewAST):
            return []

        else:
            print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
            sys.exit(1)

    #-------------------------------------------------

    def hasInnerLoop(self, stmt):
        '''Determine if there is an inner loop inside the given statement'''

        if stmt == None:
            return False
        
        if isinstance(stmt, module.loop.ast.ExpStmt):
            return False

        elif isinstance(stmt, module.loop.ast.CompStmt):
            for s in stmt.stmts:
                if self.hasInnerLoop(s):
                    return True
            return False

        elif isinstance(stmt, module.loop.ast.IfStmt):
            if self.hasInnerLoop(stmt.true_stmt):
                return True
            else:
                return self.hasInnerLoop(stmt.false_stmt)

        elif isinstance(stmt, module.loop.ast.ForStmt):
            return True

        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)
                        
        elif isinstance(stmt, module.loop.ast.NewAST):
            return False

        else:
            print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
            sys.exit(1)
