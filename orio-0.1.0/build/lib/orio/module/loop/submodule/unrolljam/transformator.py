#
# Contain the transformation procedure
#

import sys
import module.loop.ast, module.loop.ast_lib.constant_folder, module.loop.ast_lib.forloop_lib

#-----------------------------------------

class Transformator:
    '''Code transformator'''

    def __init__(self, ufactor, do_jamming, stmt):
        '''To instantiate a code transformator object'''

        self.ufactor = ufactor
        self.do_jamming = do_jamming
        self.stmt = stmt
        self.flib = module.loop.ast_lib.forloop_lib.ForLoopLib()
        self.cfolder = module.loop.ast_lib.constant_folder.ConstFolder()
        
    #----------------------------------------------------------

    def __addIdentWithExp(self, tnode, index_name, exp):
        '''Traverse the tree node and add any matching identifier with the provided expression'''

        if isinstance(exp, module.loop.ast.NumLitExp) and exp.val == 0:
            return tnode
        
        if isinstance(tnode, module.loop.ast.ExpStmt):
            if tnode.exp:
                tnode.exp = self.__addIdentWithExp(tnode.exp, index_name, exp)
            return tnode
    
        elif isinstance(tnode, module.loop.ast.CompStmt):
            tnode.stmts = [self.__addIdentWithExp(s, index_name, exp) for s in tnode.stmts]
            return tnode

        elif isinstance(tnode, module.loop.ast.IfStmt):
            tnode.test = self.__addIdentWithExp(tnode.test, index_name, exp)
            tnode.true_stmt = self.__addIdentWithExp(tnode.true_stmt, index_name, exp)
            if tnode.false_stmt:
                tnode.false_stmt = self.__addIdentWithExp(tnode.false_stmt, index_name, exp)
            return tnode

        elif isinstance(tnode, module.loop.ast.ForStmt):
            if tnode.init:
                tnode.init = self.__addIdentWithExp(tnode.init, index_name, exp)
            if tnode.test:
                tnode.test = self.__addIdentWithExp(tnode.test, index_name, exp)
            if tnode.iter:
                tnode.iter = self.__addIdentWithExp(tnode.iter, index_name, exp)
            tnode.stmt = self.__addIdentWithExp(tnode.stmt, index_name, exp)
            return tnode

        elif isinstance(tnode, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)

        elif isinstance(tnode, module.loop.ast.NumLitExp):
            return tnode

        elif isinstance(tnode, module.loop.ast.StringLitExp):
            return tnode

        elif isinstance(tnode, module.loop.ast.IdentExp):
            if tnode.name != index_name:
                return tnode
            else:
                add_exp = module.loop.ast.BinOpExp(tnode,
                                                   exp.replicate(),
                                                   module.loop.ast.BinOpExp.ADD)
                return module.loop.ast.ParenthExp(add_exp)

        elif isinstance(tnode, module.loop.ast.ArrayRefExp):
            tnode.exp = self.__addIdentWithExp(tnode.exp, index_name, exp)
            tnode.sub_exp = self.__addIdentWithExp(tnode.sub_exp, index_name, exp)
            return tnode
        
        elif isinstance(tnode, module.loop.ast.FunCallExp):
            tnode.exp = self.__addIdentWithExp(tnode.exp, index_name, exp)
            tnode.args = [self.__addIdentWithExp(a, index_name, exp) for a in tnode.args]
            return tnode

        elif isinstance(tnode, module.loop.ast.UnaryExp):
            tnode.exp = self.__addIdentWithExp(tnode.exp, index_name, exp)
            return tnode
        
        elif isinstance(tnode, module.loop.ast.BinOpExp):
            tnode.lhs = self.__addIdentWithExp(tnode.lhs, index_name, exp)
            tnode.rhs = self.__addIdentWithExp(tnode.rhs, index_name, exp)
            return tnode

        elif isinstance(tnode, module.loop.ast.ParenthExp):
            tnode.exp = self.__addIdentWithExp(tnode.exp, index_name, exp)
            return tnode        
        
        elif isinstance(tnode, module.loop.ast.NewAST):
            return tnode
        
        else:
            print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
            sys.exit(1)
    
    #-----------------------------------------

    def __jamStmts(self, stmtss):
        '''Jam/fuse statements whenever possible'''

        if len(stmtss) == 0:
            return module.loop.ast.CompStmt([])
        if len(stmtss) == 1:
            return module.loop.ast.CompStmt(stmtss[0])

        num = len(stmtss[0])
        for stmts in stmtss:
            assert(num == len(stmts)), 'internal error: unequal length of statement list'

        is_jam_valid = True
        contain_loop = False
        for i in range(num):
            s1 = None
            for stmts in stmtss:
                if s1 == None:
                    s1 = stmts[i]
                    if isinstance(s1, module.loop.ast.ForStmt):
                        contain_loop = True
                elif isinstance(s1, module.loop.ast.ForStmt):
                    s2 = stmts[i]
                    assert(isinstance(s2, module.loop.ast.ForStmt)), 'internal error: not a loop statement'
                    if not (str(s1.init) == str(s2.init) and str(s1.test) == str(s2.test) and str(s1.iter) == str(s2.iter)):
                        is_jam_valid = False
        if is_jam_valid:
            if not contain_loop:
                is_jam_valid = False

        if not is_jam_valid:
            n_stmts = []
            for stmts in stmtss:
                n_stmts.extend(stmts)
            return module.loop.ast.CompStmt(n_stmts)

        n_stmts = []
        for stmts in zip(*stmtss):
            if isinstance(stmts[0], module.loop.ast.ForStmt):
                l_stmtss = []
                for s in stmts:
                    if isinstance(s.stmt, module.loop.ast.CompStmt):
                        l_stmtss.append(s.stmt.stmts)
                    else:
                        l_stmtss.append([s.stmt])
                loop = stmts[0].replicate()
                loop.stmt = self.__jamStmts(l_stmtss)
                n_stmts.append(loop)
            else:
                n_stmts.extend(stmts)
        return module.loop.ast.CompStmt(n_stmts)
        
    #-----------------------------------------

    def transform(self):
        '''To unroll-and-jam the enclosed for-loop'''

        # get rid of compound statement that contains only a single statement
        while isinstance(self.stmt.stmt, module.loop.ast.CompStmt) and len(self.stmt.stmt.stmts) == 1:
            self.stmt.stmt = self.stmt.stmt.stmts[0]
        
        # extract for-loop structure
        for_loop_info = self.flib.extractForLoopInfo(self.stmt)
        index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
        
        # when ufactor = 1, no transformation will be applied
        if self.ufactor == 1:
            return self.flib.createForLoop(index_id, lbound_exp, ubound_exp,
                                           stride_exp, loop_body)
        
        # start generating the main unrolled loop
        # compute lower bound --> new_LB = LB
        new_lbound_exp = lbound_exp.replicate()
    
        # compute upper bound --> new_UB = UB-ST*(UF-1)
        it = module.loop.ast.BinOpExp(stride_exp.replicate(),
                                      module.loop.ast.NumLitExp(self.ufactor - 1,
                                                                module.loop.ast.NumLitExp.INT),
                                      module.loop.ast.BinOpExp.MUL)
        new_ubound_exp = module.loop.ast.BinOpExp(ubound_exp.replicate(),
                                                  it,
                                                  module.loop.ast.BinOpExp.SUB)
        new_ubound_exp = self.cfolder.fold(new_ubound_exp)
    
        # compute stride --> new_ST = UF*ST
        it = module.loop.ast.NumLitExp(self.ufactor, module.loop.ast.NumLitExp.INT)
        new_stride_exp = module.loop.ast.BinOpExp(it,
                                                  stride_exp.replicate(),
                                                  module.loop.ast.BinOpExp.MUL)
        new_stride_exp = self.cfolder.fold(new_stride_exp)
    
        # compute unrolled statements
        unrolled_stmt_seqs = []
        for i in range(0, self.ufactor):
            s = loop_body.replicate()
            it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
            increment_exp = module.loop.ast.BinOpExp(it,
                                                     stride_exp.replicate(),
                                                     module.loop.ast.BinOpExp.MUL)
            increment_exp = self.cfolder.fold(increment_exp)
            ns = self.__addIdentWithExp(s, index_id.name, increment_exp)
            ns = self.cfolder.fold(ns)
            if isinstance(ns, module.loop.ast.CompStmt):
                unrolled_stmt_seqs.append(ns.stmts)
            else:
                unrolled_stmt_seqs.append([ns])

        # compute the unrolled loop body by jamming/fusing the unrolled statements
        if self.do_jamming:
            unrolled_loop_body = self.__jamStmts(unrolled_stmt_seqs)
        else:
            unrolled_stmts = reduce(lambda x,y: x+y, unrolled_stmt_seqs, [])
            unrolled_loop_body = module.loop.ast.CompStmt(unrolled_stmts)
            
        # generate the main unrolled loop
        main_loop = self.flib.createForLoop(index_id, new_lbound_exp, new_ubound_exp,
                                            new_stride_exp, unrolled_loop_body)
        
        # generate the clean-up loop
        cleanup_loop = self.flib.createForLoop(index_id, None, ubound_exp,
                                               stride_exp, loop_body)
        
        # generate the transformed statement
        transformed_stmt = module.loop.ast.CompStmt([main_loop, cleanup_loop])

        # return the transformed statement
        return transformed_stmt


