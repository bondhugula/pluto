#
# Contain the transformation procedure
#

import sys
import module.loop.ast, module.loop.ast_lib.constant_folder, module.loop.ast_lib.forloop_lib
import module.loop.ast_lib.common_lib, semant

#-----------------------------------------

class Transformator:
    '''Code transformator'''

    __ivar_suffix = 't'
    __newlb_prefix = 'newlb_'
    __newub_prefix = 'newub_'
    __int_min = '-2147483648'
    __int_max = '2147483647'

    #---------------------------------------------------------

    def __init__(self, loops, ufactors, stmt):
        '''To instantiate a code transformator object'''

        self.loops = loops
        self.ufactors = ufactors
        self.stmt = stmt

        self.ufactor_map = dict(zip(self.loops, self.ufactors))
        self.itvar_map = dict(zip(self.loops, [l+self.__ivar_suffix for l in self.loops]))
        self.itvar_map_rev = dict(zip([l+self.__ivar_suffix for l in self.loops], self.loops))

        self.flib = module.loop.ast_lib.forloop_lib.ForLoopLib()
        self.cfolder = module.loop.ast_lib.constant_folder.ConstFolder()
        self.clib = module.loop.ast_lib.common_lib.CommonLib()

        # to always perform statement jamming
        self.use_jam = True
        
        # to use full unrolling (or just use an intra-tile loop) for handling unjammable statements
        self.use_full_unroll = False

        # to use non-rectangular register tiling
        self.use_nonrect_tile = True

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
            if tnode.name == index_name:
                add_exp = module.loop.ast.BinOpExp(tnode,
                                                   exp.replicate(),
                                                   module.loop.ast.BinOpExp.ADD)
                return module.loop.ast.ParenthExp(add_exp)
            else:
                return tnode

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

    #----------------------------------------------------------

    def __createIntraTileLoop(self, intratile_iname, intertile_iname,
                              ufactor, stride_exp, loop_body):
        '''Create the intra tile loop'''

        # create the iteration variable
        index_id = module.loop.ast.IdentExp(intratile_iname)

        # compute lower bound --> LB' = IT
        lbound_exp = module.loop.ast.IdentExp(intertile_iname)

        # compute upper bound --> UB' = IT+ST*(UF-1)
        it = module.loop.ast.BinOpExp(stride_exp.replicate(),
                                      module.loop.ast.NumLitExp(ufactor - 1,
                                                                module.loop.ast.NumLitExp.INT),
                                      module.loop.ast.BinOpExp.MUL)
        ubound_exp = module.loop.ast.BinOpExp(module.loop.ast.IdentExp(intertile_iname),
                                              it,
                                              module.loop.ast.BinOpExp.ADD)
        ubound_exp = self.cfolder.fold(ubound_exp)

        # compute stride --> ST' = ST
        stride_exp = stride_exp.replicate()

        # generate the intra-tile loop
        intratile_loop = self.flib.createForLoop(index_id, lbound_exp, ubound_exp,
                                                 stride_exp, loop_body.replicate())

        # return the intra-tile loop
        return intratile_loop

    #----------------------------------------------------------

    def __createUnrolledLoops(self, for_loop_info, ufactor):
        '''Create the main unrolled loop and the cleanup loop'''

        # unpack the original for-loop information
        index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info

        # create loop indexes for both the main and cleanup loops
        mloop_index_id = module.loop.ast.IdentExp(self.itvar_map[index_id.name])
        cloop_index_id = index_id

        # check for-loop structure
        if lbound_exp == None or ubound_exp == None or stride_exp == None:
            print ('error: the loop to be unrolled/jammed must NOT have empty bounds/stride ' +
                   'expressions')
            sys.exit(1)

        # compute lower bound --> new_LB = LB
        new_lbound_exp = lbound_exp.replicate()

        # compute upper bound --> new_UB = UB-ST*(UF-1)
        it = module.loop.ast.BinOpExp(stride_exp.replicate(),
                                      module.loop.ast.NumLitExp(ufactor - 1,
                                                                module.loop.ast.NumLitExp.INT),
                                      module.loop.ast.BinOpExp.MUL)
        new_ubound_exp = module.loop.ast.BinOpExp(ubound_exp.replicate(),
                                                  it,
                                                  module.loop.ast.BinOpExp.SUB)
        new_ubound_exp = self.cfolder.fold(new_ubound_exp)
        
        # compute stride --> new_ST = UF*ST
        it = module.loop.ast.NumLitExp(ufactor, module.loop.ast.NumLitExp.INT)
        new_stride_exp = module.loop.ast.BinOpExp(it,
                                                  stride_exp.replicate(),
                                                  module.loop.ast.BinOpExp.MUL)
        new_stride_exp = self.cfolder.fold(new_stride_exp)

        # compute the loop body
        mloop_body = self.clib.replaceIdent(loop_body.replicate(),
                                            cloop_index_id.name, mloop_index_id.name)

        # generate the main unrolled loop
        main_loop = self.flib.createForLoop(mloop_index_id, new_lbound_exp, new_ubound_exp,
                                            new_stride_exp, mloop_body)

        # generate the cleanup loop
        cleanup_loop = self.flib.createForLoop(cloop_index_id, mloop_index_id, ubound_exp,
                                               stride_exp, loop_body)
        
        # return the information about the tiled loops
        return (main_loop, cleanup_loop)

    #----------------------------------------------------------

    def __transformLoop(self, stmt, outer_unrolls):
        '''To traverse the abstract syntax tree to unroll/jam each encountered for-loop'''
        
        if isinstance(stmt, module.loop.ast.ExpStmt):

            # unroll this statement based on the given sequence of unroll factors
            stmts = [stmt]
            outer_unrolls_rev = outer_unrolls[:]
            outer_unrolls_rev.reverse()
            for iname, ufactor, stride_exp in outer_unrolls_rev:

                # skip unroll factor of 1
                if ufactor == 1:
                    continue

                # unroll this statement
                n_stmts = []
                for i in range(0, ufactor):
                    it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
                    increment_exp = module.loop.ast.BinOpExp(it,
                                                             stride_exp.replicate(),
                                                             module.loop.ast.BinOpExp.MUL)
                    increment_exp = self.cfolder.fold(increment_exp)
                    n_stmts += [self.__addIdentWithExp(s.replicate(), iname, increment_exp) 
                                for s in stmts]
                stmts = n_stmts

            # return the unrolled statement
            if len(stmts) == 1:
                return stmts[0]
            return module.loop.ast.CompStmt(stmts)
        
        elif isinstance(stmt, module.loop.ast.CompStmt):

            # unroll/jam this compound statement
            stmts = []
            for s in stmt.stmts:

                # unroll/jam the inner statement
                is_comp_before = isinstance(s, module.loop.ast.CompStmt)
                ns = self.__transformLoop(s, outer_unrolls)
                is_comp_after = isinstance(ns, module.loop.ast.CompStmt)

                # combine the unrolled inner statement
                if not is_comp_before and is_comp_after:
                    stmts += ns.stmts
                else:
                    stmts.append(ns)
            stmt.stmts = stmts

            # return the unrolled/jammed statement
            return stmt

        elif isinstance(stmt, module.loop.ast.IfStmt):

            # separate the bound and free unroll factors
            bound_unrolls = []
            free_unrolls = []
            for u in outer_unrolls:
                iname, ufactor, stride_exp = u
                if ufactor == 1:
                    continue
                if self.clib.containIdentName(stmt.test, iname):
                    bound_unrolls.append(u)
                else:
                    free_unrolls.append(u)

            # unroll this statement using the free unroll factors
            stmt.true_stmt = self.__transformLoop(stmt.true_stmt, free_unrolls)
            if stmt.false_stmt:
                stmt.false_stmt = self.__transformLoop(stmt.false_stmt, free_unrolls)

            # unroll this statement using the bound unroll factors (i.e. unjammable statement)
            stmts = [stmt]
            bound_unrolls.reverse()
            for iname, ufactor, stride_exp in bound_unrolls:

                # skip unroll factor of 1
                if ufactor == 1:
                    continue

                # to fully unroll the unjammable statement
                if self.use_full_unroll:
                    n_stmts = []
                    for i in range(0, ufactor):
                        it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
                        increment_exp = module.loop.ast.BinOpExp(it,
                                                                 stride_exp.replicate(),
                                                                 module.loop.ast.BinOpExp.MUL)
                        increment_exp = self.cfolder.fold(increment_exp)
                        n_stmts += [self.__addIdentWithExp(s.replicate(), iname, increment_exp) 
                                    for s in stmts]
                    stmts = n_stmts

                # to create an intra-tile loop (not to fully unroll the unjammable statement)
                else:
                    intratile_iname = self.itvar_map_rev[iname]
                    loop_body = self.clib.replaceIdent(stmts[0], iname, intratile_iname)
                    nstmt = self.__createIntraTileLoop(intratile_iname, iname, ufactor,
                                                       stride_exp, loop_body)
                    stmts = [nstmt]

            # return the unrolled statement
            if len(stmts) == 1:
                return stmts[0]
            return module.loop.ast.CompStmt(stmts)

        elif isinstance(stmt, module.loop.ast.ForStmt):

            # extract for-loop structure
            i_for_loop_info = self.flib.extractForLoopInfo(stmt)
            i_index_id, i_lbound_exp, i_ubound_exp, i_stride_exp, i_loop_body = i_for_loop_info

            # get the unroll factor of this loop
            i_ufactor = 1
            if i_index_id.name in self.ufactor_map:
                i_ufactor = self.ufactor_map[i_index_id.name]

            # separate the bound and free unroll factors
            bound_unrolls = []
            free_unrolls = []
            n_outer_unrolls = []
            for u in outer_unrolls:
                iname, ufactor, stride_exp = u
                if ufactor == 1:
                    continue
                if (self.clib.containIdentName(i_lbound_exp, iname) or
                    self.clib.containIdentName(i_ubound_exp, iname) or
                    self.clib.containIdentName(i_stride_exp, iname)):
                    bound_unrolls.append(u)
                else:
                    free_unrolls.append(u)
                n_outer_unrolls.append(u)
            outer_unrolls = n_outer_unrolls
            
            # get the original and extended free/outer unroll factors
            orig_free_unrolls = free_unrolls[:]
            ext_free_unrolls = free_unrolls
            orig_outer_unrolls = outer_unrolls[:]
            ext_outer_unrolls = outer_unrolls            
            if i_ufactor != 1:
                intertile_iname = self.itvar_map[i_index_id.name]
                ext_free_unrolls.append((intertile_iname, i_ufactor, i_stride_exp))
                ext_outer_unrolls.append((intertile_iname, i_ufactor, i_stride_exp))

            # if unroll/jam is not needed to transform this loop (i.e. unroll factor == 1)
            if i_ufactor == 1:
                stmt.stmt = self.__transformLoop(stmt.stmt, orig_free_unrolls) 
                bound_unrolls.reverse()
                if self.use_full_unroll:
                    stmts = [stmt]
                    for iname, ufactor, stride_exp in bound_unrolls:
                        n_stmts = []
                        for i in range(0, ufactor):
                            it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
                            increment_exp = module.loop.ast.BinOpExp(it,
                                                                     stride_exp.replicate(),
                                                                     module.loop.ast.BinOpExp.MUL)
                            increment_exp = self.cfolder.fold(increment_exp)
                            for s in stmts:
                                n_stmts.append(self.__addIdentWithExp(s.replicate(), iname,
                                                                      increment_exp)) 
                        stmts = n_stmts
                    if len(stmts) == 1:
                        return stmts[0]
                    return module.loop.ast.CompStmt(stmts)
                else:
                    nstmt = stmt
                    for iname, ufactor, stride_exp in bound_unrolls:
                        intratile_iname = self.itvar_map_rev[iname]
                        loop_body = self.clib.replaceIdent(nstmt, iname, intratile_iname)
                        nstmt = self.__createIntraTileLoop(intratile_iname, iname, ufactor,
                                                           stride_exp, loop_body)
                    return nstmt

            # to apply the simple (rectangular) register tiling to this loop
            elif len(bound_unrolls) == 0 or (not self.use_nonrect_tile):
                main_loop, cleanup_loop = self.__createUnrolledLoops(i_for_loop_info, i_ufactor)
                main_loop.stmt = self.__transformLoop(main_loop.stmt, ext_free_unrolls)
                cleanup_loop.stmt = self.__transformLoop(cleanup_loop.stmt, orig_free_unrolls)
                bound_unrolls.reverse()
                if self.use_full_unroll:
                    stmts = [[main_loop, cleanup_loop]]
                    for iname, ufactor, stride_exp in bound_unrolls:
                        n_stmts = []
                        for i in range(0, ufactor):
                            it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
                            increment_exp = module.loop.ast.BinOpExp(it,
                                                                     stride_exp.replicate(),
                                                                     module.loop.ast.BinOpExp.MUL)
                            increment_exp = self.cfolder.fold(increment_exp)
                            for ss in stmts:
                                n_stmts.append([self.__addIdentWithExp(s.replicate(), iname,
                                                                       increment_exp)
                                                for s in ss])
                        stmts = n_stmts
                    stmts = reduce(lambda x,y: x+y, stmts, [])
                    return module.loop.ast.CompStmt(stmts)
                else:
                    nstmt = module.loop.ast.CompStmt([main_loop, cleanup_loop])
                    for iname, ufactor, stride_exp in bound_unrolls:
                        intratile_iname = self.itvar_map_rev[iname]
                        loop_body = self.clib.replaceIdent(nstmt, iname, intratile_iname)
                        nstmt = self.__createIntraTileLoop(intratile_iname, iname, ufactor,
                                                           stride_exp, loop_body)
                    return nstmt

            # to apply the non-rectangular register tiling to this loop
            else:
                
                # reverse the bound unroll factors
                bound_unrolls_rev = bound_unrolls[:]
                bound_unrolls_rev.reverse()

                # the new lower/upper bounds
                new_lb_id = module.loop.ast.IdentExp(self.__newlb_prefix + i_index_id.name)
                new_ub_id = module.loop.ast.IdentExp(self.__newub_prefix + i_index_id.name)

                # get bound unroll factors for the lower/upper bounds
                lb_bound_unrolls_rev = []
                ub_bound_unrolls_rev = []
                for u in bound_unrolls_rev:
                    iname, ufactor, stride_exp = u
                    if self.clib.containIdentName(i_lbound_exp, iname):
                        lb_bound_unrolls_rev.append(u)
                    if self.clib.containIdentName(i_ubound_exp, iname):
                        ub_bound_unrolls_rev.append(u)

                # create initialization statements for the upper/lower bounds
                if len(lb_bound_unrolls_rev) == 0:
                    a_val = i_lbound_exp.replicate()
                else:
                    a_val = module.loop.ast.NumLitExp(self.__int_min,
                                                      module.loop.ast.NumLitExp.INT)
                lb_asgn = module.loop.ast.BinOpExp(new_lb_id.replicate(),
                                                   a_val,
                                                   module.loop.ast.BinOpExp.EQ_ASGN)
                lb_init = module.loop.ast.ExpStmt(lb_asgn)
                if len(ub_bound_unrolls_rev) == 0:
                    a_val = i_ubound_exp.replicate()
                else:
                    a_val = module.loop.ast.NumLitExp(self.__int_max,
                                                      module.loop.ast.NumLitExp.INT)
                ub_asgn = module.loop.ast.BinOpExp(new_ub_id.replicate(),
                                                   a_val,
                                                   module.loop.ast.BinOpExp.EQ_ASGN)
                ub_init = module.loop.ast.ExpStmt(ub_asgn)
                
                # create the loop for finding the new lower/upper bounds
                lhs = new_lb_id.replicate()
                rhs = module.loop.ast.FunCallExp(module.loop.ast.IdentExp('max'),
                                                 [new_lb_id.replicate(), i_lbound_exp.replicate()])
                lb_asgn = module.loop.ast.BinOpExp(lhs, rhs, module.loop.ast.BinOpExp.EQ_ASGN)
                lb_stmt = module.loop.ast.ExpStmt(lb_asgn)
                if len(lb_bound_unrolls_rev) == 0:
                    lb_stmt = None
                lhs = new_ub_id.replicate()
                rhs = module.loop.ast.FunCallExp(module.loop.ast.IdentExp('min'),
                                                 [new_ub_id.replicate(), i_ubound_exp.replicate()])
                ub_asgn = module.loop.ast.BinOpExp(lhs, rhs, module.loop.ast.BinOpExp.EQ_ASGN)
                ub_stmt = module.loop.ast.ExpStmt(ub_asgn)
                if len(ub_bound_unrolls_rev) == 0:
                    ub_stmt = None
                stmts = []
                if lb_stmt: stmts += [lb_stmt]
                if ub_stmt: stmts += [ub_stmt]
                if len(stmts) == 0:
                    print 'internal error: loop-bound initialization loop cannot be empty'
                    sys.exit(1)
                if len(stmts) == 1:
                    bound_stmt = stmts[0]
                else:
                    bound_stmt = module.loop.ast.CompStmt(stmts)
                for iname, ufactor, stride_exp in bound_unrolls_rev:
                    intratile_iname = self.itvar_map_rev[iname]
                    loop_body = self.clib.replaceIdent(bound_stmt, iname, intratile_iname)
                    bound_stmt = self.__createIntraTileLoop(intratile_iname, iname, ufactor,
                                                            stride_exp, loop_body)

                # create the statement for computing the new lower/upper bounds
                bound_stmt = module.loop.ast.CompStmt([lb_init, ub_init, bound_stmt])

                # create the main and cleanup loops
                c_for_loop_info = (i_index_id, new_lb_id, new_ub_id, i_stride_exp, i_loop_body)
                main_loop, cleanup_loop = self.__createUnrolledLoops(c_for_loop_info, i_ufactor)
                main_loop.stmt = self.__transformLoop(main_loop.stmt, ext_outer_unrolls)
                cleanup_loop.stmt = self.__transformLoop(cleanup_loop.stmt, orig_outer_unrolls)

                # create the prologue and epilogue loops
                ub_exp = module.loop.ast.BinOpExp(new_lb_id.replicate(),
                                                  i_stride_exp.replicate(),
                                                  module.loop.ast.BinOpExp.SUB)
                prol_loop = self.flib.createForLoop(i_index_id, i_lbound_exp, ub_exp,
                                                    i_stride_exp, i_loop_body)
                prol_loop.stmt = self.__transformLoop(prol_loop.stmt, orig_free_unrolls)
                lb_exp = module.loop.ast.BinOpExp(new_ub_id.replicate(),
                                                  i_stride_exp.replicate(),
                                                  module.loop.ast.BinOpExp.ADD)
                epil_loop = self.flib.createForLoop(i_index_id, lb_exp, i_ubound_exp,
                                                    i_stride_exp, i_loop_body)
                epil_loop.stmt = self.__transformLoop(epil_loop.stmt, orig_free_unrolls)
                prol_stmt = None
                epil_stmt = None
                if self.use_full_unroll:
                    prol_stmts = [prol_loop]
                    epil_stmts = [epil_loop]
                    for iname, ufactor, stride_exp in bound_unrolls_rev:
                        n_prol_stmts = []
                        n_epil_stmts = []
                        for i in range(0, ufactor):
                            it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
                            increment_exp = module.loop.ast.BinOpExp(it,
                                                                     stride_exp.replicate(),
                                                                     module.loop.ast.BinOpExp.MUL)
                            increment_exp = self.cfolder.fold(increment_exp)
                            for s in prol_stmts:
                                n_prol_stmts.append(self.__addIdentWithExp(s.replicate(), iname,
                                                                           increment_exp))
                            for s in epil_stmts:
                                n_epil_stmts.append(self.__addIdentWithExp(s.replicate(), iname,
                                                                           increment_exp))
                        prol_stmts = n_prol_stmts
                        epil_stmts = n_epil_stmts
                    if len(prol_stmts) == 1:
                        prol_stmt = prol_stmts[0]
                    else:
                        prol_stmt = module.loop.ast.CompStmt(prol_stmts)
                    if len(epil_stmts) == 1:
                        epil_stmt = epil_stmts[0]
                    else:
                        epil_stmt = module.loop.ast.CompStmt(epil_stmts)
                else:
                    prol_stmt = prol_loop
                    epil_stmt = epil_loop
                    for iname, ufactor, stride_exp in bound_unrolls_rev:
                        intratile_iname = self.itvar_map_rev[iname]
                        prol_loop_body = self.clib.replaceIdent(prol_stmt, iname, intratile_iname)
                        epil_loop_body = self.clib.replaceIdent(epil_stmt, iname, intratile_iname)
                        prol_stmt = self.__createIntraTileLoop(intratile_iname, iname, ufactor,
                                                               stride_exp, prol_loop_body)
                        epil_stmt = self.__createIntraTileLoop(intratile_iname, iname, ufactor,
                                                               stride_exp, epil_loop_body)

                # combine together the loop-bound initializations, prologue loop, main loop,
                # cleanup loop, and the epilogue loop
                stmts = []
                stmts += bound_stmt.stmts
                if isinstance(prol_stmt, module.loop.ast.CompStmt):
                    stmts += prol_stmt.stmts
                else:
                    stmts += [prol_stmt]
                stmts += [main_loop, cleanup_loop]
                if isinstance(epil_stmt, module.loop.ast.CompStmt):
                    stmts += epil_stmt.stmts
                else:
                    stmts += [epil_stmt]
                return module.loop.ast.CompStmt(stmts)
            
        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)

        elif isinstance(stmt, module.loop.ast.NewAST):
            return stmt

        else:
            print 'internal error: unexpected AST type: "%s"' % stmt.__class__.__name__
            sys.exit(1)

    #----------------------------------------------------------

    def transform(self):
        '''To apply register tiling'''

        # to check the semantics of the input statement
        semant.SemanticChecker(self.stmt).check()

        # to traverse the abstract syntax tree to unroll/jam each encountered for-loop
        transformed_stmt = self.__transformLoop(self.stmt, [])

        # return the transformed statement
        return transformed_stmt


        
    
