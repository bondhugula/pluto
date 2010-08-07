#
# Contain the transformation procedure
#

import sys
import module.loop.ast, module.loop.ast_lib.common_lib, module.loop.ast_lib.forloop_lib
import module.loop.submodule.tile.tile
import module.loop.submodule.permut.permut
import module.loop.submodule.regtile.regtile
import module.loop.submodule.scalarreplace.scalarreplace
import module.loop.submodule.boundreplace.boundreplace
import module.loop.submodule.pragma.pragma
import module.loop.submodule.arrcopy.arrcopy

#-----------------------------------------

class Transformator:
    '''Code transformator'''

    def __init__(self, tiles, permuts, regtiles, scalarrep, boundrep,
                 pragma, openmp, vector, arrcopy, stmt):
        '''To instantiate a code transformator object'''

        self.tiles = tiles
        self.permuts = permuts
        self.regtiles = regtiles
        self.scalarrep = scalarrep
        self.boundrep = boundrep
        self.pragma = pragma
        self.openmp = openmp
        self.vector = vector
        self.arrcopy = arrcopy
        self.stmt = stmt

        self.counter = 1
        self.prefix = 'cbv_'
        
        self.flib = module.loop.ast_lib.forloop_lib.ForLoopLib()
        self.clib = module.loop.ast_lib.common_lib.CommonLib()

        self.tile_smod = module.loop.submodule.tile.tile.Tile()
        self.perm_smod = module.loop.submodule.permut.permut.Permut()
        self.regt_smod = module.loop.submodule.regtile.regtile.RegTile()
        self.srep_smod = module.loop.submodule.scalarreplace.scalarreplace.ScalarReplace()
        self.brep_smod = module.loop.submodule.boundreplace.boundreplace.BoundReplace()
        self.prag_smod = module.loop.submodule.pragma.pragma.Pragma()
        self.acop_smod = module.loop.submodule.arrcopy.arrcopy.ArrCopy()

    #----------------------------------------------------------

    def __tile(self, stmt, tinfo):
        '''To apply loop tiling'''
        
        if isinstance(stmt, module.loop.ast.ExpStmt):
            return stmt

        elif isinstance(stmt, module.loop.ast.CompStmt):
            stmt.stmts = [self.__tile(s, tinfo) for s in stmt.stmts]
            return stmt
            
        elif isinstance(stmt, module.loop.ast.IfStmt):
            stmt.true_stmt = self.__tile(stmt.true_stmt, tinfo)
            if stmt.false_stmt:
                stmt.false_stmt = self.__tile(stmt.false_stmt, tinfo)
            return stmt
                
        elif isinstance(stmt, module.loop.ast.ForStmt):

            # recursively transform the loop body
            stmt.stmt = self.__tile(stmt.stmt, tinfo)

            # apply tiling if this is the loop to be tiled
            for_loop_info = self.flib.extractForLoopInfo(stmt)
            index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
            lid, tsize, tindex = tinfo
            if lid == index_id.name:
                stmt = self.tile_smod.tile(tsize, tindex, stmt)

            # return this loop statement
            return stmt

        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)
                                    
        elif isinstance(stmt, module.loop.ast.NewAST):
            return stmt

        else:
            print 'internal error: unexpected AST type: "%s"' % stmt.__class__.__name__
            sys.exit(1)
        
    #----------------------------------------------------------
    
    def __insertPragmas(self, stmt, tinfo):
        '''To insert pragma directives'''
        
        if isinstance(stmt, module.loop.ast.ExpStmt):
            return stmt

        elif isinstance(stmt, module.loop.ast.CompStmt):
            nstmts = []
            for s in stmt.stmts:
                is_comp_before = isinstance(s, module.loop.ast.CompStmt)
                ns = self.__insertPragmas(s, tinfo)
                is_comp_after = isinstance(ns, module.loop.ast.CompStmt)
                if not is_comp_before and is_comp_after:
                    nstmts.extend(ns.stmts)
                else:
                    nstmts.append(ns)
            stmt.stmts = nstmts
            return stmt
            
        elif isinstance(stmt, module.loop.ast.IfStmt):
            stmt.true_stmt = self.__insertPragmas(stmt.true_stmt, tinfo)
            if stmt.false_stmt:
                stmt.false_stmt = self.__insertPragmas(stmt.false_stmt, tinfo)
            return stmt
                
        elif isinstance(stmt, module.loop.ast.ForStmt):

            # recursively transform the loop body
            stmt.stmt = self.__insertPragmas(stmt.stmt, tinfo)

            # apply tiling if this is the loop to be tiled
            for_loop_info = self.flib.extractForLoopInfo(stmt)
            index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
            lid, pragmas = tinfo
            if lid == index_id.name:
                stmt = self.prag_smod.insertPragmas(pragmas, stmt)

            # return this loop statement
            return stmt

        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)
                                    
        elif isinstance(stmt, module.loop.ast.NewAST):
            return stmt

        else:
            print 'internal error: unexpected AST type: "%s"' % stmt.__class__.__name__
            sys.exit(1)
        
    #----------------------------------------------------------

    def __replaceBoundsInsertPrags(self, stmt, pragmas):
        '''Replace loop bounds with scalars, then insert pragmas right before the loop'''

        # get the loop structure
        for_loop_info = self.flib.extractForLoopInfo(stmt)
        index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info

        # replace complex loop bounds with scalars
        decls = []
        asgns = []
        if lbound_exp and self.clib.isComplexExp(lbound_exp):
            intmd = module.loop.ast.IdentExp(self.prefix + str(self.counter))
            self.counter += 1
            stmt.init.rhs = intmd.replicate()
            decls.append(module.loop.ast.VarDecl('register int', [intmd.name]))
            asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                            lbound_exp.replicate(),
                                            module.loop.ast.BinOpExp.EQ_ASGN)
            asgn = module.loop.ast.ExpStmt(asgn)
            asgns.append(asgn)
        if ubound_exp and self.clib.isComplexExp(ubound_exp):
            intmd = module.loop.ast.IdentExp(self.prefix + str(self.counter))
            self.counter += 1
            stmt.test.rhs = intmd.replicate()
            if len(decls) > 0:
                decls[0].var_names.append(intmd.name)
            else:
                decls.append(module.loop.ast.VarDecl('register int', [intmd.name]))
            asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                            ubound_exp.replicate(),
                                            module.loop.ast.BinOpExp.EQ_ASGN)
            asgn = module.loop.ast.ExpStmt(asgn)
            asgns.append(asgn)      
        
        # generate the transformed loop
        stmt = self.prag_smod.insertPragmas(pragmas, stmt)
        stmt.stmts = decls + asgns + stmt.stmts
        
        # return this loop statement
        return stmt

    #----------------------------------------------------------
    
    def __insertOpenMPPragmas(self, stmt, tinfo):
        '''To insert OpenMP pragma directives (on outermost loops only)'''
        
        if isinstance(stmt, module.loop.ast.ExpStmt):
            return stmt

        elif isinstance(stmt, module.loop.ast.CompStmt):
            nstmts = []
            for s in stmt.stmts:
                is_comp_before = isinstance(s, module.loop.ast.CompStmt)
                ns = self.__insertOpenMPPragmas(s, tinfo)
                is_comp_after = isinstance(ns, module.loop.ast.CompStmt)
                if not is_comp_before and is_comp_after:
                    nstmts.extend(ns.stmts)
                else:
                    nstmts.append(ns)
            stmt.stmts = nstmts
            return stmt
            
        elif isinstance(stmt, module.loop.ast.IfStmt):
            stmt.true_stmt = self.__insertOpenMPPragmas(stmt.true_stmt, tinfo)
            if stmt.false_stmt:
                stmt.false_stmt = self.__insertOpenMPPragmas(stmt.false_stmt, tinfo)
            return stmt
                
        elif isinstance(stmt, module.loop.ast.ForStmt):

            # get the loop structure
            for_loop_info = self.flib.extractForLoopInfo(stmt)
            index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
            
            # check if the initialization, test, and iteration variables exist
            if lbound_exp == None or ubound_exp == None or stride_exp == None:
                return stmt

            # replace loop bounds with scalars, and then insert openMP pragmas before the loop
            pragmas, = tinfo
            stmt = self.__replaceBoundsInsertPrags(stmt, pragmas)

            # return the transformed loop
            return stmt

        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)
                                    
        elif isinstance(stmt, module.loop.ast.NewAST):
            return stmt

        else:
            print 'internal error: unexpected AST type: "%s"' % stmt.__class__.__name__
            sys.exit(1)
        
    #----------------------------------------------------------
        
    def __insertVectorPragmas(self, stmt, tinfo):
        '''To insert vectorization pragma directives (on innermost loops only)'''

        if isinstance(stmt, module.loop.ast.ExpStmt):
            return stmt

        elif isinstance(stmt, module.loop.ast.CompStmt):
            nstmts = []
            for s in stmt.stmts:
                is_comp_before = isinstance(s, module.loop.ast.CompStmt)
                ns = self.__insertVectorPragmas(s, tinfo)
                is_comp_after = isinstance(ns, module.loop.ast.CompStmt)
                if not is_comp_before and is_comp_after:
                    nstmts.extend(ns.stmts)
                else:
                    nstmts.append(ns)
            stmt.stmts = nstmts
            return stmt
        
        elif isinstance(stmt, module.loop.ast.IfStmt):
            stmt.true_stmt = self.__insertVectorPragmas(stmt.true_stmt, tinfo)
            if stmt.false_stmt:
                stmt.false_stmt = self.__insertVectorPragmas(stmt.false_stmt, tinfo)
            return stmt
    
        elif isinstance(stmt, module.loop.ast.ForStmt):

            # apply recursion on the loop body
            stmt.stmt = self.__insertVectorPragmas(stmt.stmt, tinfo)

            # no transformation if it is not the inner most loop
            if self.flib.hasInnerLoop(stmt.stmt):
                return stmt

            # replace loop bounds with scalars, and then insert openMP pragmas before the loop
            pragmas, = tinfo
            stmt = self.__replaceBoundsInsertPrags(stmt, pragmas)

            # return the transformed loop
            return stmt

        elif isinstance(stmt, module.loop.ast.TransformStmt):
            print 'internal error: unprocessed transform statement'
            sys.exit(1)
                                    
        elif isinstance(stmt, module.loop.ast.NewAST):
            return stmt

        else:
            print 'internal error: unexpected AST type: "%s"' % stmt.__class__.__name__
            sys.exit(1)
        
    #----------------------------------------------------------
    
    def __searchLoopId(self, all_lids, input_lid):
        '''
        Given all the existing loop index names, search if there is a match with the given loop ID.
        Complain necessarily if there is no match found.
        '''

        for i in input_lid[1:]:
            if i in all_lids:
                return i
        must_exist = input_lid[0]
        if must_exist:
            print ('error: no matching loop index name in input argument: %s' %
                   (tuple(input_lid[1:]), ))
            sys.exit(1)
        else:
            return None
            
    #----------------------------------------------------------

    def transform(self):
        '''To apply the composite transformations'''

        # copy the statement
        tstmt = self.stmt.replicate()

        # reset counter (for variable name generation)
        self.counter = 1
        
        # apply loop tiling
        for loop_id, tsize, tindex in self.tiles:
            all_lids = self.flib.getLoopIndexNames(tstmt)
            lid = self.__searchLoopId(all_lids, loop_id)
            if lid != None:
                tinfo = (lid, tsize, tindex)
                tstmt = self.__tile(tstmt, tinfo)

        # apply loop permutation/interchange
        for seq in self.permuts:
            tstmt = self.perm_smod.permute(seq, tstmt)

        # apply array-copy optimization
        for do_acopy, aref, suffix, dtype, dimsizes in self.arrcopy:
            if not do_acopy:
                dimsizes = [1] * len(dimsizes)
            tstmt = self.acop_smod.optimizeArrayCopy(aref, suffix, dtype, dimsizes, tstmt)

        # apply register tiling
        loops, ufactors = self.regtiles
        if len(loops) > 0:
            tstmt = self.regt_smod.tileForRegs(loops, ufactors, tstmt)

        # apply scalar replacement
        do_scalarrep, dtype, prefix = self.scalarrep
        if do_scalarrep:
            tstmt = self.srep_smod.replaceScalars(dtype, prefix, tstmt)
        
        # apply bound replacement
        do_boundrep, lprefix, uprefix = self.boundrep
        if do_boundrep:
            tstmt = self.brep_smod.replaceBounds(lprefix, uprefix, tstmt)

        # insert pragma directives
        for loop_id, pragmas in self.pragma:
            all_lids = self.flib.getLoopIndexNames(tstmt)
            lid = self.__searchLoopId(all_lids, loop_id)
            if lid != None:
                tinfo = (lid, pragmas)
                tstmt = self.__insertPragmas(tstmt, tinfo)

        # insert openmp directives (apply only on outermost loops)
        do_openmp, pragmas = self.openmp
        if do_openmp:
            tinfo = (pragmas, )
            tstmt = self.__insertOpenMPPragmas(tstmt, tinfo)

        # insert vectorization directives (apply only on innermost loops)
        do_vector, pragmas = self.vector
        if do_vector:
            tinfo = (pragmas, )
            tstmt = self.__insertVectorPragmas(tstmt, tinfo)
            
        # return the transformed statement
        return tstmt


