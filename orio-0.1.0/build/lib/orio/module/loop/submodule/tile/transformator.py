#
# Contain the transformation procedure
#

import sys
import module.loop.ast, module.loop.ast_lib.constant_folder, module.loop.ast_lib.forloop_lib

#-----------------------------------------

class Transformator:
    '''Code transformator'''

    def __init__(self, tsize, tindex, stmt):
        '''To instantiate a code transformator object'''

        self.tsize = tsize
        self.tindex = tindex
        self.stmt = stmt
        
        self.flib = module.loop.ast_lib.forloop_lib.ForLoopLib()
        self.cfolder = module.loop.ast_lib.constant_folder.ConstFolder()
        
    #----------------------------------------------------------

    def transform(self):
        '''
        To tile the given loop structure. The resulting tiled loop will look like this.
        
        for (ii=LB; ii<=UB; ii+=Ti)
          for (i=ii; i<=min(UB,ii+Ti-ST); i+=ST)
            <loop-body>

        Such tiled code avoids the uses of division/multiplication operations. However, some
        compilers (e.g. ICC -fast) cannot finish its compilation when multiple levels of tiling
        are applied.
        '''

        # get rid of compound statement that contains only a single statement
        while isinstance(self.stmt, module.loop.ast.CompStmt) and len(self.stmt.stmts) == 1:
            self.stmt = self.stmt.stmts[0]
                                
        # extract for-loop structure
        for_loop_info = self.flib.extractForLoopInfo(self.stmt)
        index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info

        # check the tile index name
        if self.tindex == index_id.name:
            print (('error:%s: the tile index name must be different from the new tiled ' +
                    'loop index name: "%s"') % (index_id.line_no, self.tindex))
            sys.exit(1)
        
        # when tile size = 1, no transformation will be applied
        if self.tsize == 1:
            return self.flib.createForLoop(index_id, lbound_exp, ubound_exp,
                                           stride_exp, loop_body)

        # evaluate the stride expression
        try:
            stride_val = eval(str(stride_exp))
        except Exception, e:
            print ('error:%s: failed to evaluate expression: "%s"' %
                   (stride_exp.line_no, stride_exp))
            print ' --> %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)
        if not isinstance(stride_val, int) or stride_val <= 0:
            print ('error:%s: loop stride size must be a positive integer: %s' %
                   (stride_exp.line_no, stride_exp))
            sys.exit(1)

        # check whether tile_size % stride == 0
        if self.tsize % stride_val != 0:
            print ('error:%s: tile size (%s) must be divisible by the stride value (%s)'
                   % (stride_exp.line_no, self.tsize, stride_val))
            sys.exit(1)

        # create the tile index name
        tindex_id = module.loop.ast.IdentExp(self.tindex)

        # for the inter-tiling loop (i.e. outer loop)
        # compute lower bound --> LB' = LB
        tile_lbound_exp = lbound_exp.replicate()
        
        # compute upper bound --> UB' = UB
        tile_ubound_exp = ubound_exp.replicate()
        
        # compute stride --> ST' = tsize
        tile_stride_exp = module.loop.ast.NumLitExp(self.tsize, module.loop.ast.NumLitExp.INT) 

        # for the intra-tile loop (i.e. inner loop)
        # compute lower bound --> LB' = tindex
        itile_lbound_exp = tindex_id.replicate()

        # compute upper bound --> UB' = min(UB, tindex+tsize-ST)
        it1 = module.loop.ast.BinOpExp(module.loop.ast.NumLitExp(self.tsize,
                                                                 module.loop.ast.NumLitExp.INT),
                                       stride_exp.replicate(),
                                       module.loop.ast.BinOpExp.SUB)
        it2 = module.loop.ast.BinOpExp(tindex_id.replicate(), it1, module.loop.ast.BinOpExp.ADD)
        it2 = self.cfolder.fold(it2)
        itile_ubound_exp = module.loop.ast.FunCallExp(module.loop.ast.IdentExp('min'),
                                                      [ubound_exp.replicate(), it2])

        # compute stride --> ST' = ST
        itile_stride_exp = stride_exp.replicate()

        # generate the transformed statement
        transformed_stmt = self.flib.createForLoop(tindex_id, tile_lbound_exp,
                                                   tile_ubound_exp, tile_stride_exp,
                                                   self.flib.createForLoop(index_id,
                                                                           itile_lbound_exp,
                                                                           itile_ubound_exp,
                                                                           itile_stride_exp,
                                                                           loop_body))

        # return the transformed statement
        return transformed_stmt
             
