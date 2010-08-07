#
# Contain the transformation procedure
#

import sys
import module.loop.ast

#-----------------------------------------

def __simpleConstantFolding(exp):
    '''Perform a simple (i.e. limited) constant folding'''
    
    if isinstance(exp, module.loop.ast.NumLitExp):
        return exp

    elif isinstance(exp, module.loop.ast.StringLitExp):
        return exp

    elif isinstance(exp, module.loop.ast.IdentExp):
        return exp

    elif isinstance(exp, module.loop.ast.ArrayRefExp):
        return exp

    elif isinstance(exp, module.loop.ast.FunCallExp):
        return exp

    elif isinstance(exp, module.loop.ast.UnaryExp):
        exp.exp = __simpleConstantFolding(exp.exp)
        if isinstance(exp.exp, module.loop.ast.NumLitExp):
            if exp.op_type == module.loop.ast.UnaryExp.PLUS:
                return exp.exp
            elif exp.op_type == module.loop.ast.UnaryExp.MINUS:
                exp.exp.val = - exp.exp.val
                return exp.exp
            else:
                return exp
        else:
            return exp

    elif isinstance(exp, module.loop.ast.BinOpExp):
        exp.lhs = __simpleConstantFolding(exp.lhs)
        exp.rhs = __simpleConstantFolding(exp.rhs)
        if exp.op_type in (module.loop.ast.BinOpExp.MUL,
                           module.loop.ast.BinOpExp.DIV,
                           module.loop.ast.BinOpExp.MOD,
                           module.loop.ast.BinOpExp.ADD,
                           module.loop.ast.BinOpExp.SUB):
            if (isinstance(exp.lhs, module.loop.ast.NumLitExp) and
                isinstance(exp.rhs, module.loop.ast.NumLitExp)):
                new_val = None
                if exp.op_type == module.loop.ast.BinOpExp.MUL:
                    new_val = exp.lhs.val * exp.rhs.val
                elif exp.op_type == module.loop.ast.BinOpExp.DIV:
                    new_val = 1.0 * exp.lhs.val / exp.rhs.val
                elif exp.op_type == module.loop.ast.BinOpExp.MOD:
                    new_val = exp.lhs.val % exp.rhs.val
                elif exp.op_type == module.loop.ast.BinOpExp.ADD:
                    new_val = exp.lhs.val + exp.rhs.val
                elif exp.op_type == module.loop.ast.BinOpExp.SUB:
                    new_val = exp.lhs.val - exp.rhs.val
                else:
                    print 'internal error: unexpected bin-op operation'
                    sys.exit(1)
                new_lit_type = None
                if isinstance(new_val, int):
                    new_lit_type = module.loop.ast.NumLitExp.INT
                elif isinstance(new_val, float):
                    new_lit_type = module.loop.ast.NumLitExp.FLOAT
                else:
                    print 'internal error: unexpected numeric type'
                    sys.exit(1)
                return module.loop.ast.NumLitExp(new_val, new_lit_type)
            else:
                return exp
        else:
            return exp

    elif isinstance(exp, module.loop.ast.ParenthExp):
        exp.exp = __simpleConstantFolding(exp.exp)
        if isinstance(exp.exp, module.loop.ast.NumLitExp):
            return exp.exp
        else:
            return exp

    else:
        print 'internal error: unexpected AST type: "%s"' % exp.__class__.__name__
        sys.exit(1)
    
#-----------------------------------------

def __makeForLoop(id, lbound, ubound, stride, loop_body):
    '''Generate a for loop:
         for (id=lbound; id<=ubound; id=id+stride)
           loop_body'''

    init_exp = None
    test_exp = None
    iter_exp = None
    if lbound:
        init_exp = module.loop.ast.BinOpExp(id.replicate(),
                                            lbound.replicate(),
                                            module.loop.ast.BinOpExp.EQ_ASGN)
    if ubound:
        test_exp = module.loop.ast.BinOpExp(id.replicate(),
                                            ubound.replicate(),
                                            module.loop.ast.BinOpExp.LE)
    if stride:
        it = module.loop.ast.BinOpExp(id.replicate(),
                                      stride.replicate(),
                                      module.loop.ast.BinOpExp.ADD)
        iter_exp = module.loop.ast.BinOpExp(id.replicate(),
                                            it,
                                            module.loop.ast.BinOpExp.EQ_ASGN)
    return module.loop.ast.ForStmt(init_exp, test_exp, iter_exp, loop_body.replicate())
    
#-----------------------------------------

def transform(for_loop_info, arg_info):
    '''Perform code transformation'''

    # extract input information
    index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
    tsize, tindex_id = arg_info
            
    # get rid of compound statement that contains only a single statement
    while isinstance(loop_body, module.loop.ast.CompStmt) and len(loop_body.stmts) == 1:
        loop_body = loop_body.stmts[0]
                
    # when tsize = 1, no transformation will be applied
    if tsize == 1:
        return __makeForLoop(index_id, lbound_exp, ubound_exp, stride_exp, loop_body)

    # for the tiling loop (i.e. outer loop)
    # compute lower bound --> LB' = LB 
    tile_lbound_exp = lbound_exp.replicate()

    # compute upper bound --> UB' = UB
    tile_ubound_exp = ubound_exp.replicate()

    # compute stride --> ST' = tsize
    tile_stride_exp = module.loop.ast.NumLitExp(tsize, module.loop.ast.NumLitExp.INT)

    # for the intra-tile loop (i.e. inner loop)
    # compute lower bound --> LB' = tindex
    itile_lbound_exp = tindex_id.replicate()

    # compute upper bound --> UB' = min(UB, tindex+tsize-ST)
    it1 = module.loop.ast.BinOpExp(module.loop.ast.NumLitExp(tsize,
                                                             module.loop.ast.NumLitExp.INT),
                                   stride_exp.replicate(),
                                   module.loop.ast.BinOpExp.SUB)
    it2 = module.loop.ast.BinOpExp(tindex_id.replicate(), it1, module.loop.ast.BinOpExp.ADD)
    it2 = __simpleConstantFolding(it2)
    itile_ubound_exp = module.loop.ast.FunCallExp(module.loop.ast.IdentExp('min'),
                                                  [ubound_exp.replicate(), it2])
    
    # compute stride --> ST' = ST
    itile_stride_exp = stride_exp.replicate()

    # generate the entire tiled loop
    transformed_stmt = __makeForLoop(tindex_id, tile_lbound_exp,
                                     tile_ubound_exp, tile_stride_exp,
                                     __makeForLoop(index_id, itile_lbound_exp,
                                                   itile_ubound_exp, itile_stride_exp,
                                                   loop_body))
    
    return transformed_stmt
                
