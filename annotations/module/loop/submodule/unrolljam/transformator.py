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

def __addIdentifierWithExp(tnode, index_name, exp):
    '''Traverse the tree node and add any matching identifier with the provided expression'''

    if isinstance(tnode, module.loop.ast.ExpStmt):
        if tnode.exp:
            tnode.exp = __addIdentifierWithExp(tnode.exp, index_name, exp)
        return tnode
    
    elif isinstance(tnode, module.loop.ast.CompStmt):
        tnode.stmts = [__addIdentifierWithExp(s, index_name, exp) for s in tnode.stmts]
        return tnode

    elif isinstance(tnode, module.loop.ast.IfStmt):
        tnode.test = __addIdentifierWithExp(tnode.test, index_name, exp)
        tnode.true_stmt = __addIdentifierWithExp(tnode.true_stmt, index_name, exp)
        if tnode.false_stmt:
            __addIdentifierWithExp(tnode.false_stmt, index_name, exp)
        return tnode

    elif isinstance(tnode, module.loop.ast.ForStmt):
        if tnode.init:
            tnode.init = __addIdentifierWithExp(tnode.init, index_name, exp)
        if tnode.test:
            tnode.test = __addIdentifierWithExp(tnode.test, index_name, exp)
        if tnode.iter:
            tnode.iter = __addIdentifierWithExp(tnode.iter, index_name, exp)
        tnode.stmt = __addIdentifierWithExp(tnode.stmt, index_name, exp)
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
        elif isinstance(exp, module.loop.ast.NumLitExp) and exp.val == 0:
            return tnode
        else:
            add_exp = module.loop.ast.BinOpExp(tnode,
                                               exp.replicate(),
                                               module.loop.ast.BinOpExp.ADD)
            return module.loop.ast.ParenthExp(add_exp)

    elif isinstance(tnode, module.loop.ast.ArrayRefExp):
        tnode.exp = __addIdentifierWithExp(tnode.exp, index_name, exp)
        tnode.sub_exp = __addIdentifierWithExp(tnode.sub_exp, index_name, exp)
        return tnode

    elif isinstance(tnode, module.loop.ast.FunCallExp):
        tnode.exp = __addIdentifierWithExp(tnode.exp, index_name, exp)
        tnode.args = [__addIdentifierWithExp(a, index_name, exp) for a in tnode.args]
        return tnode

    elif isinstance(tnode, module.loop.ast.UnaryExp):
        tnode.exp = __addIdentifierWithExp(tnode.exp, index_name, exp)
        return tnode

    elif isinstance(tnode, module.loop.ast.BinOpExp):
        tnode.lhs = __addIdentifierWithExp(tnode.lhs, index_name, exp)
        tnode.rhs = __addIdentifierWithExp(tnode.rhs, index_name, exp)
        return tnode

    elif isinstance(tnode, module.loop.ast.ParenthExp):
        tnode.exp = __addIdentifierWithExp(tnode.exp, index_name, exp)
        return tnode        

    elif isinstance(tnode, module.loop.ast.NewAST):
        return tnode

    else:
        print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
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

def __jamStmts(s1, s2):
    '''Jam/fuse statements s1 and s2 whenever possible'''

    if not s1:
        return s2
    if not s2:
        return s1

    while isinstance(s1, module.loop.ast.CompStmt) and len(s1.stmts) == 1:
        s1 = s1.stmts[0]
    while isinstance(s2, module.loop.ast.CompStmt) and len(s2.stmts) == 1:
        s2 = s2.stmts[0]

    if (isinstance(s1, module.loop.ast.ForStmt) and
        isinstance(s2, module.loop.ast.ForStmt) and
        str(s1.init) == str(s2.init) and
        str(s1.test) == str(s2.test) and
        str(s1.iter) == str(s2.iter)):
        return module.loop.ast.ForStmt(s1.init, s1.test, s1.iter, __jamStmts(s1.stmt, s2.stmt))

    elif ((not isinstance(s1, module.loop.ast.CompStmt)) and
          (not isinstance(s2, module.loop.ast.CompStmt))):
        return module.loop.ast.CompStmt([s1,s2])
        
    else:
        stmts = []
        if isinstance(s1, module.loop.ast.CompStmt):
            stmts.extend(s1.stmts)
        else:
            stmts.append(s1)
        if isinstance(s2, module.loop.ast.CompStmt):
            s = __jamStmts(stmts.pop(), s2.stmts[0])
            if isinstance(s, module.loop.ast.CompStmt):
                stmts.extend(s.stmts)
            else:
                stmts.append(s)
            stmts.extend(s2.stmts[1:])
        else:
            s = __jamStmts(stmts.pop(), s2)
            if isinstance(s, module.loop.ast.CompStmt):
                stmts.extend(s.stmts)
            else:
                stmts.append(s)
        if len(stmts) == 1:
            return stmts[0]
        else:
            return module.loop.ast.CompStmt(stmts)

#-----------------------------------------

def transform(for_loop_info, arg_info):
    '''Perform code transformation'''

    # extract input information
    index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
    ufactor, = arg_info 

    # get rid of compound statement that contains only a single statement
    while isinstance(loop_body, module.loop.ast.CompStmt) and len(loop_body.stmts) == 1:
        loop_body = loop_body.stmts[0]

    # when ufactor = 1, no transformation will be applied
    if ufactor == 1:
        return __makeForLoop(index_id, lbound_exp, ubound_exp, stride_exp, loop_body)
                
    # start generating the main unrolled loop
    # compute lower bound --> LB' = LB
    new_lbound_exp = lbound_exp.replicate()
    
    # compute upper bound --> UB' = UB-ST*(UF-1)
    it = module.loop.ast.BinOpExp(stride_exp.replicate(),
                                  module.loop.ast.NumLitExp(ufactor - 1,
                                                            module.loop.ast.NumLitExp.INT),
                                  module.loop.ast.BinOpExp.MUL)
    new_ubound_exp = module.loop.ast.BinOpExp(ubound_exp.replicate(),
                                              it,
                                              module.loop.ast.BinOpExp.SUB)
    new_ubound_exp = __simpleConstantFolding(new_ubound_exp)
    
    # compute stride --> ST' = UF*ST
    it = module.loop.ast.NumLitExp(ufactor, module.loop.ast.NumLitExp.INT)
    new_stride_exp = module.loop.ast.BinOpExp(it,
                                              stride_exp.replicate(),
                                              module.loop.ast.BinOpExp.MUL)
    new_stride_exp = __simpleConstantFolding(new_stride_exp)
    
    # compute unrolled statements
    unrolled_stmts = []
    for i, s in enumerate([loop_body.replicate() for i in range(0, ufactor)]):
        if i == 0:
            if isinstance(s, module.loop.ast.CompStmt):
                unrolled_stmts.append(s.stmts)
            else:
                unrolled_stmts.append([s])
        else:
            if i == 1:
                increment_exp = stride_exp.replicate()
            else:
                it = module.loop.ast.NumLitExp(i, module.loop.ast.NumLitExp.INT)
                increment_exp = module.loop.ast.BinOpExp(stride_exp.replicate(),
                                                         it,
                                                         module.loop.ast.BinOpExp.MUL)
                increment_exp = __simpleConstantFolding(increment_exp)
            ns = __addIdentifierWithExp(s, index_id.name, increment_exp)
            if isinstance(ns, module.loop.ast.CompStmt):
                unrolled_stmts.append(ns.stmts)
            else:
                unrolled_stmts.append([ns])

    # compute the unrolled loop body by jamming/fusing the unrolled statements
    jammed_stmts = []
    for stmts in zip(*unrolled_stmts):
        jammed_stmts.append(reduce(lambda x,y: __jamStmts(x,y), stmts, None))
    unrolled_loop_body = reduce(lambda x,y: __jamStmts(x,y), jammed_stmts, None)

    # generate the main unrolled loop
    main_loop = __makeForLoop(index_id, new_lbound_exp, new_ubound_exp, new_stride_exp,
                              unrolled_loop_body)

    # generate the clean-up loop
    cleanup_loop = __makeForLoop(index_id, None, ubound_exp, stride_exp, loop_body)

    # generate the transformer statement
    transformed_stmt = module.loop.ast.CompStmt([main_loop, cleanup_loop])

    return transformed_stmt

