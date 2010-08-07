#
# Contain the transformation procedure
#

import sys
import module.loop.ast

#-----------------------------------------

def __containIdentName(exp, id_name):
    '''Check if the given expression contains an identifier with a name match '''

    if isinstance(exp, module.loop.ast.NumLitExp):
        return False
    
    elif isinstance(exp, module.loop.ast.StringLitExp):
        return False
    
    elif isinstance(exp, module.loop.ast.IdentExp):
        return exp.name == id_name
    
    elif isinstance(exp, module.loop.ast.ArrayRefExp):
        return __containIdentName(exp.exp, id_name) or __containIdentName(exp.sub_exp, id_name)
        
    elif isinstance(exp, module.loop.ast.FunCallExp):
        has_match = reduce(lambda x,y: x or y,
                           [__containIdentName(a, id_name) for a in exp.args],
                           False)
        return __containIdentName(exp.exp, id_name) or has_match

    elif isinstance(exp, module.loop.ast.UnaryExp):
        return __containIdentName(exp.exp, id_name)
    
    elif isinstance(exp, module.loop.ast.BinOpExp):
        return __containIdentName(exp.lhs, id_name) or __containIdentName(exp.rhs, id_name)
        
    elif isinstance(exp, module.loop.ast.ParenthExp):
        return __containIdentName(exp.exp, id_name)

    elif isinstance(exp, module.loop.ast.NewAST):
        return False
        
    else:
        print 'internal error: unexpected AST type: "%s"' % exp.__class__.__name__
        sys.exit(1)

#-----------------------------------------

def __replaceRefsWithIntermediates(tnode, ref_infos):
    '''Replace array references and function calls with intermediates if they occur repeatedly'''
    
    if not tnode:
        return tnode

    if isinstance(tnode, module.loop.ast.ExpStmt):
        tnode.exp = __replaceRefsWithIntermediates(tnode.exp, ref_infos)
        return tnode
        
    elif isinstance(tnode, module.loop.ast.CompStmt):
        tnode.stmts = [__replaceRefsWithIntermediates(s, ref_infos) for s in tnode.stmts]
        return tnode
        
    elif isinstance(tnode, module.loop.ast.IfStmt):
        tnode.test = __replaceRefsWithIntermediates(tnode.test, ref_infos)
        tnode.true_stmt = __replaceRefsWithIntermediates(tnode.true_stmt, ref_infos)
        tnode.false_stmt = __replaceRefsWithIntermediates(tnode.false_stmt, ref_infos)
        return tnode
        
    elif isinstance(tnode, module.loop.ast.ForStmt):
        print 'error:%s:ScalarReplace: loop must be perfectly nested' % tnode.line_no
        sys.exit(1)

    elif isinstance(tnode, module.loop.ast.TransformStmt):
        print 'internal error: unprocessed transform statement'
        sys.exit(1)
        
    elif isinstance(tnode, module.loop.ast.NumLitExp):
        return tnode
    
    elif isinstance(tnode, module.loop.ast.StringLitExp):
        return tnode
    
    elif isinstance(tnode, module.loop.ast.IdentExp):
        return tnode
    
    elif (isinstance(tnode, module.loop.ast.ArrayRefExp) or
          isinstance(tnode, module.loop.ast.FunCallExp)):
        rkey = str(tnode)
        if rkey in ref_infos:
            freq, is_output, exp, seq, intdm = ref_infos[rkey]
            return intdm.replicate()
        else:
            return tnode
        
    elif isinstance(tnode, module.loop.ast.UnaryExp):
        tnode.exp = __replaceRefsWithIntermediates(tnode.exp, ref_infos)
        return tnode
    
    elif isinstance(tnode, module.loop.ast.BinOpExp):
        tnode.lhs = __replaceRefsWithIntermediates(tnode.lhs, ref_infos)
        tnode.rhs = __replaceRefsWithIntermediates(tnode.rhs, ref_infos)
        return tnode
        
    elif isinstance(tnode, module.loop.ast.ParenthExp):
        tnode.exp = __replaceRefsWithIntermediates(tnode.exp, ref_infos)
        return tnode

    elif isinstance(tnode, module.loop.ast.NewAST):
        return tnode
        
    else:
        print 'internal error: unexpected AST type: "%s"' % tnode.__class__.__name__
        sys.exit(1)

#-----------------------------------------

def __collectRefInfos(tnode, ref_infos, is_output):
    '''Search the occurrences of array references and function calls in the given tree node'''

    if not tnode:
        return

    if isinstance(tnode, module.loop.ast.ExpStmt):
        __collectRefInfos(tnode.exp, ref_infos, False)
        
    elif isinstance(tnode, module.loop.ast.CompStmt):
        for s in tnode.stmts:
            __collectRefInfos(s, ref_infos, False)
        
    elif isinstance(tnode, module.loop.ast.IfStmt):
        __collectRefInfos(tnode.test, ref_infos, False)
        __collectRefInfos(tnode.true_stmt, ref_infos, False)
        __collectRefInfos(tnode.false_stmt, ref_infos, False)

    elif isinstance(tnode, module.loop.ast.ForStmt):
        print 'error:%s:ScalarReplace: loop must be perfectly nested' % tnode.line_no
        sys.exit(1)

    elif isinstance(tnode, module.loop.ast.TransformStmt):
        print 'internal error: unprocessed transform statement'
        sys.exit(1)
        
    elif isinstance(tnode, module.loop.ast.NumLitExp):
        return
    
    elif isinstance(tnode, module.loop.ast.StringLitExp):
        return
    
    elif isinstance(tnode, module.loop.ast.IdentExp):
        return
    
    elif (isinstance(tnode, module.loop.ast.ArrayRefExp) or
          isinstance(tnode, module.loop.ast.FunCallExp)):
        rkey = str(tnode)
        if rkey in ref_infos:
            old_freq, old_is_output, old_exp, old_seq = ref_infos[rkey]
            ref_infos[rkey] = (old_freq + 1, is_output or old_is_output, old_exp, old_seq)
        else:
            ref_infos[rkey] = (1, is_output, tnode.replicate(), len(ref_infos))
        
    elif isinstance(tnode, module.loop.ast.UnaryExp):
        __collectRefInfos(tnode.exp, ref_infos, is_output)
        
    elif isinstance(tnode, module.loop.ast.BinOpExp):
        if tnode.op_type == module.loop.ast.BinOpExp.EQ_ASGN:
            __collectRefInfos(tnode.lhs, ref_infos, True)
        else:
            __collectRefInfos(tnode.lhs, ref_infos, False)
        __collectRefInfos(tnode.rhs, ref_infos, False)
        
    elif isinstance(tnode, module.loop.ast.ParenthExp):
        __collectRefInfos(tnode.exp, ref_infos, is_output)

    elif isinstance(tnode, module.loop.ast.NewAST):
        return
        
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

def transform(for_loop_info, arg_info):
    '''Perform code transformation'''

    # extract input information
    control_infos, loop_body = for_loop_info
    prefix, dtype = arg_info

    # get rid of compound statement that contains only a single statement
    while isinstance(loop_body, module.loop.ast.CompStmt) and len(loop_body.stmts) == 1:
        loop_body = loop_body.stmts[0]

    # collect references information
    ref_infos = {}
    __collectRefInfos(loop_body, ref_infos, False)

    # get rid of references that occur only once in the loop body
    htable = {}
    for rkey, rinfo in ref_infos.iteritems():
        freq, is_output, exp, seq = rinfo
        if freq > 1:
            htable[rkey] = rinfo
    ref_infos = htable

    # create a new intermediate for each reference
    rinfo_list = ref_infos.items()
    rinfo_list.sort(lambda x,y: cmp(x[1][3],y[1][3]))
    htable = {}
    for i, (rkey, rinfo) in enumerate(rinfo_list):
        freq, is_output, exp, seq = rinfo
        intmd = module.loop.ast.IdentExp(prefix + str(i+1))
        htable[rkey] = (freq, is_output, exp, seq, intmd)
    ref_infos = htable

    # replace repeated references with intermediates
    loop_body = __replaceRefsWithIntermediates(loop_body, ref_infos)

    # generate the transformed statement
    rev_control_infos = control_infos[:]
    rev_control_infos.reverse()
    rev_control_infos.append(None)
    transformed_stmt = loop_body
    for cinfo in rev_control_infos:
        if cinfo:
            index_id, lbound_exp, ubound_exp, stride_exp = cinfo
        decls = []
        prelude_stmts = []
        closing_stmts = []
        rinfo_list = ref_infos.items()
        rinfo_list.sort(lambda x,y: cmp(x[1][3],y[1][3]))
        htable = {}
        for rkey, rinfo in rinfo_list:
            freq, is_output, exp, seq, intmd = rinfo
            if not cinfo or __containIdentName(exp, index_id.name):
                if dtype:
                    decl = module.loop.ast.VarDecl(dtype, [intmd.replicate()])
                    decls.append(decl)
                asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                                exp.replicate(),
                                                module.loop.ast.BinOpExp.EQ_ASGN)
                asgn = module.loop.ast.ExpStmt(asgn)
                prelude_stmts.append(asgn)
                if is_output:
                    asgn = module.loop.ast.BinOpExp(exp.replicate(),
                                                    intmd.replicate(),
                                                    module.loop.ast.BinOpExp.EQ_ASGN)
                    asgn = module.loop.ast.ExpStmt(asgn)
                    closing_stmts.append(asgn)
            else:
                htable[rkey] = rinfo    
        ref_infos = htable
        if isinstance(transformed_stmt, module.loop.ast.CompStmt):
            transformed_stmt.stmts = (decls + prelude_stmts + transformed_stmt.stmts
                                      + closing_stmts)
        else:
            stmts = decls + prelude_stmts + [transformed_stmt] + closing_stmts
            if len(stmts) == 1:
                transformed_stmt = stmts[0]
            else:
                transformed_stmt = module.loop.ast.CompStmt(stmts)
        if cinfo:
            transformed_stmt = __makeForLoop(index_id, lbound_exp, ubound_exp, stride_exp,
                                             transformed_stmt)
            
    return transformed_stmt









