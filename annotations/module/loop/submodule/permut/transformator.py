#
# Contain the transformation procedure
#

import sys
import module.loop.ast

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

def transform(stmt, arg_info):
    '''Perform code transformation'''

    # extract argument information
    loop_order, = arg_info

    # get rid of compound statement that contains only a single statement
    while isinstance(stmt, module.loop.ast.CompStmt) and len(stmt.stmts) == 1:
        stmt = stmt.stmts[0]

    # insert loop order information into a hashtable
    loop_info = {}
    for index_name, is_optional in loop_order:
        loop_info[index_name] = [is_optional]

    # create loop order (get rid of all optionality information)
    loop_order = [iname for iname, opt in loop_order]

    # extract loop control information and get the loop body
    loop_body = None
    cur_stmt = stmt
    unseen_loops = loop_order[:]
    seen_loops = []
    while True:
        if isinstance(cur_stmt, module.loop.ast.CompStmt) and len(cur_stmt.stmts) == 1:
            cur_stmt = cur_stmt.stmts[0]
            continue
        
        is_optional_list = [loop_info[i][0] for i in unseen_loops]
        all_unseen_optional = reduce(lambda x,y: x and y, is_optional_list, True)
        
        if isinstance(cur_stmt, module.loop.ast.ForStmt) and not cur_stmt.init:
            print ('error:%s:Permut: a loop is assumed to have a non-empty init exp'
                   % (cur_stmt.line_no))
            sys.exit(1)

        if (isinstance(cur_stmt, module.loop.ast.ForStmt) and
            isinstance(cur_stmt.init, module.loop.ast.BinOpExp) and
            cur_stmt.init.op_type == module.loop.ast.BinOpExp.EQ_ASGN and
            isinstance(cur_stmt.init.lhs, module.loop.ast.IdentExp)):

            iname = cur_stmt.init.lhs.name
            if iname in seen_loops:
                if all_unseen_optional:
                    loop_body = cur_stmt
                    break
                else:
                    print ('error:%s: loop "%s" cannot occur repeatedly'
                           % (cur_stmt.line_no, iname))
                    sys.exit(1)
            if iname not in unseen_loops:
                if all_unseen_optional:
                    loop_body = cur_stmt
                    break
                else:
                    print ('error:%s: loop "%s" is not specified in the loop order %s'
                           % (cur_stmt.line_no, iname, tuple(loop_order)))
                    sys.exit(1)

            linfo = loop_info[iname]
            linfo.append(cur_stmt.init)
            linfo.append(cur_stmt.test)
            linfo.append(cur_stmt.iter)

            unseen_loops.remove(iname)
            seen_loops.append(iname)
            cur_stmt = cur_stmt.stmt

        else:
            if all_unseen_optional:
                loop_body = cur_stmt
                break
            else:
                unfound_loops = filter(lambda x: not loop_info[x][0], unseen_loops)
                unfound_loops = tuple(unfound_loops)
                print ('error:%s: to-be-permuted loops %s do not exist'
                       % (stmt.line_no, unfound_loops))
                sys.exit(1)

    # generate the permuted loop
    transformed_stmt = loop_body
    rev_loop_order = loop_order[:]
    rev_loop_order.reverse()
    for iname in rev_loop_order:
        linfo = loop_info[iname]
        if len(linfo) > 1:
            opt, init_exp, test_exp, iter_exp = linfo
            transformed_stmt = module.loop.ast.ForStmt(init_exp.replicate(),
                                                       test_exp.replicate(),
                                                       iter_exp.replicate(),
                                                       transformed_stmt)

    return transformed_stmt
