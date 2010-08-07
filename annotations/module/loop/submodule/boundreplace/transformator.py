#
# Contain the transformation procedure
#

import sys
import module.loop.ast

#-----------------------------------------

__counter1 = 1    # for variable naming
__counter2 = 1    # for variable naming
__data_type = None

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

def __replaceLoopBounds(stmt, iter_name):
    '''Replace complex (lower and upper) loop bounds with a scalar'''

    if not stmt:
        return (stmt, [])

    if isinstance(stmt, module.loop.ast.ExpStmt):
        return (stmt, [])

    elif isinstance(stmt, module.loop.ast.CompStmt):
        stmts = []
        bound_inits = []
        for s in stmt.stmts:
            is_comp_before = isinstance(s, module.loop.ast.CompStmt)
            (ns, nbound_inits) = __replaceLoopBounds(s, iter_name)
            is_comp_after = isinstance(ns, module.loop.ast.CompStmt)
            if not is_comp_before and is_comp_after:
                stmts.extend(ns.stmts)
            else:
                stmts.append(ns)
            bound_inits.extend(nbound_inits)
        stmt.stmts = stmts
        return (stmt, bound_inits)
            
    elif isinstance(stmt, module.loop.ast.IfStmt):
        ntrue_stmt, ntrue_bound_inits = __replaceLoopBounds(stmt.true_stmt, iter_name)
        nfalse_stmt, nfalse_bound_inits = __replaceLoopBounds(stmt.false_stmt, iter_name)
        stmt.true_stmt = ntrue_stmt
        stmt.false_stmt = nfalse_stmt
        return (stmt, ntrue_bound_inits + nfalse_bound_inits)
        
    elif isinstance(stmt, module.loop.ast.ForStmt):
        # check initialization expression
        if (stmt.init and
            (not isinstance(stmt.init, module.loop.ast.BinOpExp) or
             stmt.init.op_type != module.loop.ast.BinOpExp.EQ_ASGN or
             not isinstance(stmt.init.lhs, module.loop.ast.IdentExp))):
            print 'error:%s: initialization expression not in (ID=LB) form' % stmt.init.line_no
            sys.exit(1)
            
        # check test expression
        if (stmt.test and
            (not isinstance(stmt.test, module.loop.ast.BinOpExp) or
             stmt.test.op_type != module.loop.ast.BinOpExp.LE or
             not isinstance(stmt.test.lhs, module.loop.ast.IdentExp))):
            print 'error:%s: test expression not in (ID<=UB) form' % stmt.test.line_no
            sys.exit(1)

        # get iterator name
        niter_name = None
        if stmt.init:
            niter_name = stmt.init.lhs.name
        if stmt.test:
            if niter_name and niter_name != stmt.test.lhs.name:
                print ('error:%s: different iterator names across init and test exps'
                       % stmt.line_no)
                sys.exit(1)
            niter_name = stmt.test.lhs.name
        if not niter_name:
            print 'error:%s: both init and test exps of a loop cannot be empty' % stmt.line_no
            sys.exit(1)
            
        # recursion on the loop body
        nstmt, nbound_inits = __replaceLoopBounds(stmt.stmt, niter_name)
        stmt.stmt = nstmt
        
        # get bounds
        lbound_exp = None
        ubound_exp = None
        if stmt.init:
            lbound_exp = stmt.init.rhs.replicate()
        if stmt.test:
            ubound_exp = stmt.test.rhs.replicate()
                            
        # replace bounds if necessary
        global __counter1
        global __counter2
        bound_inits = nbound_inits[:]
        if lbound_exp and not (isinstance(lbound_exp, module.loop.ast.IdentExp) or
                               isinstance(lbound_exp, module.loop.ast.NumLitExp)):
            intmd = module.loop.ast.IdentExp('lbv' + str(__counter1))
            __counter1 += 1
            stmt.init.rhs = intmd.replicate()
            if __data_type:
                decl = module.loop.ast.VarDecl(__data_type, [intmd.replicate()])
            else:
                decl = None
            asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                            lbound_exp.replicate(),
                                            module.loop.ast.BinOpExp.EQ_ASGN)
            asgn = module.loop.ast.ExpStmt(asgn)
            bound_inits.append((decl, asgn))
        if ubound_exp and not (isinstance(ubound_exp, module.loop.ast.IdentExp) or
                               isinstance(ubound_exp, module.loop.ast.NumLitExp)):
            intmd = module.loop.ast.IdentExp('ubv' + str(__counter2))
            __counter2 += 1
            stmt.test.rhs = intmd.replicate()
            if __data_type:
                decl = module.loop.ast.VarDecl(__data_type, [intmd.replicate()])
            else:
                decl = None
            asgn = module.loop.ast.BinOpExp(intmd.replicate(),
                                            ubound_exp.replicate(),
                                            module.loop.ast.BinOpExp.EQ_ASGN)
            asgn = module.loop.ast.ExpStmt(asgn)
            bound_inits.append((decl, asgn))

        # generate transformed statement
        decls = []
        asgns = []
        nbound_inits = []
        for decl, asgn in bound_inits:
            bound_exp = asgn.exp.rhs
            if __containIdentName(bound_exp, iter_name):
                decls.append(decl)
                asgns.append(asgn)
            else:
                nbound_inits.append((decl, asgn))
        while decls.count(None) > 0:
            decls.remove(None)
        if len(decls) > 0 or len(asgns) > 0:
            nstmt = module.loop.ast.CompStmt(decls + asgns + [stmt])
        else:
            nstmt = stmt

        return (nstmt, nbound_inits)
                   
    elif isinstance(stmt, module.loop.ast.TransformStmt):
        print 'internal error: unprocessed transform statement'
        sys.exit(1)
        
    elif isinstance(stmt, module.loop.ast.NewAST):
        return (stmt, [])
        
    else:
        print 'internal error: unexpected AST type: "%s"' % stmt.__class__.__name__
        sys.exit(1)

#-----------------------------------------

def transform(stmt, arg_info):
    '''Perform code transformation'''
    
    # extract input information
    dtype, = arg_info
    
    # update global variables
    global __counter1
    global __counter2
    global __data_type
    __counter1 = 1
    __counter2 = 1
    __data_type = dtype

    # generate the transformed statement
    (transformed_stmt, bound_inits) = __replaceLoopBounds(stmt, None)
    if len(bound_inits) > 0:
        decls, asgns = zip(*bound_inits)
        decls = list(decls)
        while decls.count(None) > 0:
            decls.remove(None)
        asgns = list(asgns)
        if isinstance(transformed_stmt, module.loop.ast.CompStmt):
            transformed_stmt.stmts = decls + asgns + transformed_stmt.stmts
        else:
            transformed_stmt = module.loop.ast.CompStmt(decls + asgns + [transformed_stmt])
        
    return transformed_stmt

