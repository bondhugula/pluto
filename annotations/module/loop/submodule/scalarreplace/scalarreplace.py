#
# Scalar replacement transformation.
#
# For example, if we want to apply scalar replacement on a perfectly nested loop (shown below):
#
#  /* The original loop */
#  for I1
#   for I2
#    ...
#    for In
#    {
#      S1(I1,I2,...,In);
#      S2(I1,I2,...,In);
#      ...
#      Sn(I1,I2,...,In);
#    }
#
# we can annotate the loop as follows:
#
#  /* The annotated loop */
#  transform ScalarReplace(prefix="it", dtype="double")
#  for I1
#   for I2
#    ...
#    for In
#    {
#      S1(I1,I2,...,In);
#      S2(I1,I2,...,In);
#      ...
#      Sn(I1,I2,...,In);
#    }
#
# First, the scalar replacer identifies common references to arrays and functions, and then
# stores common references in intermediates for future reuse. The place to create
# intermediates is somewhere between the loops and determined by the referenced iterators and
# loop structure.
#
# Notes:
#  1. 'prefix' (intermediate prefix name) is the intermediate name used as a prefix, followed
#     by a generated number (starting from 1).
#  2. 'dtype' (data type) is the data type used in the intermediate variable declaration.
#  3. The loop structure is assumed to be perfectly nested. If not, code generation might
#     include errors.
#

import sys
import module.loop.ast, module.loop.submodule.submodule, transformator

#-----------------------------------------

class ScalarReplace(module.loop.submodule.submodule.SubModule):
    
    def __init__(self):
        '''To instantiate a scalar replacement submodule'''
        module.loop.submodule.submodule.SubModule.__init__(self)

    def transform(self, trans_stmt):
        '''To transform the "transform" statement'''

        # initialize argument list
        args = {'prefix' : None,   # string
                'dtype' : None}    # string

        # read all arguments
        for k in trans_stmt.kw_args:
            kname = k.lhs.name
            if kname not in args:
                print 'error:%s: unknown argument: "%s"' % (k.lhs.line_no, kname)
                sys.exit(1)
            args[kname] = k.rhs

        # extract for-loop structure information
        for_loop_info = extractPerfectlyNestedForLoopInfo(trans_stmt.stmt)

        # check arguments and update the values of some arguments accordingly
        for kname, exp in args.items():
            if not exp:
                print 'error:%s: missing argument: "%s"' % (trans_stmt.line_no, kname)
                sys.exit(1)
            if kname in ('prefix', 'dtype'):
                if not isinstance(exp, module.loop.ast.StringLitExp):
                    print ('error:%s: prefix and data-type must have string-typed values'
                           % exp.line_no)
                    sys.exit(1)
                args[kname] = str(exp.val[1:len(exp.val)-1])
        if not args['prefix']:
            print 'error:%s: prefix name cannot be an empty string' % trans_stmt.line_no
            sys.exit(1)

        # create argument information
        arg_info = [args['prefix'], args['dtype']]

        # perform transformation
        transformed_stmt = transformator.transform(for_loop_info, arg_info)

        return transformed_stmt

#-----------------------------------------

def extractPerfectlyNestedForLoopInfo(stmt):
    '''Given a for-loop statement, extract information about its loop structure'''

    # get rid of compound statement that contains only a single statement
    while isinstance(stmt, module.loop.ast.CompStmt) and len(stmt.stmts) == 1:
        stmt = stmt.stmts[0]

    # check if it is a for-loop statement
    if not isinstance(stmt, module.loop.ast.ForStmt):
        print 'error:%s: not a for-loop statement' % stmt.line_no
        sys.exit(1)

    return __extractPerfectlyNestedForLoopInfo(stmt)

def __extractPerfectlyNestedForLoopInfo(stmt):
    '''Given a for-loop statement, extract information about its loop structure'''

    # get rid of compound statement that contains only a single statement
    while isinstance(stmt, module.loop.ast.CompStmt) and len(stmt.stmts) == 1:
        stmt = stmt.stmts[0]

    # in case of a non-for-loop statement
    if not isinstance(stmt, module.loop.ast.ForStmt):
        return ([], stmt)

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

    # check iteration expression
    if (stmt.iter and
        (not isinstance(stmt.iter, module.loop.ast.BinOpExp) or
         stmt.iter.op_type != module.loop.ast.BinOpExp.EQ_ASGN or
         not isinstance(stmt.iter.lhs, module.loop.ast.IdentExp) or
         not isinstance(stmt.iter.rhs, module.loop.ast.BinOpExp) or
         stmt.iter.rhs.op_type != module.loop.ast.BinOpExp.ADD or
         not isinstance(stmt.iter.rhs.lhs, module.loop.ast.IdentExp) or
         stmt.iter.lhs.name != stmt.iter.rhs.lhs.name)):
        if (not isinstance(stmt.iter, module.loop.ast.UnaryExp) or
            stmt.iter.op_type not in (module.loop.ast.UnaryExp.POST_INC,
                                      module.loop.ast.UnaryExp.PRE_INC) or
            not isinstance(stmt.iter.exp, module.loop.ast.IdentExp)):
            print ('error:%s: iteration expression not in (ID++ or ID+=ST or ID=ID+ST) form'
                   % stmt.iter.line_no)
            sys.exit(1)

    # check if the control expressions are all empty
    if not stmt.init and not stmt.test and not stmt.iter:
        print ('error:%s: a loop with no control expression (i.e. infinite loop) is not allowed'
               % stmt.line_no)
        sys.exit(1)
    
    # check if the iterators are the same
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
        print ('error:%s: different iterator names across init, test, and iter exps'
               % stmt.line_no)
        sys.exit(1)
        
    # extract control information
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
        else:
            assert(isinstance(stmt.iter, module.loop.ast.UnaryExp)), 'internal error: not unary'
            stride_exp = module.loop.ast.NumLitExp(1, module.loop.ast.NumLitExp.INT)
    control_info = (index_id, lbound_exp, ubound_exp, stride_exp)

    # recursively extract information from the nested statement
    control_infos, loop_body = __extractPerfectlyNestedForLoopInfo(stmt.stmt)
    control_infos = [control_info] + control_infos
    
    return (control_infos, loop_body)



