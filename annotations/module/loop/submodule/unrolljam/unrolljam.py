#
# Loop unroll and jam transformation.
#
# For example, the following loop:
#
#  /* The original loop */
#  for (i = LBi; i <= UBi; i += STi)
#  {
#    S1(i);
#    for (j = LBj; j <= UBj; j += STj)
#      S2(i,j);
#    S3(i);
#  }
#
# that can be annotated as follows:
#
#  /* The annotated loop */
#  transform UnrollJam(ufactor=N)
#  for (i = LBi; i <= UBi; i += STi)
#  {
#    S1(i);
#    for (j = LBj; j <= UBj; j += STj)
#      S2(i,j);
#    S3(i);
#  }
#
# is unrolled N times, and then its unrolled body is jammed/fused together whenever possible,
# resulting in the following loop.
#
#  /* The main unrolled loop */
#  for (i = LBi; i <= UBi-(N-1)*STi; i += N*STi)
#  {
#    S1(i);
#    S1(i + STi);
#    S1(i + 2*STi);
#    ...
#    S1(i + (N-1)*STi);
#    for (j = LBj; j <= UBj; j += STj)
#    {
#      S2(i, j);
#      S2(i + STi, j);
#      S2(i + 2*STi, j);
#      ...
#      S2(i + (N-1)*STi, j);
#    }
#    S3(i);
#    S3(i + STi);
#    S3(i + 2*STi);
#    ...
#    S3(i + (N-1)*STi);
#  }
#  /* The cleanup loop */
#  for (; i <= UBi; i += STi)
#  {
#    S1(i);
#    for (j = LBj; j <= UBj; j += STj)
#      S(i,j);
#    S3(i);
#  }
#
# Notes:
#  1. 'ufactor' (unroll factor) must be a positive integer (i.e. 1,2,3,4,5,...)
#  2. This unrolling and jamming work only for POSITIVE stride values. 
#  3. A simple constant folding is applied on the new ST, the new UB, and the increment
#     expression (i.e. (N-1)*ST).
#     For example, UB-(4-1)*2 becomes UB-6.
#  4. Jamming/fusing multiple loops is performed when LBj, UBj, and STj are free of iterators
#     being unrolled.
#

import sys
import module.loop.ast, module.loop.submodule.submodule, transformator

#-----------------------------------------

class UnrollJam(module.loop.submodule.submodule.SubModule):
    
    def __init__(self):
        '''To instantiate a loop unroll-and-jam submodule'''
        module.loop.submodule.submodule.SubModule.__init__(self)

    def transform(self, trans_stmt):
        '''To transform the "transform" statement'''

        # initialize argument list
        args = {'ufactor' : None}  # integer

        # read all arguments
        for k in trans_stmt.kw_args:
            kname = k.lhs.name
            if kname not in args:
                print 'error:%s: unknown argument: "%s"' % (k.lhs.line_no, kname)
                sys.exit(1)
            args[kname] = k.rhs

        # extract for-loop structure information
        for_loop_info = extractForLoopInfo(trans_stmt.stmt)
        index_id, lbound_exp, ubound_exp, stride_exp, loop_body = for_loop_info
        
        # check arguments and update the values of some arguments accordingly
        for kname, exp in args.items():
            if not exp:
                print 'error:%s: missing argument: "%s"' % (trans_stmt.line_no, kname)
                sys.exit(1)
            if kname == 'ufactor':
                if (not isinstance(exp, module.loop.ast.NumLitExp) or
                    exp.lit_type != module.loop.ast.NumLitExp.INT or
                    exp.val <= 0):
                    print ('error:%s: unrolling factor must be a positive integer: "%"'
                           % (exp.line_no, exp))
                    sys.exit(1)
                args[kname] = int(exp.val)

        # create argument information
        arg_info = [args['ufactor']]

        # perform transformation
        transformed_stmt = transformator.transform(for_loop_info, arg_info)

        return transformed_stmt

#-----------------------------------------

def extractForLoopInfo(stmt):
    '''Given a for-loop statement, extract information about its loop structure'''

    # get rid of compound statement that contains only a single statement
    while isinstance(stmt, module.loop.ast.CompStmt) and len(stmt.stmts) == 1:
        stmt = stmt.stmts[0]

    # check if it is a for-loop statement
    if not isinstance(stmt, module.loop.ast.ForStmt):
        print 'error:%s: not a for-loop statement' % stmt.line_no
        sys.exit(1)

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
    loop_body = stmt.stmt

    return (index_id, lbound_exp, ubound_exp, stride_exp, loop_body)

