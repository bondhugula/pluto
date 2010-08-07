#
# Loop permutation/interchange transformation.
#
# For example, to permute the following loop:
#
#  /* The original loop */
#  for (i = 1; i <= p; i += 1)
#   for (j = 1; j <= q; j += 1)
#    for (k = 1; k <= r; k += 1)
#      S(i,j,k);
#
# to j-k-i order, the loop can be annotated as follows:
#
#  /* The annotated loop */
#  transform Permut (order=(j,k,i))
#  for (i = 1; i <= p; i += 1)
#   for (j = 1; j <= q; j += 1)
#    for (k = 1; k <= r; k += 1)
#      S(i,j,k);
#
# The interchanged loop can be seen below.
#
#  /* The permuted loop */
#  for (j = 1; j <= q; j += 1)
#   for (k = 1; k <= r; k += 1)
#    for (i = 1; i <= p; i += 1)
#     S(i,j,k);
#
# Notes:
#  1. 'order' (loop order) is the loop order used to permute the loop.
#  2. Iterator name is used to represent a loop. When an iterator name is contained within
#     additional parentheses, it implies that the corresponding loop may not exist.
#     For example, order=(i,j,(l),k) means that loop 'l' will be placed accordingly, if it
#     exists. Otherwise, code transformation proceeds with no errors.
#

import sys
import module.loop.ast, module.loop.submodule.submodule, transformator

#-----------------------------------------

class Permut(module.loop.submodule.submodule.SubModule):
    
    def __init__(self):
        '''To instantiate a loop permutation submodule'''
        module.loop.submodule.submodule.SubModule.__init__(self)

    def transform(self, trans_stmt):
        '''To transform the "transform" statement'''

        # initialize argument list
        args = {'order' : None}   # list of tuples (index_name, is_optional)

        # read all arguments
        for k in trans_stmt.kw_args:
            kname = k.lhs.name
            if kname not in args:
                print 'error:%s: unknown argument: "%s"' % (k.lhs.line_no, kname)
                sys.exit(1)
            args[kname] = k.rhs

        # check arguments and update the values of some arguments accordingly
        for kname, exp in args.items():
            if not exp:
                print 'error:%s: missing argument: "%s"' % (trans_stmt.line_no, kname)
                sys.exit(1)
            if kname == 'order':
                args[kname] = extractLoopOrder(exp)
        index_names = [iname for iname, opt in args['order']]
        for i in index_names:
            if index_names.count(i) > 1:
                print ('error:%s: repeated index name in loop order: "%s"'
                       % (trans_stmt.line_no, i))
                sys.exit(1)
                
        # create argument information
        arg_info = [args['order']]
        
        # perform transformation
        transformed_stmt = transformator.transform(trans_stmt.stmt, arg_info)

        return transformed_stmt

#-----------------------------------------

def extractLoopOrder(exp):
    '''Check if the given expression is a valid loop order. If yes, scan the loop order
    and return a list of tuples of an index name and a boolean value that indicates if the
    related loop is optional.
    For instance, with a given expression of (i,j,(k),l), this function will return
    [('i', False), ('j', False), ('k', True), ('l', False)]'''

    if not isinstance(exp, module.loop.ast.ParenthExp):
        print ('error:%s: not a proper form of a loop order: "%s"' %
               (exp.line_no, exp))
        sys.exit(1)

    loop_order = []
    exps_stack = []
    exps_stack.insert(0, exp.exp)
    while len(exps_stack) > 0:
        cur_exp = exps_stack.pop(0)
        if isinstance(cur_exp, module.loop.ast.IdentExp):
            loop_order.append((cur_exp.name, False))
        elif (isinstance(cur_exp, module.loop.ast.ParenthExp) and
              isinstance(cur_exp.exp, module.loop.ast.IdentExp)):
            loop_order.append((cur_exp.exp.name, True))
        elif (isinstance(cur_exp, module.loop.ast.BinOpExp) and
              cur_exp.op_type == module.loop.ast.BinOpExp.COMMA):
            exps_stack.insert(0, cur_exp.rhs)
            exps_stack.insert(0, cur_exp.lhs)
        else:
            print ('error:%s: not a proper form of a loop order: "%s"' %
                   (exp.line_no, exp))
            sys.exit(1)

    return loop_order
        
    
