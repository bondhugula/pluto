#
# Bound replacement transformation.
#
# This transformation will replace complex (lower and upper) bounds of for-loops with
# intermediates.

# For example, the following loop:
#
#  /* The original loop */
#  for (i=max(M,N); i<=min(X,Y); i+=ST)
#   S(i)
#
# can be annotated as follows:
#
#  /* The annotated loop */
#  transform BoundReplace()
#  for (i=max(M,N); i<=min(X,Y); i+=ST)
#   S(i)
#
# to replace the loop bounds, resulting in the following loop.
#
#  /* Bound-replaced loop */
#  double LB, UB;
#  LB = max(M,N);
#  UB = min(X,Y);
#  for (i=LB; i<=UB; i+=ST)
#   S(i)
#
# Notes:
#  1. 'dtype' (data type) is the data type used in the intermediate variable declaration.
#

import sys
import module.loop.ast, module.loop.submodule.submodule, transformator

#-----------------------------------------

class BoundReplace(module.loop.submodule.submodule.SubModule):
    
    def __init__(self):
        '''To instantiate a bound replacement submodule'''
        module.loop.submodule.submodule.SubModule.__init__(self)

    def transform(self, trans_stmt):
        '''To transform the "transform" statement'''

        # initialize argument list
        args = {'dtype' : None}   # string
        
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
            if kname == 'dtype':
                if not isinstance(exp, module.loop.ast.StringLitExp):
                    print ('error:%s: data-type must have string-typed values'
                           % exp.line_no)
                    sys.exit(1)
                args[kname] = str(exp.val[1:len(exp.val)-1])
                                                            
        # create argument information
        arg_info = [args['dtype']]

        # perform transformation
        transformed_stmt = transformator.transform(trans_stmt.stmt, arg_info)

        return transformed_stmt

