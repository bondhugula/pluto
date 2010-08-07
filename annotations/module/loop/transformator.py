#
# Contain the transformation procedure
#

import sys
import ast

#-----------------------------------------

__loaded_submodules = {}
__submod_directory = 'module.loop.submodule'

#-----------------------------------------

def __checkKeywordArguments(trans_stmt):
    '''Check if the transform statement has repeated keyword argument'''

    seen_knames = []
    for k in trans_stmt.kw_args:
        kname = k.lhs.name
        if kname in seen_knames:
            print 'error:%s: repeated argument names: "%s"' % (trans_stmt.line_no, kname)
            sys.exit(1)
        seen_knames.append(kname)

#-----------------------------------------

def __loadAndApplyTransformation(trans_stmt):
    '''Load the transformation submodule and then apply the transformation'''

    # load the corresponding transformation submodule
    submod_name = trans_stmt.name
    submodule = None
    if __loaded_submodules.has_key(submod_name.lower()):
        submodule = __loaded_submodules[submod_name.lower()]
    else:
        submod_fname = __submod_directory + '.' + submod_name.lower() + '.' + submod_name.lower()
        try:
            submodule = __import__(submod_fname)
            components = submod_fname.split('.')
            for c in components[1:]:
                submodule = getattr(submodule, c)
        except Exception, e:
            print ('error: submodule "%s" does not exist (failed to load "%s")' %
                   (submod_name, submod_fname))
            print ' --> cause: %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)
        try:
            getattr(submodule, submod_name)
        except:
            print 'error: no class "%s" defined in "%s"' %s (submod_name, submod_fname)
            sys.exit(1)
        __loaded_submodules[submod_name.lower()] = submodule

    # load the transformation class from the loaded submodule
    submod_class = getattr(submodule, submod_name)

    # call transformation procedure
    try:
        transformed_stmt = submod_class().transform(trans_stmt)
    except Exception, e:
        print '%s: %s' % (e.__class__.__name__, e)
        sys.exit(1)

    return transformed_stmt

#-----------------------------------------

def __transform(stmt):
    '''Traverse each statement and apply transformation procedure. The transformation
    procedure is specified in the "transform" statement.'''

    if isinstance(stmt, ast.ExpStmt):
        return stmt

    elif isinstance(stmt, ast.CompStmt):
        stmt.stmts = [__transform(s) for s in stmt.stmts]
        return stmt

    elif isinstance(stmt, ast.IfStmt):
        stmt.true_stmt = __transform(stmt.true_stmt)
        if stmt.false_stmt:
            stmt.false_stmt = __transform(stmt.false_stmt)
        return stmt

    elif isinstance(stmt, ast.ForStmt):
        stmt.stmt = __transform(stmt.stmt)
        return stmt

    elif isinstance(stmt, ast.TransformStmt):
        stmt.stmt = __transform(stmt.stmt)
        __checkKeywordArguments(stmt)
        return __loadAndApplyTransformation(stmt)

    else:
        print 'internal error: unknown statement type: %s' % stmt.__class__.__name__
        sys.exit(1)

#-----------------------------------------

def transform(stmts):
    '''Go through each statement to apply transformation procedure'''

    # apply transformation procedure on each statement
    transformed_stmts = [__transform(s) for s in stmts]
        
    return transformed_stmts

