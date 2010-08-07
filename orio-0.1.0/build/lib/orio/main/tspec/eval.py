# 
# The evaluator for the TSpec (Tuning Specifier) language
#

import StringIO, sys, tokenize

#--------------------------------------------------------------------------------

class TSpecEvaluator:
    '''The evaluator for the TSpec language'''

    def __init__(self):
        '''To instantiate a TSpec evaluator'''
        pass
    
    #----------------------------------------------------------------------------

    def __extractVars(self, code):
        '''Return all variables that are present in the given code'''

        # tokenize the given expression code
        gtoks = tokenize.generate_tokens(StringIO.StringIO(code).readline)

        # iterate over each token and replace any matching token with its corresponding value
        vnames = []
        for toknum, tokval, _, _, _ in gtoks:
            if toknum == tokenize.NAME:
                vnames.append(tokval)

        # return all the found variable names
        return vnames

    #----------------------------------------------------------------------------

    def __substituteVars(self, code, env):
        '''
        Expand any variables that exist in the given environment to their corresponding values
        '''

        # tokenize the given expression code
        gtoks = tokenize.generate_tokens(StringIO.StringIO(code).readline)

        # iterate over each token and replace any matching token with its corresponding value
        tokens = []
        for toknum, tokval, _, _, _ in gtoks:
            if toknum == tokenize.NAME and tokval in env:
                ntoks = tokenize.generate_tokens(StringIO.StringIO(str(env[tokval])).readline)
                tokens.extend(ntoks)
            else:
                tokens.append((toknum, tokval))

        # convert the tokens back to a string
        code = tokenize.untokenize(tokens)

        # remove all the leading and trailing spaces
        code = code.strip()

        # return the modified string
        return code

    #----------------------------------------------------------------------------

    def __evalArg(self, stmt, env, name_space):
        '''To evaluate the given "let" statement'''

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), (rhs, rhs_line_no) = stmt

        # check for illegal variable references
        for vname in self.__extractVars(rhs):
            try:
                eval(vname, env)
            except:
                print 'error:%s: invalid reference: "%s"' % (rhs_line_no, vname)
                sys.exit(1)

        # evaluate the RHS expression
        try:
            rhs_val = eval(rhs, env)
        except Exception, e:
            print 'error:%s: failed to evaluate the RHS expression' % rhs_line_no
            print ' --> %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)

        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), (rhs_val, rhs_line_no))

    #----------------------------------------------------------------------------

    def __evalConstraint(self, stmt, env, name_space):
        '''To evaluate the given "constraint" statement'''

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), (rhs, rhs_line_no) = stmt
        
        # substitute all environment variables with their corresponding values
        rhs = self.__substituteVars(rhs, env)

        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), (rhs, rhs_line_no))

    #----------------------------------------------------------------------------

    def __evalDecl(self, stmt, env, name_space):
        '''To evaluate the given "decl" statement'''

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), type_seq, dim_exp_seq, (rhs, rhs_line_no) = stmt
        
        # check for types
        type_names = []
        for t, l in type_seq:
            if t in type_names:
                print 'error:%s: repeated type name: "%s"' % (l, t)
                sys.exit(1)
            type_names.append(t)

        # substitute all environment variables in each dimension expression
        n_dim_exp_seq = []
        for e, l in dim_exp_seq:
            e = self.__substituteVars(e, env)
            n_dim_exp_seq.append((e, l))
        dim_exp_seq = n_dim_exp_seq

        # substitute all environment variables in the RHS expression
        if rhs and rhs != 'random':
            rhs = self.__substituteVars(rhs, env)

        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), type_seq, dim_exp_seq, (rhs, rhs_line_no))

    #----------------------------------------------------------------------------

    def __evalDef(self, stmt, env, name_space):
        '''To evaluate the given "def" statement'''

        # copy the environment and name space
        env = env.copy()
        name_space = name_space.copy()

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), stmt_seq = stmt
        
        # evaluate each statement in the definition statement body
        stmt_seq = self.__evaluate(stmt_seq, env, name_space)
        
        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), stmt_seq)

    #----------------------------------------------------------------------------

    def __evalLet(self, stmt, env, name_space):
        '''To evaluate the given "let" statement'''

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), (rhs, rhs_line_no) = stmt

        # check for illegal variable references
        for vname in self.__extractVars(rhs):
            try:
                eval(vname, env)
            except:
                print 'error:%s: invalid reference: "%s"' % (rhs_line_no, vname)
                sys.exit(1)

        # evaluate the RHS expression
        try:
            rhs_val = eval(rhs, env)
        except Exception, e:
            print 'error:%s: failed to evaluate the RHS expression' % rhs_line_no
            print ' --> %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)

        # update the environment
        env[id_name] = rhs_val

        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), (rhs_val, rhs_line_no))

    #----------------------------------------------------------------------------

    def __evalParam(self, stmt, env, name_space):
        '''To evaluate the given "param" statement'''

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), is_range, (rhs, rhs_line_no) = stmt

        # check for illegal variable references
        for vname in self.__extractVars(rhs):
            try:
                eval(vname, env)
            except:
                print 'error:%s: invalid reference: "%s"' % (rhs_line_no, vname)
                sys.exit(1)

        # evaluate the RHS expression
        try:
            rhs_val = eval(rhs, env)
        except Exception, e:
            print 'error:%s: failed to evaluate the RHS expression' % rhs_line_no
            print ' --> %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)

        # check the RHS value
        if is_range:
            if not isinstance(rhs_val, list) and not isinstance(rhs_val, tuple):
                print 'error:%s: RHS must be a list/tuple' % rhs_line_no
                sys.exit(1)
            if len(rhs_val) == 0:
                print 'error:%s: RHS must not be an empty list' % rhs_line_no
                sys.exit(1)
            etype = type(rhs_val[0])
            for e in rhs_val:
                if not isinstance(e, etype):
                    print 'error:%s: RHS must be a list of equal-typed elements' % rhs_line_no
                    sys.exit(1)
        
        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), is_range, (rhs_val, rhs_line_no))
        
    #----------------------------------------------------------------------------

    def __evalSpec(self, stmt, env, name_space):
        '''To evaluate the given "spec" statement'''

        # copy the environment and name space
        env = env.copy()
        name_space = name_space.copy()

        # unpack the statement
        keyw, line_no, (id_name, id_line_no), stmt_seq = stmt

        # evaluate each statement in the specification statement body
        stmt_seq = self.__evaluate(stmt_seq, env, name_space)

        # return the evaluated statement
        return (keyw, line_no, (id_name, id_line_no), stmt_seq)

    #----------------------------------------------------------------------------

    def __evaluate(self, stmt, env, name_space):
        '''
        To evaluate the given statement. Note that the given statement could be a statement sequence.
        '''

        # in the case of a single statement
        if isinstance(stmt, tuple):

            # get keyword, line number, and identifier name
            keyw = stmt[0]
            line_no = stmt[1]
            (id_name, id_line_no) = stmt[2]

            # check for any predefined name
            if id_name in name_space:
                print 'error:%s: name "%s" already defined' % (id_line_no, id_name)
                sys.exit(1)

            # first update the name space before evaluation (if necessary)
            if keyw in ('def', 'spec'):
                name_space[id_name] = keyw
    
            # evaluate each statement
            if keyw == 'arg':
                e = self.__evalArg(stmt, env, name_space)
            elif keyw == 'constraint':
                e = self.__evalConstraint(stmt, env, name_space)
            elif keyw == 'decl':
                e = self.__evalDecl(stmt, env, name_space)
            elif keyw == 'def':
                e = self.__evalDef(stmt, env, name_space)
            elif keyw == 'let':
                e = self.__evalLet(stmt, env, name_space)
            elif keyw == 'param':
                e = self.__evalParam(stmt, env, name_space)
            elif keyw == 'spec':
                e = self.__evalSpec(stmt, env, name_space)
            else:
                print 'internal error:%s: unrecognized TSpec statement' % line_no
                sys.exit(1)

            # update the name_space
            name_space[id_name] = keyw

            # return the evaluated statement
            return e
            
        # in the case of a sequence of statements
        elif isinstance(stmt, list):

            # evaluate each statement
            e = [self.__evaluate(s, env, name_space) for s in stmt]

            # return the evaluated statement sequence
            return e

        # unexpected input
        else:
            print 'internal error: unexpected type of TSpec statement'
            sys.exit(1)
            
    #----------------------------------------------------------------------------

    def evaluate(self, stmt_seq):
        '''To evaluate the given statement sequence'''
        return self.__evaluate(stmt_seq, {}, {})



