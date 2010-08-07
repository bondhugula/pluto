# 
# A class to encompass a collection of transformation information used by the PolySyn module
#

import sys

#--------------------------------------------------------------------------------

class TransfInfo:
    '''The transformation information'''

    def __init__(self, parallel, tiles, permut, unroll_factors, scalar_replace, vectorize,
                 profiling_code, compile_cmd, compile_opts):
        '''To instantiate a transformation information instance'''

        self.parallel = parallel
        self.tiles = tiles
        self.permut = permut
        self.unroll_factors = unroll_factors
        self.scalar_replace = scalar_replace
        self.vectorize = vectorize

        self.profiling_code = profiling_code
        self.compile_cmd = compile_cmd
        self.compile_opts = compile_opts

#--------------------------------------------------------------------------------

class TransfInfoGen:
    '''The generator for transformation information'''

    def __init__(self):
        '''To instantiate an instance for the transformation information generator'''
        pass

    #----------------------------------------------------------------------------

    def generate(self, assigns, perf_params):
        '''Generate a complete transformation information from the given assignment sequence'''

        # check that no argument names are repeated
        vnames = {}
        for _, (vname, line_no), (_, _) in assigns:
            if vname in vnames:
                print 'error:%s: argument name "%s" already defined' % (line_no, vname)
                sys.exit(1)
            vnames[vname] = None

        # expected argument names
        PARALLEL = 'parallel'
        TILES = 'tiles'
        PERMUT = 'permut'
        UFACTORS = 'unroll_factors'
        SREPLACE = 'scalar_replace'
        VECTOR = 'vectorize'
        PCODE = 'profiling_code'
        CCMD = 'compile_cmd'
        COPTS = 'compile_opts'
        
        # argument information
        parallel = None
        tiles = None
        permut = None
        unroll_factors = None
        scalar_replace = None
        vectorize = None
        profiling_code = None
        compile_cmd = None
        compile_opts = None
        
        # evaluate each argument
        for line_no, (vname, vname_line_no), (rhs, rhs_line_no) in assigns:

            # evaluate RHS expression
            try:
                rhs_val = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the RHS expression' % rhs_line_no
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
                
            # tile sizes
            if vname == TILES:
                if not isinstance(rhs_val, list) and not isinstance(rhs_val, tuple):
                    print 'error:%s: tile size value must be a tuple/list' % rhs_line_no
                    sys.exit(1)
                for t in rhs_val:
                    if not isinstance(t, int) or t <= 0:
                        print ('error:%s: a tile size must be a positive integer. found: "%s"' %
                               (rhs_line_no, t))
                        sys.exit(1)
                tiles = rhs_val
                
            # permutation
            elif vname == PERMUT:
                if not isinstance(rhs_val, list) and not isinstance(rhs_val, tuple):
                    print 'error:%s: permutation value must be a tuple/list' % rhs_line_no
                    sys.exit(1)
                for i in range(0,len(rhs_val)):
                    if i not in rhs_val:
                        print ('error:%s: not a valid permutation order. found: "%s"' %
                               (rhs_line_no, rhs_val))
                        sys.exit(1)
                permut = rhs_val

            # unroll factors
            elif vname == UFACTORS:
                if not isinstance(rhs_val, list) and not isinstance(rhs_val, tuple):
                    print 'error:%s: unroll factors value must be a tuple/list' % rhs_line_no
                    sys.exit(1)
                for t in rhs_val:
                    if not isinstance(t, int) or t <= 0:
                        print ('error:%s: unroll factor must be a positive integer. found: "%s"' %
                               (rhs_line_no, t))
                        sys.exit(1)
                unroll_factors = rhs_val

            # parallelization
            elif vname == PARALLEL:
                if not isinstance(rhs_val, bool):
                    print 'error:%s: parallelize value must be a boolean' % rhs_line_no
                    sys.exit(1)
                parallel = rhs_val

            # scalar replacement
            elif vname == SREPLACE:
                if isinstance(rhs_val, bool):
                    scalar_replace = (rhs_val, 'double')
                elif ((isinstance(rhs_val, tuple) or isinstance(rhs_val, list)) and
                      len(rhs_val) == 2 and isinstance(rhs_val[0], bool) and
                      isinstance(rhs_val[1], str)):
                    scalar_replace = rhs_val
                else:
                    print (('error:%s: scalar replacement value must be in the form of ' +
                            '(True/False, <type-string>)') % rhs_line_no)
                    sys.exit(1)

            # vectorization
            elif vname == VECTOR:
                if not isinstance(rhs_val, bool):
                    print 'error:%s: vectorize value must be a boolean' % rhs_line_no
                    sys.exit(1)
                vectorize = rhs_val

            # profiling code
            elif vname == PCODE:
                if not isinstance(rhs_val, str):
                    print 'error:%s: profiling code value must be a string' % rhs_line_no
                    sys.exit(1)
                profiling_code = rhs_val

            # compile command
            elif vname == CCMD:
                if not isinstance(rhs_val, str):
                    print 'error:%s: compile command value must be a string' % rhs_line_no
                    sys.exit(1)
                compile_cmd = rhs_val
                
            # compile options
            elif vname == COPTS:
                if not isinstance(rhs_val, str):
                    print 'error:%s: compile options value must be a string' % rhs_line_no
                    sys.exit(1)
                compile_opts = rhs_val
                
            # unrecognized argument name
            else:
                print 'error:%s: unrecognized argument name: "%s"' % (vname_line_no, vname)
                sys.exit(1)

        # if an argument is undefined, set it to its default value
        if parallel == None:
            parallel = False
        if tiles == None:
            tiles = []
        if permut == None:
            permut = []
        if unroll_factors == None:
            unroll_factors = []
        if scalar_replace == None:
            scalar_replace = (False, 'double')
        if vectorize == None:
            vectorize = False
        if profiling_code == None:
            print 'error:PolySyn: missing profiling code'
            sys.exit(1)
            
        # generate and return the transformation information
        return TransfInfo(parallel, tiles, permut, unroll_factors, scalar_replace, vectorize,
                          profiling_code, compile_cmd, compile_opts)


