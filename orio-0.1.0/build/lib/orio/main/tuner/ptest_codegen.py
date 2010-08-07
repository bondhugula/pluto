#
# The code generator for performance-testing code
#

import random, re, sys

#-----------------------------------------------------

class PerfTestCodeGen:
    '''The code generator used to produce a performance-testing code'''

    #-----------------------------------------------------
    
    def __init__(self, input_params, input_decls):
        '''To instantiate the testing code generator'''

        self.input_params = input_params
        self.input_decls = input_decls

        self.iparam_code = self.__genIParams(input_params)
        self.decl_code = self.__genDecls(input_decls)
        self.malloc_code = self.__genMAllocs(input_decls)
        self.init_code = self.__genInits(input_decls)
        
    #-----------------------------------------------------

    def generate(self, code):
        '''Generate the testing code used to get performance cost'''

        raise NotImplementedError('%s: unimplemented abstract function "generate"' %
                                  (self.__class__.__name__))
    
    #-----------------------------------------------------

    def __genIParams(self, input_params):
        '''
        Generate declaration code for:
         - input parameters (using preprocessor definitions)
        '''
        
        # generate definitions for each input parameter
        iparams = ['#define %s %s' % (pname, rhs) for pname, rhs in input_params]
        
        # generate and return the input parameter code
        iparam_code = '\n'.join(iparams)
        return iparam_code

    #-----------------------------------------------------

    def __genDecls(self, input_decls):
        '''
        Generate declaration code for:
         - declarations for the input variables
        '''

        # generate the input variable declarations
        decls = []
        for is_static, vtype, vname, vdims, rhs in input_decls:
            if len(vdims) == 0:
                decls.append('%s %s;' % (vtype, vname))
            else:
                if is_static:
                    dim_code = '[%s]' % ']['.join(vdims)
                    decls.append('%s %s%s;' % (vtype, vname, dim_code))
                else:
                    ptr_code = '*' * len(vdims)
                    decls.append('%s %s%s;' % (vtype, ptr_code, vname))

        # generate and return the declaration code
        decl_code = '\n'.join(decls)
        return decl_code

    #-----------------------------------------------------
    
    def __genMAllocs(self, input_decls):
        '''
        Generate declaration code for:
         - memory allocations for input arrays (that are dynamic arrays)
        '''

        # generate iteration variables
        max_dim = 0
        for _,_,_,vdims,_ in input_decls:
            max_dim = max(max_dim, len(vdims))
        iter_vars = map(lambda x: 'i%s' % x, range(1, max_dim+1))

        # generate code for the declaration of the iteration variables
        if len(iter_vars) == 0:
            ivars_decl_code = ''
        else:
            ivars_decl_code = 'int %s;' % ','.join(iter_vars)
        
        # generate memory allocations for dynamic input arrays
        mallocs = []
        for is_static, vtype, vname, vdims, rhs in input_decls:
            if len(vdims) > 0 and not is_static:
                for i in range(0, len(vdims)):
                    loop_code = ''
                    if i > 0:
                        ivar = iter_vars[i-1]
                        dim = vdims[i-1]
                        loop_code += (' ' * (i-1)) + \
                            'for (%s=0; %s<%s; %s++) {\n' % (ivar, ivar, dim, ivar)
                    dim_code = ''
                    if i > 0:
                        dim_code = '[%s]' % ']['.join(iter_vars[:i])
                    rhs_code = ('(%s%s) malloc((%s) * sizeof(%s%s))' %
                                (vtype, '*' * (len(vdims) - i),
                                 vdims[i], vtype, '*' * (len(vdims) - i - 1)))
                    loop_body_code = (' ' * i) + '%s%s = %s;' % (vname, dim_code, rhs_code)
                    code = loop_code + loop_body_code
                    if code:
                        mallocs.append(code)
                brace_code = '}' * (len(vdims) - 1)
                if brace_code:
                    mallocs.append(brace_code)
        
        # return an empty code if no dynamic memory allocation is needed
        if len(mallocs) == 0:
            return ''

        # generate and return the declaration code
        malloc_code = '\n'.join([ivars_decl_code] + mallocs)
        malloc_code = '  ' + re.sub('\n', '\n  ', malloc_code)
        return malloc_code

    #-----------------------------------------------------
    
    def __genInits(self, input_decls):
        '''
        Generate code for:
         - value initializations for all input variables
        '''

        # generate iteration variables
        max_dim = 0
        for _,_,_,vdims,_ in input_decls:
            max_dim = max(max_dim, len(vdims))
        iter_vars = map(lambda x: 'i%s' % x, range(1, max_dim+1))

        # generate code for the declaration of the iteration variables
        if len(iter_vars) == 0:
            ivars_decl_code = ''
        else:
            ivars_decl_code = 'int %s;' % ','.join(iter_vars)
        
        # generate array value initializations
        inits = []
        for _, vtype, vname, vdims, rhs in input_decls:

            # skip if it does not have an initial value (i.e. RHS == None)
            if rhs == None:
                continue

            # if it is a scalar
            if len(vdims) == 0:
                if rhs == 'random':
                    rhs = '(%s) %s' % (vtype, random.uniform(1, 10))
                inits.append('%s = %s;' % (vname, rhs))
                continue

            # generate array value initialization code
            rank = len(vdims)
            used_iter_vars = iter_vars[:rank]
            loop_code = ''
            for i, (ivar, dim) in enumerate(zip(used_iter_vars, vdims)):
                loop_code += (' ' * i) + 'for (%s=0; %s<%s; %s++)\n' % (ivar, ivar, dim, ivar)
            dim_code = '[%s]' % ']['.join(used_iter_vars)
            if rhs == 'random':
                sum = '(%s)' % '+'.join(used_iter_vars)
                rhs = '%s %% %s + %s' % (sum, 5, 1)
            loop_body_code = (' ' * rank) + '%s%s = %s;' % (vname, dim_code, rhs)
            inits.append(loop_code + loop_body_code)

        # return an empty code if no initialization is needed
        if len(inits) == 0:
            return ''

        # generate and return the initialization code
        init_code = '\n'.join([ivars_decl_code] + inits)
        init_code = '  ' + re.sub('\n', '\n  ', init_code)
        return init_code
    
