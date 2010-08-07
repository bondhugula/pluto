#
# The tuner class to initiate the empirical performance tuning process
#

import re, sys
import main.dyn_loader, main.tspec.tspec, ptest_codegen_basic, ptest_driver

#--------------------------------------------------

# the name of the module containing various search algorithms
SEARCH_MOD_NAME = 'main.tuner.search'

#--------------------------------------------------

class PerfTuner:
    '''The empirical performance tuner'''

    # regular expressions
    __vname_re = r'[A-Za-z_]\w*'
    __import_re = r'\s*import\s+spec\s+(' + __vname_re + r');\s*'

    #-------------------------------------------------
    
    def __init__(self, specs_map, odriver):
        '''To instantiate an empirical performance tuner object'''

        self.specs_map = specs_map
        self.odriver = odriver
        self.dloader = main.dyn_loader.DynLoader()
        self.cmd_line_opts = odriver.cmd_line_opts
        self.verbose = self.cmd_line_opts.verbose
        
    #-------------------------------------------------

    def __extractTuningInfo(self, code, line_no):
        '''Extract tuning information from the given annotation code'''

        # parse the code
        match_obj = re.match(r'^' + self.__import_re + r'$', code)

        # if the code contains a single import statement
        if match_obj:

            # get the specification name and line number
            spec_name = match_obj.group(1)
            spec_name_line_no = line_no + code[:match_obj.start(1)].count('\n')
            
            # if the tuning info is not defined
            if spec_name not in self.specs_map:
                print 'error:%s: undefined specification: "%s"' % (spec_name_line_no, spec_name)
                sys.exit(1)

            # get the tuning information from the specifications map
            tinfo = self.specs_map[spec_name]

            # return the tuning information
            return tinfo

        # if the tuning specification is hardcoded into the given code
        else:

            # parse the specification code to get the tuning information
            tinfo = main.tspec.tspec.TSpec().parseSpec(code, line_no)

            # return the tuning information
            return tinfo
        
    #-------------------------------------------------

    def __listAllCombinations(self, seqs):
        '''
        Enumerate all combinations of the given sequences.
          e.g. input: [['a','b'],[1,2]] --> [['a',1],['a',2],['b',1],['b',2]]
        '''
        
        # the base case
        if len(seqs) == 0:
            return []
        
        # the recursive step
        trailing_combs = self.__listAllCombinations(seqs[1:])
        if trailing_combs == []:
            trailing_combs = [[]]
        combs = []
        for i in seqs[0]:
            for c in trailing_combs:
                combs.append([i] + c)
                
        # return the combinations
        return combs
    
    #-------------------------------------------------

    def __getProblemSizes(self, iparam_params, iparam_constraints):
        '''Return all valid problem sizes'''

        # combine the input parameter constraints
        iparam_constraint = 'True'
        for vname, rhs in iparam_constraints:
            iparam_constraint += ' and (%s)' % rhs

        # compute all possible combinations of problem sizes
        prob_sizes = []
        pnames, pvalss = zip(*iparam_params)
        for pvals in self.__listAllCombinations(pvalss):
            prob_sizes.append(zip(pnames, pvals))

        # exclude all invalid problem sizes
        n_prob_sizes = []
        for p in prob_sizes:
            try:
                is_valid = eval(iparam_constraint, dict(p))
            except Exception, e:
                print 'error:%s: failed to evaluate the input parameter constraint expression'
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
            if is_valid:
                n_prob_sizes.append(p)
        prob_sizes = n_prob_sizes

        # check if the new problem sizes is empty
        if len(prob_sizes) == 0:
            print ('error: no valid problem sizes exist. please check the input parameter ' +
                   'constraints')
            sys.exit(1)
        
        # return all possible combinations of problem sizes
        return prob_sizes

    #-------------------------------------------------

    def __buildCoordSystem(self, perf_params):
        '''Return information about the coordinate systems that represent the search space'''

        # get the axis names and axis value ranges
        axis_names = []
        axis_val_ranges = []
        for pname, prange in perf_params:
            axis_names.append(pname)
            
            # remove duplications and then perform sorting
            n_prange = []
            for r in prange:
                if r not in n_prange:
                    n_prange.append(r)
            prange = n_prange
            prange.sort()
            axis_val_ranges.append(prange)

        # return the axis names and the axis value ranges
        return (axis_names, axis_val_ranges)
        
    #-------------------------------------------------

    def tune(self, module_body_code, line_no, cfrags):
        '''
        Perform empirical performance tuning on the given annotated code. And return the best
        optimized code variant.
        '''
        
        # extract the tuning information specified from the given annotation
        tinfo = self.__extractTuningInfo(module_body_code, line_no)

        # create a performance-testing code generator for each distinct problem size
        ptcodegens = []
        for prob_size in self.__getProblemSizes(tinfo.iparam_params, tinfo.iparam_constraints):
            c = ptest_codegen_basic.PerfTestCodeGenBasic(prob_size, tinfo.ivar_decls,
                                                         tinfo.ivar_decl_file, tinfo.ivar_init_file,
                                                         tinfo.ptest_skeleton_code_file)
            ptcodegens.append(c)
            
        # create the performance-testing driver
        ptdriver = ptest_driver.PerfTestDriver(tinfo.build_cmd, tinfo.build_opts,
                                               tinfo.pcount_method, tinfo.pcount_reps, 
                                               tinfo.driver_src, tinfo.run_cmd, self.verbose)

        # get the axis names and axis value ranges to represent the search space
        axis_names, axis_val_ranges = self.__buildCoordSystem(tinfo.pparam_params)

        # combine the performance parameter constraints
        pparam_constraint = 'True'
        for vname, rhs in tinfo.pparam_constraints:
            pparam_constraint += ' and (%s)' % rhs

        # dynamically load the search engine class
        class_name = tinfo.search_algo
        mod_name = '.'.join([SEARCH_MOD_NAME, class_name.lower(), class_name.lower()])
        search_class = self.dloader.loadClass(mod_name, class_name)

        # convert the search time limit (from minutes to seconds) and get the total number of
        # search runs
        search_time_limit = 60 * tinfo.search_time_limit
        search_total_runs = tinfo.search_total_runs

        # get the search-algorithm-specific arguments
        search_opts = dict(tinfo.search_opts)

        # perform the performance tuning for each distinct problem size
        optimized_code_seq = []
        for ptcodegen in ptcodegens:

            if self.verbose:
                print '\n----- begin empirical tuning for problem size -----'
                iparams = ptcodegen.input_params[:]
                iparams.sort(lambda x,y: cmp(x[0],y[0]))
                for pname, pvalue in iparams:
                    print ' %s = %s' % (pname, pvalue)

            # create the search engine
            search_eng = search_class(cfrags, axis_names, axis_val_ranges, pparam_constraint,
                                      search_time_limit, search_total_runs, search_opts,
                                      self.cmd_line_opts, ptcodegen, ptdriver, self.odriver)

            # search for the best performance parameters
            best_perf_params, best_perf_cost = search_eng.search()

            # print the best performance parameters
            if self.verbose:
                print '----- the obtained best performance parameters -----'
                pparams = best_perf_params.items()
                pparams.sort(lambda x,y: cmp(x[0],y[0]))
                for pname, pvalue in pparams:
                    print ' %s = %s' % (pname, pvalue)
        
            # generate the optimized code using the obtained best performance parameters
            cur_optimized_code_seq = self.odriver.optimizeCodeFrags(cfrags, best_perf_params)

            # check the optimized code sequence
            if len(cur_optimized_code_seq) != 1:
                print 'internal error: the empirically optimized code cannot be multiple versions'
                sys.exit(1)
            
            # get the optimized code
            optimized_code, _ = cur_optimized_code_seq[0]

            # insert comments into the optimized code to include information about 
            # the best performance parameters and the input problem sizes
            iproblem_code = ''
            iparams = ptcodegen.input_params[:]
            iparams.sort(lambda x,y: cmp(x[0],y[0]))
            for pname, pvalue in iparams:
                if pname == '__builtins__':
                    continue
                iproblem_code += '  %s = %s \n' % (pname, pvalue)
            pparam_code = ''
            pparams = best_perf_params.items()
            pparams.sort(lambda x,y: cmp(x[0],y[0]))
            for pname, pvalue in pparams:
                if pname == '__builtins__':
                    continue
                pparam_code += '  %s = %s \n' % (pname, pvalue)
            info_code = '\n\n/**-- (Generated by Orio) \n'
            info_code += 'Best performance cost: \n'
            info_code += '  %f \n' % best_perf_cost
            info_code += 'Tuned for specific problem sizes: \n'
            info_code += iproblem_code
            info_code += 'Best performance parameters: \n'
            info_code += pparam_code
            info_code += '--**/\n\n'
            optimized_code = info_code + optimized_code

            # store the optimized for this problem size
            optimized_code_seq.append((optimized_code, ptcodegen.input_params[:]))

        # return the optimized code
        return optimized_code_seq


