#
# The optimization driver used to initiate the optimization process
#

import sys
import code_frag, dyn_loader, tuner.tuner

#----------------------------------------------------------------

# the name of the performance-tuning annotation
PTUNE_NAME = 'PerfTuning'

# the name of the module containing various code transformations
TMOD_NAME = 'module'

#----------------------------------------------------------------

class OptDriver:
    '''The optimization driver whose function is to initiate the optimization process'''

    def __init__(self, specs_map, cmd_line_opts):
        '''To instantiate an optimization driver'''
        self.specs_map = specs_map
        self.cmd_line_opts = cmd_line_opts
        self.ptuner = tuner.tuner.PerfTuner(specs_map, self)
        self.dloader = dyn_loader.DynLoader()

    #-------------------------------------------------------------

    def __checkPerfTuningAnns(self, cfrags):
        '''
        To ensure that:
        1. Nested performance-tuning annotation does not exist
        2. No more than one performance-tuning annotation exists
        '''

        # iterate over all code fragments
        num_ptune_anns = 0
        for cf in cfrags:

            # if a top-level annotation code region
            if isinstance(cf, code_frag.AnnCodeRegion) and cf.leader_ann.mod_name == PTUNE_NAME:

                # increment the total number of performance-tuning annotations
                num_ptune_anns += 1

                # if the total number of top-level performance-tuning annotations is more than one
                if num_ptune_anns > 1:
                    print (('error:%s: the total number of performance-tuning annotations cannot ' +
                            'be more than one') % cf.leader_ann.mod_name_line_no)
                    sys.exit(1)

                # iterate over all nested code fragments
                nested_cfrags = cf.cfrags[:]
                while len(nested_cfrags) > 0:
                    nested_cf = nested_cfrags.pop(0)

                    # if the nested code fragment is a performance-tuning annotation
                    if (isinstance(nested_cf, code_frag.LeaderAnn) and
                        nested_cf.mod_name == PTUNE_NAME):
                        print (('error:%s: performance-tuning annotations must be defined at ' +
                                'top level and cannot be nested') % nested_cf.line_no)
                        sys.exit(1)
                                                
                    # if the nested code fragment is an annotation code region
                    if isinstance(nested_cf, code_frag.AnnCodeRegion):
                        nested_cfrags.append(nested_cf.leader_ann)
                        nested_cfrags.extend(nested_cf.cfrags)
                        nested_cfrags.append(nested_cf.trailer_ann)
                    
    #-------------------------------------------------------------

    def __optimizeCodeFrag(self, cfrag, perf_params):
        '''Apply optimization described in the annotations to the given code fragment'''

        # apply no optimizations to non-annotation code fragment
        if isinstance(cfrag, code_frag.NonAnn):
            return [(cfrag.code, [])]

        # optimize annotated code region
        elif isinstance(cfrag, code_frag.AnnCodeRegion):

            # initiate empirical performance tuning
            if cfrag.leader_ann.mod_name == PTUNE_NAME:

                # apply empirical performance tuning
                optimized_code_seq = self.ptuner.tune(cfrag.leader_ann.mod_code,
                                                      cfrag.leader_ann.mod_code_line_no,
                                                      cfrag.cfrags)

            # initiate code transformation and generation
            else:

                # recursively apply optimizations to the annotation body
                optimized_body_code_seq = self.optimizeCodeFrags(cfrag.cfrags, perf_params)

                # check the optimized body code sequence
                if len(optimized_body_code_seq) != 1:
                    print ('internal error: the optimized annotation body code cannot ' +
                           'be multiple versions')
                    sys.exit(1)

                # get the optimized body code
                optimized_body_code, _ = optimized_body_code_seq[0]

                # dynamically load the transformation module class
                class_name = cfrag.leader_ann.mod_name
                mod_name = '.'.join([TMOD_NAME, class_name.lower(), class_name.lower()])
                mod_class = self.dloader.loadClass(mod_name, class_name)

                # apply code transformations
                try:
                    transformator = mod_class(perf_params,
                                              cfrag.leader_ann.mod_code,
                                              optimized_body_code,
                                              self.cmd_line_opts,
                                              cfrag.leader_ann.mod_code_line_no,
                                              cfrag.leader_ann.indent_size)
                    optimized_code = transformator.transform()
                except Exception, e:
                    print ('error:%s: encountered an error when transforming annotation "%s"' %
                           (cfrag.leader_ann.mod_name_line_no, cfrag.leader_ann.mod_name))
                    print ' --> %s: %s' % (e.__class__.__name__, e)
                    sys.exit(1)

                # create the optimized code sequence
                optimized_code_seq = [(optimized_code, [])]

            # prepend the leader annotation and append the trailer annotation to each of
            # the optimized code
            leader_code = cfrag.leader_ann.code
            trailer_code = cfrag.trailer_ann.code
            optimized_code_seq = [((leader_code + c + trailer_code), i) \
                                  for c, i in optimized_code_seq]

            # return the optimized code sequence
            return optimized_code_seq
        
        # unexpected type of code fragment
        else:
            print 'internal error: unexpected type of code fragment'
            sys.exit(1)
        
    #-------------------------------------------------------------

    def optimizeCodeFrags(self, cfrags, perf_params, is_top_level = False):
        '''Apply optimizations specified in the annotations to each code fragment'''

        # check for the validity of performance tuning annotations
        if is_top_level:
            self.__checkPerfTuningAnns(cfrags)

        # apply optimizations on each code fragment
        optimized_code_seq = None
        for cf in cfrags:
            cur_seq = self.__optimizeCodeFrag(cf, perf_params)
            if optimized_code_seq == None:
                optimized_code_seq = cur_seq
            else:
                optimized_code_seq = [((c1 + c2), (i1 + i2)) \
                                      for c1, i1 in optimized_code_seq for c2, i2 in cur_seq]

        # return the optimized code
        return optimized_code_seq

