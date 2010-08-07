#
# A class used to perform syntactic transformation
#

import os, re, shlex, sys

#---------------------------------------------------------

class SynTransformator:
    '''The syntactic transformator'''

    def __init__(self, verbose, permut, unroll_factors, scalar_replace, vectorize):
        '''To instantiate a syntactic transformator instance'''

        self.verbose = verbose
        self.permut = permut
        self.unroll_factors = unroll_factors
        self.scalar_replace = scalar_replace
        self.vectorize = vectorize

    #---------------------------------------------------------

    def __executeAnnTool(self, code):
        '''To execute the annotation tool to perform syntactic code transformations'''
        
        # the used filenames
        ifname = '_orio_unroll_code.i'
        ofname = '_orio_unroll_code.o'
        
        # write the given code into a file
        try:
            f = open(ifname, 'w')
            f.write(code)
            f.close()
        except:
            print 'error: cannot write to file: %s' % ifname
            sys.exit(1)
            
        # execute Orio
        cmd = 'orcc -o %s %s' % (ofname, ifname)
        if self.verbose: print ' running command:\n\t%s\n' % cmd
        try:
            os.system(cmd)
        except:
            print 'error: failed to execute command: "%s"' % cmd
            sys.exit(1)
            
        # read the generated code
        try:
            f = open(ofname, 'r')
            transformed_code = f.read()
            f.close()
        except:
            print 'error: cannot read file: %s' % ofname
            sys.exit(1)
            
        # delete the used files
        for fname in [ifname, ofname]:
            try:
                os.unlink(fname)
            except:
                print 'error: cannot delete files: %s' % fname
                sys.exit(1)
                
        # return the transformed code
        return transformed_code

    #---------------------------------------------------------
    
    def __containsVars(self, code, vars):
        '''To determine if the given code contains the given variable name'''

        for token in shlex.shlex(code):
            for v in vars:
                if token == v:
                    return True
        return False
        
    #---------------------------------------------------------
    
    def __getHotspotPermut(self, loop_nest, permut, code):
        '''To find the right permutation/loop-order for the given hotspot code'''
        
        # if no permutation is needed
        if permut == None or len(permut) == 0:
            return ''

        # update the loop nest and the loop permutation
        if len(permut) > len(loop_nest):
            lpermut = permut[:]
            for i in range(0, len(lpermut)-len(loop_nest)):
                lpermut.remove(i)
            min_i = min(lpermut)
            lpermut = [loop_nest[i-min_i] for i in lpermut]
        else:
            loop_nest = loop_nest[-len(permut):]
            lpermut = [loop_nest[i] for i in permut]

        # get the lower and upper bounds of each loop
        bounds = []
        for lid in lpermut:
            reg_exp = r'for\s*\(\s*%s\s*=\s*(.*?)\s*;\s*%s\s*(<=|<|>|>=)\s*(.*?)\s*;' % (lid, lid)
            m = re.search(reg_exp, code)
            lb = m.group(1)
            ub = m.group(3)
            bounds.append((lb, ub))

        # check the validity of the loop permutation
        valid = True
        for i in range(0, len(lpermut)):
            lb, ub = bounds[i]
            vars = lpermut[i+1:]
            if self.__containsVars(lb, vars) or self.__containsVars(ub, vars):
                valid = False
                break
        
        # if invalid permutation
        if not valid:
            return ''

        # get the permutation parameter code
        permut_code = 'permut = [%s]' % lpermut
        return permut_code

    #---------------------------------------------------------
    
    def __generateAnnotations(self, code):
        '''To generate code-transformation annotations around the hotspot codes'''

        # regular expressions for the hotspot tags
        hotspot_head_tag_re = r'/\*@\s*hotspot\s+begin\s*(.*?)\s*@\*/'
        hotspot_tail_tag_re = r'/\*@\s*hotspot\s+end\s*@\*/'

        # insert code-transformation annotations into the code
        while True:

            # find the hotspot's head tag
            head_m = re.search(hotspot_head_tag_re, code)

            # stop if no hotspot's head tag is found
            if not head_m:
                break
            
            # get the hotspot's loop nest
            lnest = eval(head_m.group(1))

            # find the hotspot's tail tag
            tail_m = re.search(hotspot_tail_tag_re, code)

            # if no hotspot's tail tag cannot be found
            if not tail_m:
                print 'internal-error:polysyn: hotspot tail tag cannot be found'
                sys.exit(1)

            # find the hotspot loop nest code
            hotspot_code = code[head_m.end():tail_m.start()]

            # get the loop permutation of the hotspot
            permut_code = self.__getHotspotPermut(lnest, self.permut, hotspot_code)

            # choose the unroll factors for the hotspot loop nest accordingly
            mlength = min(len(lnest), len(self.unroll_factors))
            if mlength == 0:
                regtile_code = ''
            else:
                lnest = lnest[-mlength:]
                ufactors = self.unroll_factors[:mlength]
                regtile_code = 'regtile = (%s,%s)' % (lnest, ufactors)

            # get the scalar-replacement code
            srep_code = "scalarreplace = (%s, '%s')" % self.scalar_replace

            # get the vectorization code
            vec_code = ''
            if self.vectorize:
                vec_code = "vector = (True, ['ivdep','vector always'])"
                
            # get all performance parameters
            pparams = []
            if permut_code: pparams.append(permut_code)
            if regtile_code: pparams.append(regtile_code)
            if srep_code: pparams.append(srep_code)
            if vec_code: pparams.append(vec_code)

            # generate the code-transformation annotation
            ann_code = ''
            ann_code += '/*@ begin Loop(\n'
            ann_code += 'transform Composite(\n'
            ann_code += ',\n  '.join(pparams)
            ann_code += ')\n'
            ann_code += hotspot_code + '\n'
            ann_code += ') @*/\n'
            ann_code += '/*@ end @*/\n'

            # replace the marked hotspot code with the annotated code
            code = code[:head_m.start()] + ann_code  + code[tail_m.end():]

        # return the annotated code
        return code

    #---------------------------------------------------------
    
    def transform(self, code):
        '''To perform syntactic transformation'''

        # get the annotated code
        annotated_code = self.__generateAnnotations(code)

        # execute the annotation tool to perform the syntactic transformations
        transformed_code = self.__executeAnnTool(annotated_code)

        # return the transformed code
        return transformed_code

