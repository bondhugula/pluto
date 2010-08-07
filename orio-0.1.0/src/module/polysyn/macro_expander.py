#
# A class used to expand statements that use macro definitions
# In the Pluto-generated code, all statements are represented as Sn(...), which
# refers to some macro definitions. To make syntactic transformations easier, so we need to
# replace all those statements with the corresponding original statements.
#

import os, re, sys

#---------------------------------------------------------

class MacroExpander:
    '''The macro expander used to substitute to the original statements'''

    def __init__(self):
        '''To instantiate a macro expander instance'''
        pass

    #---------------------------------------------------------
    
    def __extractMacroDefs(self, code):
        '''To extract all macro definitions in the form of Sn(...)'''

        # regular expression
        macro_def_re = r'#\s*define\s+S\d+\s*\((.|\n)*?\).*\n'

        # find all macro references in the code
        match_iter = re.finditer(macro_def_re, code)

        # get the macro definitions
        macro_defs = []
        macro_def_positions = []
        for match_obj in match_iter:
            macro_defs.append(match_obj.group(0))
            macro_def_positions.append((match_obj.start(), match_obj.end()))

        # return the macro definitions code
        return macro_defs, macro_def_positions

    #---------------------------------------------------------

    def __rewriteMacroDefs(self, macro_defs):
        '''To rewrite "{<statement>;}" to "<statement>" '''

        rexp = r'\{([^;]*?);\}'
        n_macro_defs = []
        for md in macro_defs:
            match_obj = re.search(rexp, md)
            if match_obj:
                stmt_code = match_obj.group(1)
                md = re.sub(rexp, stmt_code, md)
                n_macro_defs.append(md)
        return n_macro_defs
            
    #---------------------------------------------------------
    
    def __extractMacroRefs(self, code):
        '''To extract all macro references in the form of Sn(...)'''

        # regular expressions
        macro_ref_header_re = r'[^A-Za-z0-9_](S\d+\s*\()'

        # get the macro references
        macro_refs = []
        macro_ref_positions = []
        last_pos = 0
        while True:

            # find the next macro reference
            match_obj = re.search(macro_ref_header_re, code)

            # if no match
            if not match_obj:
                break

            # find the next matching closing parenthesis
            num_open_parenth = 1
            match_pos = None
            sub_code = code[match_obj.end(1):]
            for i,s in enumerate(sub_code):
                if s == '(':
                    num_open_parenth += 1
                elif s == ')':
                    num_open_parenth -= 1
                    if num_open_parenth == 0:
                        match_pos = i
                        break
            if match_pos == None:
                print 'error: no matching ")" exists for "%s"' % match_obj.group(1)
                sys.exit(1)

            # get the macro reference and its position
            ref = match_obj.group(1) + sub_code[:match_pos+1]
            start_pos = last_pos + match_obj.start(1)
            end_pos = last_pos + match_obj.end(1) + match_pos + 1
            macro_refs.append(ref)
            macro_ref_positions.append((start_pos, end_pos))

            # update the code and the last end position
            code = code[(end_pos-last_pos):]
            last_pos = end_pos

        # return the macro references
        return macro_refs, macro_ref_positions

   #---------------------------------------------------------
    
    def __substituteMacro(self, macro_defs_code, macro_ref_code):
        '''
        To substitute the given macro reference based on the given macro definitions.
        We perform such substitution by executing 'gcc -E'.
        '''

        # regular expressions
        directive_re = r'#.*?\n'

        # generate the source code
        code = ''
        code += macro_defs_code
        code += '\n\n'
        code += macro_ref_code
        code += '\n'

        # set the used file names
        src_fname = '_orio_macro_expand.c'
        out_fname = '_orio_macro_expand.out'
        
        # write the code into a C file
        try:
            f = open(src_fname, 'w')
            f.write(code)
            f.close()
        except:
            print 'error: failed to write code into file: %s' % src_fname
            sys.exit(1)

        # execute 'gcc -E' to perform the substitution
        cmd  = 'gcc -E %s -o %s' % (src_fname, out_fname)
        try:
            os.system(cmd)
        except Exception, e:
            print 'error: failed to execute the testing code: "%s"' % cmd
            print ' --> %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)

        # read the output file
        try:
            f = open(out_fname, 'r')
            out_code = f.read()
            f.close()
        except:
            print 'error: failed to read from file: %s' % out_fname
            sys.exit(1)

        # delete all the generated files
        try:
            os.unlink(src_fname)
            os.unlink(out_fname)
        except:
            print 'error: failed to delete files: %s and %s' % (src_fname, out_fname)
            sys.exit(1)

        # get the substitution code
        subs_code = re.sub(directive_re, '', out_code).strip()

        # return the substitution code
        return subs_code

    #---------------------------------------------------------

    def replaceStatements(self, pluto_code):
        '''
        To take the Pluto-generated code as input, and expand all macro references in loop bodies.
        Macros to be expanded are always in the form Sn(...).
        '''

        # extract all relevant macro definitions
        macro_defs, macro_def_positions = self.__extractMacroDefs(pluto_code)

        # rewrite "{<statement>;}" to "<statement>"
        macro_defs = self.__rewriteMacroDefs(macro_defs)

        # delete all relevant macro definitions
        rev_macro_def_positions = macro_def_positions[:]
        rev_macro_def_positions.reverse()
        for start_pos, end_pos in rev_macro_def_positions:
            pluto_code = pluto_code[:start_pos] + pluto_code[end_pos:]

        # extract all relevant macro references
        macro_refs, macro_ref_positions = self.__extractMacroRefs(pluto_code)

        # find a substitution for each macro reference
        macro_defs_code = reduce(lambda x,y: x+y, macro_defs, '')
        macro_subs = [self.__substituteMacro(macro_defs_code, r) for r in macro_refs]
        
        # replace all macro references in the code with the substitutions
        rev_macro_ref_positions = macro_ref_positions[:]
        rev_macro_subs = macro_subs[:]
        rev_macro_ref_positions.reverse()
        rev_macro_subs.reverse()
        for sub, (start_pos, end_pos) in zip(rev_macro_subs, rev_macro_ref_positions):
            pluto_code = pluto_code[:start_pos] + sub + pluto_code[end_pos:]

        # return the Pluto-generated code with substitutions applied
        return pluto_code
