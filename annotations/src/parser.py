#
# To parse the source code and divide it into smaller code fragments
#

import re, sys
import code_frag, main

#----------------------------------------

# begin and end keywords
__begin_kword = r'BEGIN'
__end_kword = r'END'
__begin_re = __begin_kword.lower()
__end_re = __end_kword.lower()

#----------------------------------------
# For C/C++ source 
#----------------------------------------

__open_delimC = r'/\*@'
__close_delimC = r'@\*/'
__whitespaceC = r'\s'
__any_charC = r'(.|\n)'

__annot_patternC = (__open_delimC + __any_charC + r'*?' + __close_delimC)
__leader_annotC = (__open_delimC + __whitespaceC + r'*' + __begin_re +
                   __whitespaceC + r'+(\w+)' + __whitespaceC + r'*\((' +
                   __any_charC + r'*?)\)' + __whitespaceC + r'*' + __close_delimC)
__trailer_annotC = (__open_delimC + __whitespaceC + r'*' + __end_re +
                    __whitespaceC + r'*' + __close_delimC)

#----------------------------------------

def __isLeaderAnnotC(s):
    '''Check if the given string is a leader annotation'''
    return re.match('^' + __leader_annotC + '$', s)

def __isTrailerAnnotC(s):
    '''Check if the given string is a trailer annotation'''
    return re.match('^' + __trailer_annotC + '$', s)

def __isAnnotC(s):
    '''Check if the given string is an annotation'''
    return re.match('^' + __annot_patternC + '$', s)

def __getAnnotPatternC():
    '''Check if the given string has an annotation pattern'''
    return __annot_patternC

def __getModuleNameC(s):
    '''Return the annotation module name'''
    m = re.match('^' + __leader_annotC + '$', s)
    return m.group(1)

def __getModuleBodyC(s):
    '''Return the module body block'''
    m = re.match('^' + __leader_annotC + '$', s)
    return m.group(2)

#----------------------------------------
# For FORTRAN source
#----------------------------------------

__comment_headerF = r'\n[cC\*!]'
__open_delimF = __comment_headerF + r'\[--'
__close_delimF = r'--\]'
__whitespaceF = r'([ \t\r\f\v]|' + __comment_headerF + ')' 
__any_charF = r'(.|' + __comment_headerF + ')'

__annot_patternF = (__open_delimF + __any_charF + '*?' + __close_delimF + r'.*')
__leader_annotF = (__open_delimF + __whitespaceF + r'*' + __begin_re +
                   __whitespaceF + r'+(\w+)' + __whitespaceF + r'*\((' +
                   __any_charF + '*?)\)' + __whitespaceF + r'*' + __close_delimF + r'.*')
__trailer_annotF = (__open_delimF + __whitespaceF + r'*' + __end_re +
                    __whitespaceF + r'*' + __close_delimF + r'.*')

#----------------------------------------

def __isLeaderAnnotF(s):
    '''Check if the given string is a leader annotation (Fortran)'''
    return re.match('^' + __leader_annotF + '$', s)

def __isTrailerAnnotF(s):
    '''Check if the given string is a trailer annotation (Fortran)'''
    return re.match('^' + __trailer_annotF + '$', s)

def __isAnnotF(s):
    '''Check if the given string is an annotation (Fortran)'''
    return re.match('^' + __annot_patternF + '$', s)

def __getAnnotPatternF():
    '''Check if the given string has an annotation pattern (Fortran)'''
    return __annot_patternF

def __getModuleNameF(s):
    '''Return the annotation module name (Fortran)'''
    m = re.match('^' + __leader_annotF + '$', s)
    return m.group(3)

def __getModuleBodyF(s):
    '''Return the module body block (Fortran)'''
    m = re.match('^' + __leader_annotF + '$', s)
    return m.group(5)

def __getCommentHeaderF():
    '''Return the regular expression of Fortran comment header'''
    return __comment_headerF

#----------------------------------------

def __identifyAnnotRegion(code_frags):
    '''Identify the annotated code regions in the given sequence of code fragments'''

    def __findTrailerAnnot(code_frags):
        '''Return the index position of the trailer annotation that matches with
        the first leader annotation'''
        
        leader_seen_num = 0
        i = 0
        for cf in code_frags:
            if isinstance(cf, code_frag.LeaderAnnotation):
                leader_seen_num += 1    
            elif isinstance(cf, code_frag.TrailerAnnotation):
                leader_seen_num -= 1
                if leader_seen_num == 0:
                    return i
            i += 1
        return -1

    def __identifyAnnotRegionAcc(code_frags, acc):
        '''Identify the annotated code regions in the given sequence of code fragments
        (Accumulator Passing Style)'''

        # base case
        if len(code_frags) == 0:
            return acc

        # if the first code fragment is a leader annotation
        elif isinstance(code_frags[0], code_frag.LeaderAnnotation):
            i = __findTrailerAnnot(code_frags)
            if i == -1:
                print ('error:%s: no matching trailer annotation exists' %
                       code_frags[0].line_no)
                sys.exit(1)
            annot_body = __identifyAnnotRegion(code_frags[1:i])
            acc.append(code_frag.AnnotatedCodeRegion(code_frags[0], annot_body,
                                                     code_frags[i]))
            return __identifyAnnotRegionAcc(code_frags[i+1:], acc)

        # if the first code fragment is a trailer annotation
        elif isinstance(code_frags[0], code_frag.TrailerAnnotation):
            print 'error:%s: no matching leader annotation exists' % code_frags[0].line_no
            sys.exit(1)
        
        # if the first code fragment is a non annotation
        elif isinstance(code_frags[0], code_frag.NonAnnotation):
            acc.append(code_frags[0])
            return __identifyAnnotRegionAcc(code_frags[1:], acc)

        # unknown type of code fragment
        else:
            print 'internal error: unexpected type of code fragment'
            sys.exit(1)

    return __identifyAnnotRegionAcc(code_frags, [])

#-----------------------------------------

def __createCodeFragment(code, prev_code_frag, lang):
    '''Analyze the given code and then create the right code fragment'''

    # compute the indentation and line number
    if prev_code_frag:
        indent = ''
        if lang == main.C_CPP:
            pos = prev_code_frag.code.rfind('\n')
            if pos >= 0:
                indent = ' ' * (len(prev_code_frag.code) - pos - 1)
        line_no = prev_code_frag.line_no + prev_code_frag.code.count('\n')
    else:
        indent = ''
        line_no = 0

    # generate the right code fragment
    if lang == main.C_CPP:

        # leader annotation
        if __isLeaderAnnotC(code):
            module_name = __getModuleNameC(code)
            module_body = __getModuleBodyC(code)
            return code_frag.LeaderAnnotation(code, indent, line_no,
                                              module_name, module_body)

        # trailer annotation
        elif __isTrailerAnnotC(code):
            return code_frag.TrailerAnnotation(code, indent, line_no)            

        # unknown annotation structure
        elif __isAnnotC(code):
            print 'error:%s: unknown form of annotation' % line_no
            sys.exit(1)

        # non annotation
        else:
            return code_frag.NonAnnotation(code, indent, line_no)

    elif lang == main.FORTRAN:

        # leader annotation
        if (__isLeaderAnnotF(code)):
            module_name = __getModuleNameF(code)
            module_body = re.sub(__getCommentHeaderF(), '\n',
                                 __getModuleBodyF(code))
            return code_frag.LeaderAnnotation(code, indent, line_no,
                                              module_name, module_body)

        # trailer annotation
        elif __isTrailerAnnotF(code):
            return code_frag.TrailerAnnotation(code, indent, line_no)

        # unknown annotation structure
        elif __isAnnotF(code):
            print 'error:%s: unknown form of annotation' % line_no
            sys.exit(1)

        # non annotation
        else:
            return code_frag.NonAnnotation(code, indent, line_no)

    else:
        print 'internal error: unknown source language'
        sys.exit(1)        
                
#-----------------------------------------

def parse(src_code, lang):
    '''Divide the source code into smaller code fragments'''

    # check for an empty source code
    if not src_code:
        return []

    # append newlines at the beginning and end of the source file
    src_code = '\n' + src_code + '\n'    

    # get all annotation matches
    annot_matches = []
    if lang == main.C_CPP:
        annot_matches = re.finditer(__getAnnotPatternC(), src_code)
    elif lang == main.FORTRAN:
        annot_matches = re.finditer(__getAnnotPatternF(), src_code)
    else:
        print 'internal error: unknown source language'
        sys.exit(1)

    # create a sequence of code fragments from the source code
    code_frags = []
    prev_code_frag = None
    prev_end_pos = 0
    for m in annot_matches:
        start_pos, end_pos = m.span()
        assert(prev_end_pos <= start_pos), '"prev_end_pos < start_pos" must hold'
        if prev_end_pos < start_pos:
            cf = __createCodeFragment(src_code[prev_end_pos:start_pos],
                                      prev_code_frag, lang)
            code_frags.append(cf)
            prev_code_frag = cf
        cf = __createCodeFragment(src_code[start_pos:end_pos], prev_code_frag, lang)
        code_frags.append(cf)
        prev_code_frag = cf
        prev_end_pos = end_pos
    if prev_end_pos <= len(src_code) - 1:
        code_frags.append(__createCodeFragment(src_code[prev_end_pos:],
                                               prev_code_frag, lang))

    # group the corresponding code fragments into a single annotation region
    code_frags = __identifyAnnotRegion(code_frags)

    return code_frags
