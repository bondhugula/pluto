#
# Parser to extract annotations from the source code
#

import re, sys
import code_frag

#----------------------------------------

class AnnParser:
    '''The parser used for annotations extraction'''

    # regular expressions
    __vname_re = r'[A-Za-z_]\w*'
    __any_re = r'(.|\n)'
    __ann_re = r'/\*@' + __any_re + r'*?@\*/'
    __leader_ann_re = (r'/\*@\s*begin\s+(' + __vname_re + r')\s*\(\s*(' + __any_re +
                       r'*?)\s*\)\s*@\*/')
    __trailer_ann_re = r'/\*@\s*end\s*@\*/'
    __non_indent_char_re = r'[^ \t]'
    
    #----------------------------------------

    def __init__(self):
        '''To instantiate the annotation parser'''
        pass
    
    #----------------------------------------

    def __getIndentSizeFrom(self, code):
        '''
        Compute the indentation size based on the given code (i.e. count the number of spaces from
        the end of the given code, until a non-space character or a newline is found)
        '''

        indent_size = 0
        for i in range(len(code)-1, -1, -1):
            if re.match(self.__non_indent_char_re, code[i]):
                break
            else:
                indent_size += 1
        return indent_size

    #----------------------------------------

    def __markAnnCodeRegions(self, code_seq):
        '''Mark all annotation code regions by encapsulating each code region with a list'''

        # initialize the code sequence with annotation code regions
        marked_code_seq = []

        # iterate over all codes
        for i in range(0, len(code_seq)):
            code, code_line_no, indent_size, is_ann = code_seq[i]

            # if an annotation
            if is_ann:

                # if a leader annotation
                if re.match(self.__leader_ann_re, code):

                    # find the index position of a matching trailer annotation
                    trailer_ipos = -1
                    leaders_seen = 1
                    for j in range(i+1, len(code_seq)):
                        t_code, t_code_line_no, t_indent_size, t_is_ann = code_seq[j]
                        if t_is_ann:
                            if re.match(self.__leader_ann_re, t_code):
                                leaders_seen += 1
                            else:
                                leaders_seen -= 1
                                if leaders_seen == 0:
                                    trailer_ipos = j
                                    break

                    # if no matching trailer annotations
                    if trailer_ipos == -1:
                        print 'error:%s: no matching trailer annotation exists' % code_line_no
                        sys.exit(1)

                    # apply recursions on the annotation body and the trailing code sequence
                    body_code_seq = self.__markAnnCodeRegions(code_seq[i+1:trailer_ipos])
                    trailing_code_seq = self.__markAnnCodeRegions(code_seq[trailer_ipos+1:])

                    # return the code sequence with annotation code regions
                    marked_code_seq.append([code_seq[i], body_code_seq, code_seq[trailer_ipos]])
                    return marked_code_seq + trailing_code_seq

                # if a trailer annotation
                else:
                    print 'error:%s: no matching leader annotation exists' % code_line_no
                    sys.exit(1)

            # if a non-annotation
            else:
                marked_code_seq.append(code_seq[i])

        # return the code sequence with annotation code regions
        return marked_code_seq

    #----------------------------------------

    def __getCodeSeq(self, code, line_no):
        '''
        Parse the code and return a code sequence that consists of non-annotations and
        annotation code regions.
        A code region is denoted as a list, which its first element is a leader annotation, and
        its last element is a trailer annotation. And anything in between the two annotations
        is another code sequence.
        '''

        # initialize the code seq
        code_seq = []

        # divide the code into a code sequence of non-annotations and annotations
        while True:

            # find the next annotation in the code
            match_obj = re.search(self.__ann_re, code)

            # if nothing matches
            if not match_obj:
                indent_size = 0
                if len(code_seq) > 0:
                    indent_size = self.__getIndentSizeFrom(code_seq[-1][0])
                if code != '':
                    code_seq.append((code, line_no, indent_size, False))
                break

            # get the leading non-annotation
            non_ann = code[:match_obj.start()]
            non_ann_line_no = line_no
            non_ann_indent_size = 0
            if len(code_seq) > 0:
                non_ann_indent_size = self.__getIndentSizeFrom(code_seq[-1][0])

            # insert the non-annotation into the code sequence
            if non_ann != '':
                code_seq.append((non_ann, non_ann_line_no, non_ann_indent_size, False))

            # get the matching annotation
            ann = code[match_obj.start():match_obj.end()]
            ann_line_no = line_no + code[:match_obj.start()].count('\n')
            ann_indent_size = 0
            if len(code_seq) > 0:
                ann_indent_size = self.__getIndentSizeFrom(code_seq[-1][0])

            # insert the matching annotation into the code sequence
            code_seq.append((ann, ann_line_no, ann_indent_size, True))

            # an unrecognized form of annotation
            if not re.match(self.__leader_ann_re, ann) and not re.match(self.__trailer_ann_re, ann):
                print 'error:%s: unrecognized form of annotation' % ann_line_no
                sys.exit(1)
                    
            # update the code and line number
            line_no += code[:match_obj.end()].count('\n')
            code = code[match_obj.end():]

        # mark all annotation code regions
        code_seq = self.__markAnnCodeRegions(code_seq)

        # return the code sequence
        return code_seq

    #----------------------------------------

    def __getModuleInfo(self, code, line_no):
        '''
        Given the leader annotation code, return the module name, the module code,
        and their corresponding starting line number.
        '''
        
        # parse the given code
        match_obj = re.match(self.__leader_ann_re, code)

        # if not a match
        if not match_obj:
            print 'error:%s: not a leader annotation code' % line_no
            sys.exit(1)

        # create the module info
        mname = match_obj.group(1)
        mname_line_no = line_no + code[:match_obj.start(1)].count('\n')
        mcode = match_obj.group(2)
        mcode_line_no = line_no + code[:match_obj.start(2)].count('\n')
        mod_info = (mname, mname_line_no, mcode, mcode_line_no)

        # return the module info
        return mod_info

    #----------------------------------------

    def __convertToCodeFragment(self, code):
        '''Convert the given code into a code fragment object'''

        # if a code region (indicated by a list)
        if isinstance(code, list):

            # assert that the code list has exactly three elements
            if len(code) != 3:
                print 'internal error: the code list must have a length of three '
                sys.exit(1)

            # get all three elements
            leader, leader_line_no, leader_indent_size, leader_is_ann = code[0]
            body_code_seq = code[1]
            trailer, trailer_line_no, trailer_indent_size, trailer_is_ann = code[2]

            # create a leader-annotation code-fragment
            mname, mname_line_no, mcode, mcode_line_no = self.__getModuleInfo(leader, leader_line_no)
            leader_cfrag = code_frag.LeaderAnn(leader, leader_line_no, leader_indent_size,
                                               mname, mname_line_no, mcode, mcode_line_no)

            # apply recursions on the annotation's body code sequence
            cfrags = map(self.__convertToCodeFragment, body_code_seq)
            
            # create a trailer-annotation code-fragment
            trailer_cfrag = code_frag.TrailerAnn(trailer, trailer_line_no, trailer_indent_size)

            # create a code-region code-fragment
            return code_frag.AnnCodeRegion(leader_cfrag, cfrags, trailer_cfrag)

        # a non-annotation
        else:

            # create a non-annotation code fragment
            code, line_no, indent_size, is_ann = code
            cfrag = code_frag.NonAnn(code, line_no, indent_size)

            # check if the given code is an annotation
            if is_ann:
                print 'internal error:%s: unexpected annotation' % line_no
                sys.exit(1)
            
            # return the code fragment
            return cfrag

    #----------------------------------------

    def removeAnns(self, code):
        '''Remove all annotations from the given code'''
        return re.sub(self.__ann_re, '', code)

    #----------------------------------------

    def parse(self, code, line_no = 1):
        '''Parse the code and return a sequence of code fragments'''

        # parse the code to obtain the code sequence
        code_seq = self.__getCodeSeq(code, line_no)

        # convert the code sequence to a sequence of code fragments
        cfrags = map(self.__convertToCodeFragment, code_seq)

        # return the sequence of code fragments
        return cfrags


