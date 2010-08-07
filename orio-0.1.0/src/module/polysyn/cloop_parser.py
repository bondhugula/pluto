#
# A simple C parser used to extract loop nests that include the hotspot
#

import re, sys

#----------------------------------------------------------------------

class CLoopParser:
    '''A simple C parser for loops extractions'''

    def __init__(self):
        '''To instantiate a C loop parser instance'''
        
        self.hotspot_time_ratio = 0.60
        self.max_num_hotspots = 1

    #-------------------------------------------------------------------

    def __findClosingChar(self, open_char, close_char, code, init_num_open_chars):
        '''Return the position of the matching closing character in the given code'''
        
        num_open_chars = init_num_open_chars
        match_pos = -1
        for i,s in enumerate(code):
            if s == open_char:
                num_open_chars += 1
            elif s == close_char:
                num_open_chars -= 1
                if num_open_chars == 0:
                    match_pos = i
                    break
        return match_pos

    #-------------------------------------------------------------------

    def __parseForLoop(self, linestr, line_num):
        '''Parse the given line string for a for loop'''
        
        m1 = re.match(r'^\s*for\s*\(', linestr)
        if not m1:
            return None
        pos = self.__findClosingChar('(', ')', linestr[m1.end():], 1)
        if pos < 0:
            print 'internal-error:cloop_parser: for-loop is not declared in one line'
            sys.exit(1)
        m2 = re.match(r'\s*\{\s*$', linestr[m1.end():][pos+1:])
        if not m2:
            print 'internal-error:cloop_parser: for-loop does not have an opening curly brace'
            sys.exit(1)
        m3 = re.match(r'^\s*(\w+)\s*=', linestr[m1.end():][:pos])
        if not m3:
            print 'internal-error:cloop_parser: init-expression of for loop cannot be empty'
            sys.exit(1)
        loop_id = m3.group(1)
        return ('forloop', linestr, line_num, loop_id)

    #-------------------------------------------------------------------

    def __parseIfStmt(self, linestr, line_num):
        '''Parse the given line string for an if statement'''
        
        m1 = re.match(r'^\s*if\s*\(', linestr)
        if not m1:
            return None
        pos = self.__findClosingChar('(', ')', linestr[m1.end():], 1)
        if pos < 0:
            print 'internal-error:cloop_parser: if-statement is not declared in one line'
            sys.exit(1)
        m2 = re.match(r'\s*\{\s*$', linestr[m1.end():][pos+1:])
        if not m2:
            print 'internal-error:cloop_parser: if-statement does not have an opening curly brace'
            sys.exit(1)
        return ('ifstmt', linestr, line_num)

    #-------------------------------------------------------------------

    def __parseClosingBrace(self, linestr, line_num):
        '''Parse the given line string for a closing curly brace'''

        m = re.match(r'^\s*\}\s*$', linestr)
        if not m:
            return None
        return ('closebrace', linestr, line_num)

    #-------------------------------------------------------------------

    def __parseMacroStmt(self, linestr, line_num):
        '''Parse the given line string for a statement containing a macro statement Sn(...)'''
        
        m1 = re.match(r'^\s*S\d+\s*\(', linestr)
        if not m1:
            return None
        pos = self.__findClosingChar('(', ')', linestr[m1.end():], 1)
        if pos < 0:
            print 'internal-error:cloop_parser: Sn(...)-statement is not declared in one line'
            sys.exit(1)
        m2 = re.match(r'\s*;\s*$', linestr[m1.end():][pos+1:])
        if not m2:
            print 'internal-error:cloop_parser: Sn(...)-statement does not end with a semicolon'
            sys.exit(1)
        return ('macrostmt', linestr, line_num)
    
    #-------------------------------------------------------------------

    def __parseDirective(self, linestr, line_num):
        '''Parse the given line string for a directive'''
        
        m = re.match(r'^\s*#.*?$', linestr)
        if not m:
            return None
        return ('directive', linestr, line_num)

    #-------------------------------------------------------------------

    def __parseAssignStmt(self, linestr, line_num):
        '''Parse the given line string for an assignment statement'''

        m = re.match(r'^\s*\w+\s*\=.*?;\s*$', linestr)
        if not m:
            return None
        return ('assignment', linestr, line_num)

    #-------------------------------------------------------------------

    def __parseElseStmt(self, linestr, line_num):
        '''Parse the given line string for an else statement'''
        
        m = re.match(r'^\s*\}\s*else\s*\{\s*$', linestr)
        if not m:
            return None
        return ('else', linestr, line_num)

    #-------------------------------------------------------------------

    def __parseWhitespace(self, linestr, line_num):
        '''Parse the given line string for whitespaces'''
        
        m = re.match(r'^\s*$', linestr)
        if not m:
            return None
        return ('space', linestr, line_num)

    #-------------------------------------------------------------------

    def __parseCodeLine(self, linestr, line_num):
        '''Parse the given line string'''
        return (self.__parseForLoop(linestr, line_num) or
                self.__parseIfStmt(linestr, line_num) or
                self.__parseClosingBrace(linestr, line_num) or
                self.__parseMacroStmt(linestr, line_num) or
                self.__parseDirective(linestr, line_num) or
                self.__parseAssignStmt(linestr, line_num) or
                self.__parseElseStmt(linestr, line_num) or
                self.__parseWhitespace(linestr, line_num))

    #-------------------------------------------------------------------

    def __focusToCloogCode(self, pluto_code, hotspots_info):
        '''Narrow focus to the Cloog code only'''

        # find the opening and closing tags of the Cloog code
        open_tag_re = r'/\*\s*Generated\s+from\s+PLuTo.*?\*/'
        close_tag_re = r'/\*\s*End\s+of\s+CLooG\s+code\s*\*/'
        open_m = re.search(open_tag_re, pluto_code)
        close_m = re.search(close_tag_re, pluto_code)
        if (not open_m) or (not close_m):
            print ('error:polysyn: cannot find the opening and closing tags for the Cloog code')
            sys.exit(1)

        # get the Cloog code
        cloog_code = pluto_code[open_m.end():close_m.start()]

        # get the hotspots only related to the Cloog code
        start_line_no = pluto_code[:open_m.start()].count('\n') + 1
        end_line_no = pluto_code[:close_m.start()].count('\n') + 1
        cloog_hotspots_info = []
        for t,l in hotspots_info:
            if start_line_no <= l <= end_line_no:
                l = l - start_line_no + 1
                cloog_hotspots_info.append((t,l))

        # update the Pluto code (with the Cloog code being eliminated and replaced with a new tag)
        tag = '/*@ cloog code @*/'
        pluto_code = pluto_code[:open_m.end()] + tag + pluto_code[close_m.start():]

        # return the Cloog code and Cloog hotspots information, and the Pluto code
        return (cloog_code, cloog_hotspots_info, pluto_code)

    #-------------------------------------------------------------------

    def __insertCloogCode(self, pluto_code, cloog_code):
        '''Insert the Cloog code back into the Pluto code'''

        tag_re = r'/\*@\s*cloog\s+code\s*@\*/'
        pluto_code = re.sub(tag_re, cloog_code, pluto_code)
        return pluto_code

    #-------------------------------------------------------------------

    def __createScopes(self, parsed_clines):
        '''To find the scope of each for-loop and if statement'''

        nparsed_clines = []
        i = 0
        while True:
            if i == len(parsed_clines):
                break
            pcl = parsed_clines[i]

            # for loop
            if pcl[0] == 'forloop':
                head = pcl
                body = []
                tail = None
                num_openings = 1
                while True:
                    i += 1
                    if i == len(parsed_clines):
                        print 'internal-error:cloop_parser: cannot find a matching closing brace'
                        sys.exit(1)
                    pcl = parsed_clines[i]
                    if pcl[0] in ('forloop', 'ifstmt'):
                        num_openings += 1
                        body.append(pcl)
                    elif pcl[0] == 'closebrace':
                        num_openings -= 1
                        if num_openings == 0:
                            tail = pcl
                            break
                        else:
                            body.append(pcl)
                    else:
                        body.append(pcl)
                nparsed_clines.append([head, self.__createScopes(body), tail])
                i += 1

            # if statement
            elif pcl[0] == 'ifstmt':
                head = pcl
                true_stmt = []
                middle = None
                false_stmt = []
                tail = None
                body = true_stmt
                num_openings = 1
                while True:
                    i += 1
                    if i == len(parsed_clines):
                        print 'internal-error:cloop_parser: cannot find a matching closing brace'
                        sys.exit(1)
                    pcl = parsed_clines[i]
                    if pcl[0] in ('forloop', 'ifstmt'):
                        num_openings += 1
                        body.append(pcl)
                    elif pcl[0] == 'else':
                        if num_openings == 1:
                            body = false_stmt
                            middle = pcl
                        else:
                            body.append(pcl)
                    elif pcl[0] == 'closebrace':
                        num_openings -= 1
                        if num_openings == 0:
                            tail = pcl
                            break
                        else:
                            body.append(pcl)
                    else:
                        body.append(pcl)
                nparsed_clines.append([head,
                                       self.__createScopes(true_stmt),
                                       middle,
                                       self.__createScopes(false_stmt),
                                       tail])
                i += 1

            # other code
            else:
                nparsed_clines.append(pcl)
                i += 1

        return nparsed_clines

    #-------------------------------------------------------------------

    def __findInnerLoopNests(self, clines):
        '''To find all inner loop nests'''

        inner_loop_nests = []
        for cl in clines:
            
            if isinstance(cl, list):
                
                if cl[0][0] == 'forloop':
                    head, body, tail = cl
                    _, _, head_l, loop_id = head
                    _, _, tail_l = tail
                    inests = self.__findInnerLoopNests(body)
                    if inests:
                        inests = [([(loop_id, head_l, tail_l, body)] + n) for n in inests]
                    else:
                        inests = [[(loop_id, head_l, tail_l, body)]]
                    inner_loop_nests.extend(inests)

                elif cl[0][0] == 'ifstmt':
                    head, true_stmt, middle, false_stmt, tail = cl
                    true_inests = self.__findInnerLoopNests(true_stmt)
                    false_inests = self.__findInnerLoopNests(false_stmt)
                    inner_loop_nests.extend(true_inests)
                    inner_loop_nests.extend(false_inests)

                else:
                    print 'internal-error:cloop_parser: unrecognized scope: "%s"' % cl
                    sys.exit(1)

        return inner_loop_nests

    #-------------------------------------------------------------------

    def __findHotspotInnerLoopNest(self, clines):
        '''To find the inner loop nest of the hotspot'''

        hotspot_inests = self.__findInnerLoopNests(clines)

        if len(hotspot_inests) > 1:
            print ('internal-error:cloop_parser: hotspot must not contain two (or more) ' +
                   'separate subloop nests')
            sys.exit(1)

        if len(hotspot_inests) == 1:
            return hotspot_inests[0]
        return []

    #-------------------------------------------------------------------

    def __isParallelLoop(self, line_num, clines):
        '''To determine whether the loop specified at the given line number is parallelized'''

        prev_cl = None
        for cl in clines:

            if isinstance(cl, list):
                
                if cl[0][0] == 'forloop':
                    head, body, tail = cl
                    _, _, head_l, _ = head
                    _, _, tail_l = tail
                    if line_num == head_l:
                        return (prev_cl and not isinstance(prev_cl, list) and
                                prev_cl[1].startswith('#pragma omp parallel'))
                    elif head_l < line_num < tail_l:
                        return self.__isParallelLoop(line_num, body)
        
                elif cl[0][0] == 'ifstmt':
                    head, true_stmt, middle, false_stmt, tail = cl
                    _, _, head_l = head
                    middle_l = None
                    if middle:
                        _, _, middle_l = middle
                    _, _, tail_l = tail
                    if middle_l:
                        if head_l < line_num < middle_l:
                            return self.__isParallelLoop(line_num, true_stmt)
                        elif middle_l < line_num < tail_l:
                            return self.__isParallelLoop(line_num, false_stmt)
                    elif head_l < line_num < tail_l:
                        return self.__isParallelLoop(line_num, true_stmt)

                else:
                    print 'internal-error:cloop_parser: unrecognized scope: "%s"' % cl
                    sys.exit(1)

            prev_cl = cl

        print 'internal-error:cloop_parser: code at line number "%s" cannot be found' % line_num
        sys.exit(1)

    #-------------------------------------------------------------------

    def __findOuterLoopNest(self, line_num, clines, outer_loop_nest):
        '''To find the outer loops surrounding the code at the given line number'''

        for cl in clines:
            
            if isinstance(cl, list):

                if cl[0][0] == 'forloop':
                    head, body, tail = cl
                    _, _, head_l, loop_id = head
                    _, _, tail_l = tail
                    n_outer_loop_nest = outer_loop_nest[:] + [(loop_id, head_l, tail_l, body)]
                    if line_num == head_l or line_num == tail_l:
                        return n_outer_loop_nest
                    elif head_l < line_num < tail_l:
                        return self.__findOuterLoopNest(line_num, body, n_outer_loop_nest)
                
                elif cl[0][0] == 'ifstmt':
                    head, true_stmt, middle, false_stmt, tail = cl
                    _, _, head_l = head
                    middle_l = None
                    if middle:
                        _, _, middle_l = middle
                    _, _, tail_l = tail
                    n_outer_loop_nest = outer_loop_nest[:]
                    if line_num == head_l or line_num == tail_l:
                        return n_outer_loop_nest
                    elif middle_l:
                        if line_num == middle_l:
                            print 'internal-error:cloop_parser: hotspot cannot be at else keyword'
                            sys.exit(1)
                        if head_l < line_num < middle_l:
                            return self.__findOuterLoopNest(line_num, true_stmt, n_outer_loop_nest)
                        elif middle_l < line_num < tail_l:
                            return self.__findOuterLoopNest(line_num, false_stmt, n_outer_loop_nest)
                    elif head_l < line_num < tail_l:
                        return self.__findOuterLoopNest(line_num, true_stmt, n_outer_loop_nest)

                else:
                    print 'internal-error:cloop_parser: unrecognized scope: "%s"' % cl
                    sys.exit(1)

            else:
                _, _, l = cl
                if line_num == l:
                    return outer_loop_nest

        print 'internal-error:cloop_parser: code at line number "%s" cannot be found' % line_num
        sys.exit(1)

    #-------------------------------------------------------------------

    def __findHotspotLoopNest(self, line_num, clines):
        '''
        To find the loop nest that includes the hotspot.
        Hotspot restrictions:
        - cannot be a parallelized loop (i.e. preceeded by an OpenMP directive)
        - cannot have a loop that has two (or more) separate subloop nests, such as:
          for (...) {
            for (...) {...}
            ...
            for (...) {...}
          }
        - can be imperfectly nested
        '''

        # get the outer loops that surround the hotspot
        oloop_nest = self.__findOuterLoopNest(line_num, clines, [])

        # remove all parallelized loops
        n_oloop_nest = []
        for lnest in oloop_nest:
            lid, slnum, elnum, lbody = lnest
            if self.__isParallelLoop(slnum, clines):
                n_oloop_nest = []
            else:
                n_oloop_nest.append(lnest)
        oloop_nest = n_oloop_nest

        # remove loops that have two (or more) separate subloop nests
        n_oloop_nest = []
        for lnest in oloop_nest:
            lid, slnum, elnum, lbody = lnest
            iloop_nests = self.__findInnerLoopNests(lbody)
            if len(iloop_nests) <= 1:
                n_oloop_nest.append(lnest)
        oloop_nest = n_oloop_nest

        # extend the loop nest to include more inner loops
        if len(oloop_nest) > 0:
            lid, slnum, elnum, lbody = oloop_nest[-1]
            loop_nest = oloop_nest + self.__findHotspotInnerLoopNest(lbody)
        else:
            loop_nest = oloop_nest

        # remove the loop body from the elements of the loop nest
        loop_nest = [(lid, slnum, elnum) for lid, slnum, elnum, lbody in loop_nest]

        # return the hotspot's loop nest
        return loop_nest

    #-------------------------------------------------------------------
    
    def __markHotspots(self, clines, hotspot_lnests):
        '''To insert header and trailer annotations into each hotspot loop nest'''

        clines = [list(cl) for cl in clines]
        for lnest in hotspot_lnests:
            loop_ids = [l for l,_,_ in lnest]
            _, slnum, elnum = lnest[0]
            clines[slnum-1][1] = ('/*@ hotspot begin %s @*/\n' % loop_ids) + clines[slnum-1][1]
            clines[elnum-1][1] = clines[elnum-1][1] + '\n/*@ hotspot end @*/' 
        return clines

    #-------------------------------------------------------------------
    
    def getHotspotLoopNests(self, pluto_code, hotspots_info):
        '''To get a sequence of hotspot loops'''

        # narrow focus on the Cloog code only
        cloog_code, hotspots_info, pluto_code = self.__focusToCloogCode(pluto_code, hotspots_info)

        # check the hotspots information
        if len(hotspots_info) == 0:
            print 'error:polysyn: cannot find any hotspots in the Cloog code'
            sys.exit(1)

        # prune out some hotspot information
        best_hotspot = hotspots_info[0]
        n_hotspots_info = [best_hotspot]
        max_time_percent, _ = best_hotspot
        for time_percent, line_num in hotspots_info[1:]:
            if time_percent >= self.hotspot_time_ratio * max_time_percent:
                n_hotspots_info.append((time_percent, line_num))
        hotspots_info = n_hotspots_info

        # remove all comments from the Cloog code
        comment_line_re = r'//.*?\n'
        cloog_code = re.sub(comment_line_re, '', cloog_code)
        comment_re = r'/\*.*?\*/'
        while True:
            m = re.search(comment_re, cloog_code)
            if not m:
                break
            num_newlines = m.group().count('\n')
            replacement = '\n' * num_newlines
            cloog_code = cloog_code[:m.start()] + replacement + cloog_code[m.end():]
        
        # parse each line of the Cloog code
        clines = []
        for i, cl in enumerate(cloog_code.split('\n')):
            d = self.__parseCodeLine(cl, i+1)
            if not d:
                print 'internal-error:cloop_parser: unrecognized line of code: "%s"' % cl
                sys.exit(1)
            clines.append(d)

        # find the scope of each for-loop and if statement
        original_clines = clines
        clines = self.__createScopes(clines[:])

        # find all hotspot loop nests
        hotspot_lnests = []
        for time_percent, line_num in hotspots_info:
            lnest = self.__findHotspotLoopNest(line_num, clines)
            if lnest not in hotspot_lnests:
                hotspot_lnests.append(lnest)

        # limit the number of hotspots to be optimized
        hotspot_lnests = hotspot_lnests[:min(len(hotspot_lnests), self.max_num_hotspots)]

        # remove empty hotspot loop nests
        hotspot_lnests = filter(lambda x: len(x) > 0, hotspot_lnests) 

        # check the hotspot loop nests
        if len(hotspot_lnests) == 0:
            print 'error:polysyn: hotspot has no loops to be optimized'
            sys.exit(1)

        # insert header and trailer annotations into each hotspot loop nest
        original_clines = self.__markHotspots(original_clines, hotspot_lnests)

        # update the Cloog code
        original_clines = [cl[1] for cl in original_clines]
        cloog_code = '\n'.join(original_clines)

        # insert the updated Cloog code back into the Pluto code
        pluto_code = self.__insertCloogCode(pluto_code, cloog_code)

        # return the updated Pluto code
        return pluto_code

        



