#
# TSpec (Tuning Specifier) class
# 

import sys
import eval, parser, tune_info

#-----------------------------------------------

class TSpec:
    '''The TSpec (Tuning Specifier)'''

    def __init__(self):
        '''To instantiate a TSpec instance'''
        pass

    #-----------------------------------------------

    def parseSpec(self, spec_code, line_no):
        '''To parse the given specification body code and to return its tuning information'''

        # parse and evaluate the specification statement
        stmt_seq = parser.TSpecParser().parseSpec(spec_code, line_no)
        stmt_seq = eval.TSpecEvaluator().evaluate(stmt_seq)

        # generate tuning information
        tinfo = tune_info.TuningInfoGen().generate(stmt_seq)
        
        # return the tuning information
        return tinfo
        
    #-----------------------------------------------

    def parseProgram(self, prog_code):
        '''To parse the entire tuning specification code and to return its tuning information'''

        # parse and evaluate the entire code
        stmt_seq = parser.TSpecParser().parseProgram(prog_code)
        stmt_seq = eval.TSpecEvaluator().evaluate(stmt_seq)

        # create a generator for performance tuning information
        tinfo_gen = tune_info.TuningInfoGen()

        # generate tuning information for each specification statement and insert it into
        # a specifications mapping
        specs_map = {}
        for s in stmt_seq:
            if s[0] == 'spec':
                _, _, (sname, _), sseq = s
                specs_map[sname] = tinfo_gen.generate(sseq)

        # return the specifications mapping
        return specs_map
        
