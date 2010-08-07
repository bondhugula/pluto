#
# The main file of the Orio tool
#

import os, sys

#----------------------------------------------

# source language types
C_CPP = 1
FORTRAN = 2

#----------------------------------------------

def start(argv, lang):
    '''The main starting procedure'''

    # check for Fortran source, which is not supported yet now
    if lang == FORTRAN:
        print 'error: no support for Fortran yet'
        sys.exit(1)

    # include the annotation tool in the Python's search path
    abs_path = os.path.abspath(__file__)
    cur_dir = os.path.dirname(abs_path)
    base_dir,_ = os.path.split(cur_dir)
    sys.path.append(base_dir)
    
    # import other required Python packages
    import ann_parser, cmd_line_opts, opt_driver, tspec.tspec

    # get the command line
    cline_opts = cmd_line_opts.CmdParser().parse(argv)

    # need to be verbose?
    verbose = cline_opts.verbose

    if verbose: print '\n====== START ORIO ======'

    # read source code
    if verbose: print '\n----- begin reading the source file: %s -----' % cline_opts.src_filename
    try:
        f = open(cline_opts.src_filename, 'r')
        src_code = f.read()
        f.close()
    except:
        print 'error: cannot open file for reading: %s' % cline_opts.src_filename
        sys.exit(1)
    if verbose: print '----- finish reading the source file -----'

    # obtain the mapping for performance tuning specifications
    specs_map = {}
    if cline_opts.spec_filename:
        if verbose: print ('\n----- begin reading the tuning specification file: %s -----' %
                           cline_opts.spec_filename)
        try:
            f = open(cline_opts.spec_filename, 'r')
            tspec_prog = f.read()
            f.close()
        except:
            print 'error: cannot open file for reading: %s' % cline_opts.spec_filename
            sys.exit(1)
        specs_map = tspec.tspec.TSpec().parseProgram(tspec_prog)
        if verbose: print '----- finish reading the tuning specification -----'

    # parse the source code and return a sequence of code fragments
    if verbose: print '\n----- begin parsing annotations -----'
    cfrags = ann_parser.AnnParser().parse(src_code)
    if verbose: print '----- finish parsing annotations -----'

    # perform optimizations based on information specified in the annotations
    if verbose: print '\n----- begin optimizations -----'
    odriver = opt_driver.OptDriver(specs_map, cline_opts)
    optimized_code_seq = odriver.optimizeCodeFrags(cfrags, {}, True)
    if verbose: print '----- finish optimizations -----'

    # remove all annotations from output
    if cline_opts.erase_annot:
        if verbose: print '\n----- begin removing annotations from output-----'
        optimized_code_seq = [[ann_parser.AnnParser().removeAnns(c), i] \
                              for c, i in optimized_code_seq]
        if verbose: print '----- finish removing annotations from output-----'

    # write output
    if verbose: print '\n----- begin writing the output file(s) -----'
    for optimized_code, input_params in optimized_code_seq:
        out_filename = cline_opts.out_filename
        if len(optimized_code_seq) > 1:
            path_name, ext = os.path.splitext(cline_opts.out_filename)
            suffix = ''
            for pname, pval in input_params:
                suffix += '_%s_%s' % (pname, pval)
            out_filename = ('%s%s' % (path_name, suffix)) + ext
        if verbose: print '--> writing output to: %s' % out_filename
        try:
            f = open(out_filename, 'w')
            f.write(optimized_code)
            f.close()
        except:
            print 'error: cannot open file for writing: %s' % cline_opts.out_filename
            sys.exit(1)
    if verbose: print '----- finish writing the output file(s) -----'

    if verbose: print '\n====== END ORIO ======'
