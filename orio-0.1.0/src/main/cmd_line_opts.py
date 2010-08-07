#
# Contain a parser to extract the command line options, and a class definition for
# command line options
#

import getopt, os, sys

#----------------------------------------------

# the usage message
USAGE_MSG = '''
description: compile shell for Orio

usage: %s [options] <ifile> 
  <ifile>   input file containing the annotated code

options:
  -e, --erase-annot              remove annotations from the output
  -h, --help                     display this message
  -o <file>, --output=<file>     place the output to <file>
  -s <file>, --spec=<file>       read tuning specifications from <file>
  -v, --verbose                  verbosely show details of the results of the running program
''' % os.path.basename(sys.argv[0])

#----------------------------------------------

class CmdLineOpts:
    '''The command line options'''

    def __init__(self, src_filename, out_filename, spec_filename, verbose, erase_annot):
        '''To instantiate an object that represents the command line options'''

        self.src_filename = src_filename       # the name of the source file
        self.out_filename = out_filename       # the name of the output file
        self.spec_filename = spec_filename     # the name of the tuning specification file
        self.verbose = verbose                 # show details of the results of the running program
        self.erase_annot = erase_annot         # do we need to remove annotations from the output?

#----------------------------------------------

class CmdParser:
    '''Parser for command line options'''
    
    def __init__(self):
        '''To instantiate the command line option parser'''
        pass

    #----------------------------------------------

    def parse(self, argv):
        '''To extract the command line options'''

        # variables to represents the command line options
        src_filename = None      
        out_filename = None      
        spec_filename = None
        verbose = False
        erase_annot = False      

        # get all options
        try:
            opts, args = getopt.getopt(argv[1:],
                                       'eho:s:v',
                                       ['erase-annot', 'help', 'output=', 'spec=', 'verbose'])
        except Exception, e:
            print 'error: %s' % e
            print USAGE_MSG
            sys.exit(1)

        # evaluate all options
        for opt, arg in opts:
            if opt in ('-e', '--erase-annot'):
                erase_annot = True
            elif opt in ('-h', '--help'):
                print USAGE_MSG
                sys.exit(1)
            elif opt in ('-o', '--output'):
                out_filename = arg
            elif opt in ('-s', '--spec'):
                spec_filename = arg
            elif opt in ('-v', '--verbose'):
                verbose = True

        # check on the arguments
        if len(args) < 1:
            print 'error: missing arguments'
            print USAGE_MSG
            sys.exit(1)
        if len(args) > 1:
            print 'error: too many arguments'
            print USAGE_MSG
            sys.exit(1)

        # get the source filename
        src_filename = args[0]

        # create the output filename
        if not out_filename:
            dirs, fname = os.path.split(src_filename)
            out_filename = os.path.join(dirs, '_' + fname)

        # check if the source file is readable
        try:
            f = open(src_filename, 'r')
            f.close()
        except:
            print 'error: cannot open file for reading: %s' % src_filename
            sys.exit(1)

        # check if the tuning specification file is readable
        if spec_filename:
            try:
                f = open(spec_filename, 'r')
                f.close()
            except:
                print 'error: cannot open file for reading: %s' % spec_filename
                sys.exit(1)

        # create an object for the command line options
        cline_opts = CmdLineOpts(src_filename, out_filename, spec_filename, verbose, erase_annot)

        # return the command line option object
        return cline_opts

