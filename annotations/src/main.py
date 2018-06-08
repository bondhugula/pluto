#
# The main file of the annotation tool
#

import getopt, os, sys
import code_frag, parser

#----------------------------------------------

# source language
C_CPP = 0
FORTRAN = 1
PYTHON = 2

#----------------------------------------------

# global variables
src_filename = None
out_filename = None

#----------------------------------------------

# transformation modules
__loaded_modules = {}
__mod_directory = 'module'

#----------------------------------------------

# usage message
usage_msg = '''
description: compile shell for the annotation software system

usage: %s [options] <ifile> [<ofile>]
  <ifile>          the input file containing the annotated code
  <ofile>          the output file (optional)

options:
  -h, --help       print this message
''' % os.path.basename(sys.argv[0])

#----------------------------------------------

def __parseCommandLineParameters(argv):
    '''Parse command line parameters'''

    global src_filename 
    global out_filename
    
    opts, args = getopt.getopt(argv[1:], 'h', ['help'])
        
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print usage_msg
            sys.exit(1)
                
    if len(args) < 1:
        print 'error: missing arguments'
        print usage_msg
        sys.exit(1)
        
    src_filename = args[0]
        
    if len(args) > 1:
        out_filename = args[1]

    if len(args) > 2:
        print 'error: too many arguments'
        print usage_msg
        sys.exit(1)
        
    if not out_filename:
        dirs, fname = os.path.split(src_filename)
        out_filename = os.path.join(dirs, '_' + fname)

    try:
        f = open(src_filename, 'r')
        f.close()
    except:
        print 'error: cannot open file for reading: %s' % src_filename
        sys.exit(1)
        
    try:
        need_delete = os.path.exists(out_filename)
        f = open(out_filename, 'w')
        f.close()
        if need_delete:
            os.unlink(out_filename)
    except:
        print 'error: cannot open file for writing: %s' % out_filename
        sys.exit(1)
            
#----------------------------------------------

def __transform(c_frag, lang):
    '''Transform the given code fragment'''

    # do not apply any transformation on a non-annotation code fragment
    if isinstance(c_frag, code_frag.NonAnnotation):
        return code_frag.TransformedCode(c_frag.code)

    # transform an annotated code region
    elif isinstance(c_frag, code_frag.AnnotatedCodeRegion):

        # recursively transform the annotation body
        annot_body_code = ''
        for cf in c_frag.annot_body:
            transformed_code_frag = __transform(cf, lang)
            annot_body_code += transformed_code_frag.code

        # load the corresponding transformation module 
        mod_name = c_frag.leader_annot.module_name
        module = None
        if __loaded_modules.has_key(mod_name.lower()):
            module = __loaded_modules[mod_name.lower()]
        else:
            mod_fname = __mod_directory + '.' + mod_name.lower() + '.' + mod_name.lower()
            try:
                module = __import__(mod_fname)
                components = mod_fname.split('.')
                for c in components[1:]:
                    module = getattr(module, c)
            except Exception, e:
                print ('error: module "%s" does not exist (failed to load "%s")' %
                       (mod_name, mod_fname))
                print ' --> cause: %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)
            try:
                getattr(module, mod_name)
            except:
                print 'error: no class "%s" defined in "%s"' %s (mod_name, mod_fname)
                sys.exit(1)
            __loaded_modules[mod_name.lower()] = module

        # load the transformation class from the loaded module
        mod_class = getattr(module, mod_name)

        # call transformation procedure
        leader_annot_info = (c_frag.leader_annot.code,
                             c_frag.leader_annot.indent,
                             c_frag.leader_annot.line_no,
                             c_frag.leader_annot.module_name,
                             c_frag.leader_annot.module_body)
        try:
            transformed_code = mod_class().transform(leader_annot_info,
                                                     annot_body_code,
                                                     c_frag.trailer_annot.code,
                                                     lang)
        except Exception, e:
            print '%s: %s' % (e.__class__.__name__, e)
            sys.exit(1)

        # return the transformed code fragment
        return code_frag.TransformedCode(transformed_code)
                
    # unexpected type of code fragment
    else:
        print 'internal error: unexpected type of code fragment'
        sys,exit(1)

#----------------------------------------------

def main(argv, lang):
    '''The main procedure'''

    # get command line parameters
    __parseCommandLineParameters(argv)

    # check for unsupported source languages
    if lang == FORTRAN or lang == PYTHON:
        print 'error: Fortran and Python are not yet supported'
        sys.exit(1)
    
    # read source code
    try:
        f = open(src_filename, 'r')
        src_code = f.read()
        f.close()
    except:
        print 'error: cannot open file for reading: %s' % src_filename
        sys.exit(1)
    
    # parse the source code and divide it into smaller code fragments
    code_frags = parser.parse(src_code, lang)

    # transform each code fragment
    transformed_code_frags = map(lambda x: __transform(x, lang), code_frags)

    # get the transformed code
    transformed_code = ''
    for cf in transformed_code_frags:
        transformed_code += cf.code

    # write to the output file
    try:
        f = open(out_filename, 'w')
        f.write(transformed_code)
        f.close()
    except:
        print 'error: cannot open file for reading: %s' % src_filename
        sys.exit(1)
    


