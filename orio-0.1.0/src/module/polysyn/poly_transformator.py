#
# A class used to perform polyhedral-based transformation by using the available tool,
# called Pluto. The used polyhedral transformation is loop tiling and automatic parallelization
# for multicore platforms.
#

import glob, re, os, sys

#---------------------------------------------------------

class PolyTransformator:
    '''The polyhedral transformator'''

    def __init__(self, verbose, parallel, tiles):
        '''To instantiate a polyhedral transformator instance'''

        self.verbose = verbose
        self.tiles = tiles
        self.parallel = parallel
        
    #---------------------------------------------------------
    
    def __plutoTransform(self, code):
        '''Use Pluto to perform polyhedral transformations'''

        # check if Pluto has been correctly installed
        if os.popen('polycc').read() == '':
            print 'error: Pluto is not installed. Cannot use "polycc" command.'
            sys.exit(1)

        # check loop tiling
        use_tiling = True
        if len(self.tiles) == 0:
            use_tiling = False

        # write the tile sizes into "tile.sizes" file
        if use_tiling:
            ts_fname = 'tile.sizes'
            content = ''
            for t in self.tiles:
                content += '%s\n' % t
            try:
                f = open(ts_fname, 'w')
                f.write(content)
                f.close()
            except:
                print 'error: cannot write to file: %s' % ts_fname
                sys.exit(1)
                
        # write the annotation body code into a file
        fname = '_orio_polysyn.c'
        try:
            f = open(fname, 'w')
            f.write(code)
            f.close()
        except:
            print 'error: cannot open file for writing: %s' % fname
            sys.exit(1)

        # create the Pluto command
        cmd = 'polycc %s --noprevector' % fname
        if self.parallel:
            cmd += ' --parallel'
        if use_tiling:
            cmd += ' --tile --l2tile'

        # execute Pluto
        if self.verbose: print ' running command:\n\t%s\n' % cmd 
        try:
            os.system(cmd)
        except:
            print 'error: failed to run command: %s' % cmd
            sys.exit(1)
   
        # delete unneeded files
        path_name, ext = os.path.splitext(fname)
        removed_fnames = [fname] + glob.glob(path_name + '.kernel.*')
        if use_tiling:
            removed_fnames += [ts_fname]
        for f in removed_fnames:
            try:
                os.unlink(f)
            except:
                print 'error: failed to remove file: %s' % f
                sys.exit(1)

        # get the Pluto-generated code
        plutogen_fnames = glob.glob(path_name + '.*' + ext)
        if len(plutogen_fnames) != 1:
            print 'error: failed to generate Pluto-transformed code'
            sys.exit(1)
        plutogen_fname = plutogen_fnames[0]
        try:
            f = open(plutogen_fname, 'r')
            pluto_code = f.read()
            f.close()
        except:
            print 'error: cannot open file for writing: %s' % fname
            sys.exit(1)
            
        # delete the Pluto-generated file
        try:
            os.unlink(plutogen_fname)
        except:
            print 'error: failed to remove file: %s' % plutogen_fname
            sys.exit(1)

        # return the Pluto-generated code
        return pluto_code

    #---------------------------------------------------------
    
    def transform(self, code):
        '''To perform loop tiling and parallelization using Pluto'''
        
        # use Pluto to perform polyhedral transformations
        pluto_code = self.__plutoTransform(code)
        
        # return the Pluto-generated code
        return pluto_code
    
    
