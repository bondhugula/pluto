#
# Loop transformation submodule that implements a combination of various loop transformations.
#

import sys
import module.loop.submodule.submodule, transformator
import module.loop.submodule.tile.tile
import module.loop.submodule.permut.permut
import module.loop.submodule.regtile.regtile
import module.loop.submodule.scalarreplace.scalarreplace
import module.loop.submodule.boundreplace.boundreplace
import module.loop.submodule.pragma.pragma
import module.loop.submodule.arrcopy.arrcopy

#---------------------------------------------------------------------

class Composite(module.loop.submodule.submodule.SubModule):
    '''The composite loop transformation submodule'''
    
    def __init__(self, perf_params = None, transf_args = None, stmt = None):
        '''To instantiate a composite loop transformation submodule'''
        
        module.loop.submodule.submodule.SubModule.__init__(self, perf_params, transf_args, stmt)

        # transformation submodules
        self.tile_smod = module.loop.submodule.tile.tile.Tile()
        self.perm_smod = module.loop.submodule.permut.permut.Permut()
        self.regt_smod = module.loop.submodule.regtile.regtile.RegTile()
        self.srep_smod = module.loop.submodule.scalarreplace.scalarreplace.ScalarReplace()
        self.brep_smod = module.loop.submodule.boundreplace.boundreplace.BoundReplace()
        self.prag_smod = module.loop.submodule.pragma.pragma.Pragma()
        self.acop_smod = module.loop.submodule.arrcopy.arrcopy.ArrCopy()

    #-----------------------------------------------------------------

    def __readTransfArgs(self, perf_params, transf_args):
        '''Process the given transformation arguments'''

        # all expected argument names
        TILE = 'tile'
        PERMUT = 'permut'
        REGTILE = 'regtile'
        SCALARREP = 'scalarreplace'
        BOUNDREP = 'boundreplace'
        PRAGMA = 'pragma'
        OPENMP = 'openmp'
        VECTOR = 'vector'
        ARRCOPY = 'arrcopy'

        # all expected transformation arguments
        tiles = ([], None)
        permuts = ([], None)
        regtiles = (([],[]), None)
        scalarrep = (False, None)
        boundrep = (False, None)
        pragma = ([], None)
        openmp = ((False, ''), None)
        vector = ((False, ''), None)
        arrcopy = ([], None)

        # iterate over all transformation arguments
        for aname, rhs, line_no in transf_args:
            
            # evaluate the RHS expression
            try:
                rhs = eval(rhs, perf_params)
            except Exception, e:
                print 'error:%s: failed to evaluate the argument expression: %s' % (line_no, rhs)
                print ' --> %s: %s' % (e.__class__.__name__, e)
                sys.exit(1)

            # update transformation arguments
            if aname == TILE:
                tiles = (rhs, line_no)
            elif aname == PERMUT:
                permuts = (rhs, line_no)
            elif aname == REGTILE:
                regtiles = (rhs, line_no)
            elif aname == SCALARREP:
                scalarrep = (rhs, line_no)
            elif aname == BOUNDREP:
                boundrep = (rhs, line_no)
            elif aname == PRAGMA:
                pragma = (rhs, line_no)
            elif aname == OPENMP:
                openmp = (rhs, line_no)
            elif aname == VECTOR:
                vector = (rhs, line_no)
            elif aname == ARRCOPY:
                arrcopy = (rhs, line_no)

            # unknown argument name
            else:
                print 'error:%s: unrecognized transformation argument: "%s"' % (line_no, aname)
                sys.exit(1)

        # check semantics of the transformation arguments
        (tiles, permuts, regtiles, scalarrep, boundrep,
         pragma, openmp, vector, arrcopy) = self.checkTransfArgs(tiles, permuts, regtiles,
                                                                 scalarrep, boundrep, pragma,
                                                                 openmp, vector, arrcopy)

        # return information about the transformation arguments
        return (tiles, permuts, regtiles, scalarrep, boundrep, pragma, openmp, vector, arrcopy)

    #-----------------------------------------------------------------

    def checkTransfArgs(self, tiles, permuts, regtiles, scalarrep, boundrep, pragma,
                        openmp, vector, arrcopy):
        '''Check the semantics of the given transformation arguments'''
        
        # evaluate arguments for loop tiling
        rhs, line_no = tiles
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: tile argument must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        targs = []
        for e in rhs:
            if (not isinstance(e, list) and not isinstance(e, tuple)) or len(e) != 3:
                print (('error:%s: element of tile argument must be in the form of ' +
                        '(<loop-id>,<tsize>,<tindex>): %s') % (line_no, e))
                sys.exit(1)
            loop_id, tsize, tindex = e
            loop_id = self.__convertLoopId(loop_id, line_no)
            tsize, tindex = self.tile_smod.checkTransfArgs((tsize, line_no), (tindex, line_no))
            targs.append((loop_id, tsize, tindex))
        tiles = targs

        # evaluate arguments for loop permutation/interchange
        rhs, line_no = permuts
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: permutation argument must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        for e in rhs:
            seq, = self.perm_smod.checkTransfArgs((e, line_no))
        permuts = rhs

        # evaluate arguments for register tiling
        rhs, line_no = regtiles
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: register-tiling argument must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        if len(rhs) != 2:
            print (('error:%s: register-tiling argument must be in the form of ' +
                    '(<loop-ids>,<ufactors>): %s') % (line_no, rhs))
            sys.exit(1)
        loops, ufactors = rhs
        loops, ufactors = self.regt_smod.checkTransfArgs((loops, line_no), (ufactors, line_no))
        regtiles = (loops, ufactors)
        
        # evaluate arguments for scalar replacement
        rhs, line_no = scalarrep
        if isinstance(rhs, bool) or rhs == 0 or rhs == 1:
            scalarrep = (rhs, None, None)
        else:
            if ((not isinstance(rhs, list) and not isinstance(rhs, tuple)) or len(rhs) < 1 or
                len(rhs) > 3 or (not isinstance(rhs[0], bool) and rhs[0] != 0 and rhs[0] != 1)):
                print (('error:%s: scalar replacement argument must be in the form of ' +
                        '((True|False),<dtype>,<prefix>): %s') % (line_no, rhs))
                sys.exit(1)
            do_scalarrep = rhs[0]
            dtype = None
            prefix = None
            if len(rhs) >= 2:
                dtype = rhs[1]
            if len(rhs) >= 3:
                prefix = rhs[2]
            dtype, prefix = self.srep_smod.checkTransfArgs((dtype, line_no), (prefix, line_no))
            scalarrep = (do_scalarrep, dtype, prefix)

        # evaluate arguments for bound replacement
        rhs, line_no = boundrep
        if isinstance(rhs, bool) or rhs == 0 or rhs == 1:
            boundrep = (rhs, None, None)
        else:
            if ((not isinstance(rhs, list) and not isinstance(rhs, tuple)) or len(rhs) < 1 or
                len(rhs) > 3 or (not isinstance(rhs[0], bool) and rhs[0] != 0 and rhs[0] != 1)):
                print (('error:%s: bound replacement argument must be in the form of ' +
                        '((True|False),<lprefix>,<uprefix>): %s') % (line_no, rhs))
                sys.exit(1)
            do_boundrep = rhs[0]
            lprefix = None
            uprefix = None
            if len(rhs) >= 2:
                lprefix = rhs[1]
            if len(rhs) >= 3:
                uprefix = rhs[2]
            lprefix, uprefix = self.brep_smod.checkTransfArgs((lprefix, line_no), (uprefix, line_no))
            boundrep = (do_boundrep, lprefix, uprefix)

        # evaluate arguments for pragma directives
        rhs, line_no = pragma
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: pragma argument must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        targs = []
        for e in rhs:
            if (not isinstance(e, list) and not isinstance(e, tuple)) or len(e) != 2:
                print (('error:%s: element of pragma directive argument must be in the form of ' +
                        '(<loop-id>,<pragma-strings>): %s') % (line_no, e))
                sys.exit(1)
            loop_id, pragmas = e
            loop_id = self.__convertLoopId(loop_id, line_no)
            pragmas, = self.prag_smod.checkTransfArgs((pragmas, line_no))
            targs.append((loop_id, pragmas))
        pragma = targs

        # evaluate arguments for openmp pragma directive
        rhs, line_no = openmp
        if ((not isinstance(rhs, list) and not isinstance(rhs, tuple)) or len(rhs) != 2 or
            not isinstance(rhs[0], bool)):
            print (('error:%s: element of openmp pragma directive argument must be in the form of ' +
                    '((True|False),<pragma-strings>): %s') % (line_no, rhs))
            sys.exit(1)
        do_openmp, pragmas = rhs
        pragmas, = self.prag_smod.checkTransfArgs((pragmas, line_no))
        openmp = do_openmp, pragmas
        
        # evaluate arguments for vectorization pragma directive
        rhs, line_no = vector
        if ((not isinstance(rhs, list) and not isinstance(rhs, tuple)) or len(rhs) != 2 or
            not isinstance(rhs[0], bool)):
            print (('error:%s: element of vectorization pragma directive argument must be in ' +
                    'the form of ((True|False),<pragma-strings>): %s') % (line_no, rhs))
            sys.exit(1)
        do_vector, pragmas = rhs
        pragmas, = self.prag_smod.checkTransfArgs((pragmas, line_no))
        vector = do_vector, pragmas

        # evaluate arguments for array-copy optimization
        rhs, line_no = arrcopy
        if not isinstance(rhs, list) and not isinstance(rhs, tuple):
            print 'error:%s: array-copy argument must be a list/tuple: %s' % (line_no, rhs)
            sys.exit(1)
        targs = []
        for e in rhs:
            if ((not isinstance(e, list) and not isinstance(e, tuple)) or len(e) > 5 or
                len(e) < 3 or not isinstance(e[0], bool)):
                print (('error:%s: element of tile argument must be in the form of ' +
                        '((True|False),<array-ref-str>,<dim-sizes>,<suffix>,<dtype>): %s') %
                       (line_no, e))
                sys.exit(1)
            dtype = None
            suffix = None
            if len(e) == 3:
                do_acopy, aref, dimsizes = e
            elif len(e) == 4:
                do_acopy, aref, dimsizes, suffix = e
            else:
                do_acopy, aref, dimsizes, suffix, dtype = e
            (aref, suffix,
             dtype, dimsizes)= self.acop_smod.checkTransfArgs((aref, line_no), (suffix, line_no),
                                                              (dtype, line_no), (dimsizes, line_no))
            targs.append((do_acopy, aref, suffix, dtype, dimsizes))
        arrcopy = targs

        # return information about the transformation arguments
        return (tiles, permuts, regtiles, scalarrep, boundrep, pragma, openmp, vector, arrcopy)

    #-----------------------------------------------------------------

    def applyTransf(self, tiles, permuts, regtiles, scalarrep, boundrep,
                    pragma, openmp, vector, arrcopy, stmt):
        '''To apply a sequence of transformations'''

        # perform the composite transformations
        t = transformator.Transformator(tiles, permuts, regtiles, scalarrep,
                                        boundrep, pragma, openmp, vector, arrcopy, self.stmt)
        transformed_stmt = t.transform()

        # return the transformed statement
        return transformed_stmt
        
    #-----------------------------------------------------------------

    def __convertLoopId(self, lid, line_no):
        '''
        Convert the loop ID to a list: [True/False, id1, id2, id3, ...].
        The 'True' boolean value indicates that at least one of the loop ID must exist in the
        statement body. A 'False' value means that it is OK if no loop IDs exist in the statement
        body.
        The sequence of IDs imply that "apply optimizations on id1 (if exist), if not, apply
        optimizations on id2 (if exist), and so on and so forth".
        '''

        # check if the loop ID is well-formed
        if isinstance(lid, str):
            pass
        elif (isinstance(lid, tuple) or isinstance(lid, list)) and len(lid) > 0:
            for i in lid:
                if not isinstance(i, str):
                    print 'error:%s: loop ID must be a string: %s' % (line_no, i)
                    sys.exit(1)
        else:
            print 'error:%s: invalid loop ID representation: %s' % (line_no, lid)
            sys.exit(1)            

        # create the loop ID abstraction
        lids = []
        if isinstance(lid, str):
            lids.append(True)
            lids.append(lid)
        elif (isinstance(lid, tuple) or isinstance(lid, list)) and len(lid) > 0:
            lids.append(isinstance(lid, tuple))
            lids.extend(lid)
        else:
            print 'internal error: incorrect representation of the loop IDs'
            sys.exit(1)
        return lids

    #-----------------------------------------------------------------

    def transform(self):
        '''To apply various loop transformations'''

        # read all transformation arguments
        args_info = self.__readTransfArgs(self.perf_params, self.transf_args)
        (tiles, permuts, regtiles, scalarrep,
         boundrep, pragma, openmp, vector, arrcopy) = args_info

        # perform all transformations
        transformed_stmt = self.applyTransf(tiles, permuts, regtiles, scalarrep, boundrep,
                                            pragma, openmp, vector, arrcopy, self.stmt)

        # return the transformed statement
        return transformed_stmt



