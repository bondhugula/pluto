#
# Implementation of the exhaustive search algorithm 
#

import sys, time
import main.tuner.search.search

#-----------------------------------------------------

class Exhaustive(main.tuner.search.search.Search):
    '''The search engine that uses an exhaustive search approach'''

    def __init__(self, cfrags, axis_names, axis_val_ranges, constraint, time_limit, total_runs,
                 search_opts, cmd_line_opts, ptcodegen, ptdriver, odriver):
        '''To instantiate an exhaustive search engine'''

        main.tuner.search.search.Search.__init__(self, cfrags, axis_names, axis_val_ranges,
                                                 constraint, time_limit, total_runs, search_opts,
                                                 cmd_line_opts, ptcodegen, ptdriver, odriver)

        # read all algorithm-specific arguments
        self.__readAlgoArgs()

        # complain if the total number of search runs is defined (i.e. exhaustive search
        # only needs to be run once)
        if self.total_runs > 1:
            print ('error: the total number of %s search runs must be one (or can be undefined)' %
                   self.__class__.__name__)
            sys.exit(1)
            
    #--------------------------------------------------
        
    def __readAlgoArgs(self):
        '''To read all algorithm-specific arguments'''
                
        for vname, rhs in self.search_opts.iteritems():
            print ('error: unrecognized %s algorithm-specific argument: "%s"' %
                   (self.__class__.__name__, vname))
            sys.exit(1)

    #--------------------------------------------------

    def __getNextCoord(self, coord):
        '''
        Return the next neighboring coordinate to be considered in the search space.
        Return None if all coordinates in the search space have been visited.
        '''
        next_coord = coord[:]
        for i in range(0, self.total_dims):
            ipoint = next_coord[i]
            iuplimit = self.dim_uplimits[i]
            if ipoint < iuplimit-1:
                next_coord[i] += 1
                break
            else:
                next_coord[i] = 0
                if i == self.total_dims - 1:
                    return None
        return next_coord

    #--------------------------------------------------

    def searchBestCoord(self):
        '''
        To explore the search space and retun the coordinate that yields the best performance
        (i.e. minimum performance cost).
        '''

        if self.verbose: print '\n----- begin exhaustive search -----'
        
        # record the best coordinate and its best performance cost
        best_coord = None
        best_perf_cost = self.MAXFLOAT
        
        # start the timer
        start_time = time.time()

        # start from the origin coordinate (i.e. [0,0,...])
        coord = [0] * self.total_dims 

        # evaluate every coordinate in the search space
        while True:

            # determine the performance cost of the current coordinate
            perf_cost = self.getPerfCost(coord)

            if self.verbose: print 'coordinate: %s, cost: %s' % (coord, perf_cost)

            # compare to the best result so far
            if perf_cost < best_perf_cost:
                best_coord = coord
                best_perf_cost = perf_cost
                if self.verbose: print '>>>> best coordinate found: %s, cost: %s' % (coord,perf_cost)

            # check if the time is up
            if self.time_limit > 0 and (time.time()-start_time) > self.time_limit:
                break

            # move to the next coordinate in the search space
            coord = self.__getNextCoord(coord)

            # check if all coordinates have been visited
            if coord == None:
                break

        # compute the total search time
        search_time = time.time() - start_time

        if self.verbose: print '----- end exhaustive search -----'
        if self.verbose: print '----- begin summary -----'
        if self.verbose: print ' best coordinate: %s, cost: %s' % (best_coord, best_perf_cost)
        if self.verbose: print ' total search time: %.2f seconds' % search_time
        if self.verbose: print '----- end summary -----'
        
        # return the best coordinate
        return best_coord
            

