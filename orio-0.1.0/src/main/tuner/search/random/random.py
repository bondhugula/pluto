#
# Implementation of the random search algorithm
#

import sys, time
import main.tuner.search.search

#-----------------------------------------------------

class Random(main.tuner.search.search.Search):
    '''
    The search engine that uses a random search approach, enhanced with a local search that finds
    the best neighboring coordinate.

    Below is a list of algorithm-specific arguments used to steer the search algorithm.
      local_distance            the distance number used in the local search to find the best
                                neighboring coordinate located within the specified distance
    '''

    # algorithm-specific argument names
    __LOCAL_DIST = 'local_distance'       # default: 0
    
    #--------------------------------------------------
    
    def __init__(self, cfrags, axis_names, axis_val_ranges, constraint, time_limit, total_runs, 
                 search_opts, cmd_line_opts, ptcodegen, ptdriver, odriver):
        '''To instantiate a random search engine'''

        main.tuner.search.search.Search.__init__(self, cfrags, axis_names, axis_val_ranges,
                                                 constraint, time_limit, total_runs, search_opts,
                                                 cmd_line_opts, ptcodegen, ptdriver, odriver)

        # set all algorithm-specific arguments to their default values
        self.local_distance = 0

        # read all algorithm-specific arguments
        self.__readAlgoArgs()
        
        # complain if both the search time limit and the total number of search runs are undefined
        if self.time_limit <= 0 and self.total_runs <= 0:
            print (('error: %s search requires either (both) the search time limit or (and) the ' +
                    'total number of search runs to be defined') % self.__class__.__name__)
            sys.exit(1)

    #--------------------------------------------------
    
    def __readAlgoArgs(self):
        '''To read all algorithm-specific arguments'''

        # check for algorithm-specific arguments
        for vname, rhs in self.search_opts.iteritems():

            # local search distance
            if vname == self.__LOCAL_DIST:
                if not isinstance(rhs, int) or rhs < 0:
                    print ('error: %s argument "%s" must be a positive integer or zero'
                           % (self.__class__.__name__, vname))
                    sys.exit(1)
                self.local_distance = rhs

            # unrecognized algorithm-specific argument
            else:
                print ('error: unrecognized %s algorithm-specific argument: "%s"' %
                       (self.__class__.__name__, vname))
                sys.exit(1)

    #--------------------------------------------------

    def __initRandomCoord(self, coord_records):
        '''Randomly initialize a coordinate in the search space'''

        # check if all coordinates have been explored
        if len(coord_records) >= self.space_size:
            return None

        # randomly pick a coordinate that has never been explored before
        while True:
            coord = self.getRandomCoord()
            if str(coord) not in coord_records:
                coord_records[str(coord)] = None
                return coord
    
    #--------------------------------------------------
    
    def searchBestCoord(self):
        '''
        To explore the search space and retun the coordinate that yields the best performance
        (i.e. minimum performance cost).
        '''

        if self.verbose: print '\n----- begin random search -----'
        
        # initialize a storage to remember all coordinates that have been explored
        coord_records = {}

        # record the best coordinate and its best performance cost
        best_coord = None
        best_perf_cost = self.MAXFLOAT

        # record the number of runs
        runs = 0
        
        # start the timer
        start_time = time.time()

        # execute the randomized search method
        while True:

            # randomly pick a coordinate in the search space
            coord = self.__initRandomCoord(coord_records)

            # if all coordinates in the search space have been explored
            if coord == None:
                break
            
            # get the performance cost of the current coordinate
            perf_cost = self.getPerfCost(coord)
            old_perf_cost = perf_cost
            
            if self.verbose: print '(run %s) coordinate: %s, cost: %s' % (runs+1, coord, perf_cost)
            
            # perform a local search on the randomly picked coordinate
            coord, perf_cost = self.searchBestNeighbor(coord, self.local_distance)

            # if the neighboring coordinate has a better performance cost
            if perf_cost < old_perf_cost:
                if self.verbose: print ('--> better neighbor found: %s, cost: %s' %
                                        (coord, perf_cost))

            # compare to the best result so far
            if perf_cost < best_perf_cost:
                best_coord = coord
                best_perf_cost = perf_cost
                if self.verbose: print '>>>> best coordinate found: %s, cost: %s' % (coord,perf_cost)
                
            # increment the number of runs
            runs += 1

            # check if the time is up
            if self.time_limit > 0 and (time.time()-start_time) > self.time_limit:
                break

            # check if the maximum limit of runs is reached
            if self.total_runs > 0 and runs >= self.total_runs:
                break

        # compute the total search time
        search_time = time.time() - start_time
        
        if self.verbose: print '----- end random search -----'
        if self.verbose: print '----- begin summary -----'
        if self.verbose: print ' best coordinate: %s, cost: %s' % (best_coord, best_perf_cost)
        if self.verbose: print ' total search time: %.2f seconds' % search_time
        if self.verbose: print ' total completed runs: %s' % runs
        if self.verbose: print '----- end summary -----'
        
        # return the best coordinate
        return best_coord
            
