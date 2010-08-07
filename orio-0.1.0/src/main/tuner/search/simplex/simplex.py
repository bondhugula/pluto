#
# Implementation of the Nelder-Mead Simplex algorithm
#
# The detailed algorithm is described in the following paper.
#   "Convergence Properties of the Nelder-Mead Simplex Method in Low Dimensions"
#   by Jeffrey C. Lagarias
#

import random, sys, time
import main.tuner.search.search

#-----------------------------------------------------

class Simplex(main.tuner.search.search.Search):
    '''
    The search engine that uses the Nelder-Mead Simplex algorithm, enhanced with a local search
    that finds the best neighboring coordinate.

    Below is a list of algorithm-specific arguments used to steer the search algorithm.
      local_distance            the distance number used in the local search to find the best
                                neighboring coordinate located within the specified distance
      reflection_coef           the reflection coefficient
      expansion_coef            the expansion coefficient
      contraction_coef          the contraction coefficient
      shrinkage_coef            the shrinkage coefficient
    '''

    # algorithm-specific argument names
    __LOCAL_DIST = 'local_distance'       # default: 0
    __REFL_COEF = 'reflection_coef'       # default: 1.0
    __EXP_COEF = 'expansion_coef'         # default: 2.0
    __CONT_COEF = 'contraction_coef'      # default: 0.5
    __SHRI_COEF = 'shrinkage_coef'        # default: 0.5

    #-----------------------------------------------------

    def __init__(self, cfrags, axis_names, axis_val_ranges, constraint, time_limit, total_runs,
                 search_opts, cmd_line_opts, ptcodegen, ptdriver, odriver):
        '''To instantiate a Nelder-Mead simplex search engine'''
        
        main.tuner.search.search.Search.__init__(self, cfrags, axis_names, axis_val_ranges,
                                                 constraint, time_limit, total_runs, search_opts,
                                                 cmd_line_opts, ptcodegen, ptdriver, odriver)

        # other private class variables
        self.__simplex_size = self.total_dims + 1

        # set all algorithm-specific arguments to their default values
        self.local_distance = 0
        self.refl_coefs = [1.0]
        self.exp_coefs = [2.0]
        self.cont_coefs = [0.5]
        self.shri_coef = 0.5

        # read all algorithm-specific arguments
        self.__readAlgoArgs()

        # complain if both the search time limit and the total number of search runs are undefined
        if self.time_limit <= 0 and self.total_runs <= 0:
            print (('error: %s search requires either (both) the search time limit or (and) the ' +
                    'total number of search runs to be defined') % self.__class__.__name__)
            sys.exit(1)

    #-----------------------------------------------------

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

            # reflection coefficient
            elif vname == self.__REFL_COEF:
                if isinstance(rhs, int) or isinstance(rhs, float):
                    rhs = [rhs]
                if isinstance(rhs, list):
                    for n in rhs:
                        if (not isinstance(n, int) and not isinstance(n, float)) or n <= 0:
                            print ('error: %s argument "%s" must be number(s) greater than zero'
                                   % (self.__class__.__name__, vname))
                            sys.exit(1)
                    self.refl_coefs = rhs
                else:
                    print ('error: %s argument "%s" must be number(s) greater than zero'
                           % (self.__class__.__name__, vname))
                    sys.exit(1)
                        
            # expansion coefficient
            elif vname == self.__EXP_COEF:
                if isinstance(rhs, int) or isinstance(rhs, float):
                    rhs = [rhs]
                if isinstance(rhs, list):
                    for n in rhs:
                        if (not isinstance(n, int) and not isinstance(n, float)) or n <= 1:
                            print ('error: %s argument "%s" must be number(s) greater than one'
                                   % (self.__class__.__name__, vname))
                            sys.exit(1)
                    self.exp_coefs = rhs
                else:
                    print ('error: %s argument "%s" must be number(s) greater than one'
                           % (self.__class__.__name__, vname))
                    sys.exit(1)
            
            # contraction coefficient
            elif vname == self.__CONT_COEF:
                if isinstance(rhs, int) or isinstance(rhs, float):
                    rhs = [rhs]
                if isinstance(rhs, list):
                    for n in rhs:
                        if (not isinstance(n, int) and not isinstance(n, float)) or n <= 0 or n >= 1:
                            print ('error: %s argument "%s" must be number(s) between zero and one'
                                   % (self.__class__.__name__, vname))
                            sys.exit(1)
                    self.cont_coefs = rhs
                else:
                    print ('error: %s argument "%s" must be number(s) between zero and one'
                           % (self.__class__.__name__, vname))
                    sys.exit(1)
            
            # shrinkage coefficient
            elif vname == self.__SHRI_COEF:
                if (not isinstance(rhs, int) and not isinstance(rhs, float)) or rhs <= 0 or rhs >= 1:
                    print ('error: %s argument "%s" must be a single number between zero and one'
                           % (self.__class__.__name__, vname))
                    sys.exit(1)
                self.shri_coef = rhs
                
            # unrecognized algorithm-specific argument
            else:
                print ('error: unrecognized %s algorithm-specific argument: "%s"' %
                       (self.__class__.__name__, vname))
                sys.exit(1)
        
    #-----------------------------------------------------

    def __checkSearchSpace(self):
        '''Check the size of the search space if it is valid for this search method'''

        # is the search space size too small?
        # Nelder-Mead requires to initialize a simplex that has N+1 vertices, where N is the
        # number of dimensions
        if self.space_size < self.__simplex_size:
            print (('error: the search space is too small for %s algorithm. ' +
                    'please use the exhaustive search.') % self.__class__.__name__)
            sys.exit(1)

    #-----------------------------------------------------

    def __initRandomSimplex(self, simplex_records):
        '''Randomly initialize a simplex in the search space'''

        # remove some simplex records, if necessary
        max_num_records = 100000
        if len(simplex_records) > max_num_records:
            for i in range(i, int(max_num_records*0.05)):
                simplex_records.popitem()

        # randomly create a new simplex that has never been used before
        while True:

            # randomly pick (N+1) vertices to form a simplex, where N is the number of dimensions
            simplex = []
            while True:
                coord = self.getRandomCoord()
                if coord not in simplex:
                    simplex.append(coord)
                    if len(simplex) == self.__simplex_size:
                        break

            # check if the new simplex has never been used before
            simplex.sort()
            if str(simplex) not in simplex_records:
                simplex_records[str(simplex)] = None
                return simplex
            
    #-----------------------------------------------------

    def __getCentroid(self, coords):
        '''Return a centroid coordinate'''
        total_coords = len(coords)
        centroid = coords[0]
        for c in coords[1:]:
            centroid = self.addCoords(centroid, c)
        centroid = self.mulCoords((1.0/total_coords), centroid)
        return centroid

    def __getReflection(self, coord, centroid):
        '''Return a reflection coordinate'''
        sub_coord = self.subCoords(centroid, coord)
        return map(lambda x: self.addCoords(centroid, self.mulCoords(x, sub_coord)),
                   self.refl_coefs)
    
    def __getExpansion(self, coord, centroid):
        '''Return an expansion coordinate'''
        sub_coord = self.subCoords(coord, centroid)
        return map(lambda x: self.addCoords(centroid, self.mulCoords(x, sub_coord)),
                   self.exp_coefs)
    
    def __getContraction(self, coord, centroid):
        '''Return a contraction coordinate'''
        sub_coord = self.subCoords(coord, centroid)
        return map(lambda x: self.addCoords(centroid, self.mulCoords(x, sub_coord)),
                   self.cont_coefs)

    def __getShrinkage(self, coord, rest_coords):
        '''Return a shrinkage simplex'''
        return map(lambda x: self.addCoords(coord, self.mulCoords(self.shri_coef,
                                                                  self.subCoords(x, coord))),
                   rest_coords)
    
    #-----------------------------------------------------

    def searchBestCoord(self):
        '''To search the coordinate that yields the best performance parameters'''

        if self.verbose: print '\n----- begin simplex search -----'

        # check if the size of the search space is valid for this search
        self.__checkSearchSpace()

        # initialize a storage to remember all initial simplexes that have been explored
        simplex_records = {}

        # record the global best coordinate and its performance cost
        best_global_coord = None
        best_global_perf_cost = self.MAXFLOAT
        
        # record the number of runs
        runs = 0
        
        # start the timer
        start_time = time.time()

        # execute the Nelder-Mead Simplex method
        while True:
            
            # list of the last several moves (used for termination criteria)
            last_simplex_moves = []

            # randomly initialize a simplex in the search space
            simplex = self.__initRandomSimplex(simplex_records)

            if self.verbose: print '\n(run %s) initial simplex: %s' % (runs+1, simplex)

            # get the performance cost of each coordinate in the simplex
            perf_costs = map(self.getPerfCost, simplex)

            while True:

                # sort the simplex coordinates in an increasing order of performance costs
                sorted_simplex_cost = zip(simplex, perf_costs)
                sorted_simplex_cost.sort(lambda x,y: cmp(x[1],y[1]))

                # unbox the coordinate-cost tuples
                simplex, perf_costs = zip(*sorted_simplex_cost)
                simplex = list(simplex)
                perf_costs = list(perf_costs)
                
                if self.verbose: print '-> simplex: %s' % simplex

                # check if the time is up
                if self.time_limit > 0 and (time.time()-start_time) > self.time_limit:
                    break
                
                # termination criteria: a loop is present
                if str(simplex) in last_simplex_moves:
                    if self.verbose: print '-> converged with simplex: %s' % simplex
                    break

                # record the last several simplex moves (used for the termination criteria)
                last_simplex_moves.append(str(simplex))
                while len(last_simplex_moves) > 10:
                    last_simplex_moves.pop(0)
                
                # best coordinate
                best_coord = simplex[0]
                best_perf_cost = perf_costs[0]

                # worst coordinate
                worst_coord = simplex[len(simplex)-1]
                worst_perf_cost = perf_costs[len(perf_costs)-1]

                # 2nd worst coordinate
                second_worst_coord = simplex[len(simplex)-2]
                second_worst_perf_cost = perf_costs[len(perf_costs)-2]

                # calculate centroid
                centroid = self.__getCentroid(simplex[:len(simplex)-1])

                # reflection
                refl_coords = self.__getReflection(worst_coord, centroid)
                refl_perf_costs = map(self.getPerfCost, refl_coords)
                refl_perf_cost = min(refl_perf_costs)
                ipos = refl_perf_costs.index(refl_perf_cost)
                refl_coord = refl_coords[ipos]

                # the replacement of the worst coordinate
                next_coord = None
                next_perf_cost = None
            
                # if cost(best) <= cost(reflection) < cost(2nd_worst)
                if best_perf_cost <= refl_perf_cost < second_worst_perf_cost:
                    next_coord = refl_coord
                    next_perf_cost = refl_perf_cost
                    if self.verbose: print '--> reflection to %s' % next_coord 

                # if cost(reflection) < cost(best)
                elif refl_perf_cost < best_perf_cost:

                    # expansion
                    exp_coords = self.__getExpansion(refl_coord, centroid)
                    exp_perf_costs = map(self.getPerfCost, exp_coords)
                    exp_perf_cost = min(exp_perf_costs)
                    ipos = exp_perf_costs.index(exp_perf_cost)
                    exp_coord = exp_coords[ipos]

                    # if cost(expansion) < cost(reflection)
                    if exp_perf_cost < refl_perf_cost:
                        next_coord = exp_coord
                        next_perf_cost = exp_perf_cost
                        if self.verbose: print '--> expansion to %s' % next_coord 
                    else:
                        next_coord = refl_coord
                        next_perf_cost = refl_perf_cost
                        if self.verbose: print '--> reflection to %s' % next_coord 
                        
                # if cost(reflection) < cost(worst)
                elif refl_perf_cost < worst_perf_cost:

                    # outer contraction
                    cont_coords = self.__getContraction(refl_coord, centroid)
                    cont_perf_costs = map(self.getPerfCost, cont_coords)
                    cont_perf_cost = min(cont_perf_costs)
                    ipos = cont_perf_costs.index(cont_perf_cost)
                    cont_coord = cont_coords[ipos]
                    
                    # if cost(contraction) < cost(reflection)
                    if cont_perf_cost < refl_perf_cost:
                        next_coord = cont_coord
                        next_perf_cost = cont_perf_cost
                        if self.verbose: print '--> outer contraction to %s' % next_coord 

                # if cost(reflection) >= cost(worst)
                else:
                
                    # inner contraction
                    cont_coords = self.__getContraction(worst_coord, centroid)
                    cont_perf_costs = map(self.getPerfCost, cont_coords)
                    cont_perf_cost = min(cont_perf_costs)
                    ipos = cont_perf_costs.index(cont_perf_cost)
                    cont_coord = cont_coords[ipos]

                    # if cost(contraction) < cost(worst)
                    if cont_perf_cost < worst_perf_cost:
                        next_coord = cont_coord
                        next_perf_cost = cont_perf_cost
                        if self.verbose: print '--> inner contraction to %s' % next_coord 

                # if shrinkage is needed
                if next_coord == None and next_perf_cost == None:

                    # shrinkage
                    simplex = self.__getShrinkage(best_coord, simplex)
                    perf_costs = map(self.getPerfCost, simplex)
                    if self.verbose: print '--> shrinkage on %s' % best_coord 
                    
                # replace the worst coordinate with the better coordinate
                else:
                    simplex.pop()
                    perf_costs.pop()
                    simplex.append(next_coord)
                    perf_costs.append(next_perf_cost)
                
            # get the best simplex coordinate and its performance cost
            best_simplex_coord = simplex[0]
            best_simplex_perf_cost = perf_costs[0]
            old_best_simplex_perf_cost = best_simplex_perf_cost

            if self.verbose: print ('-> best simplex coordinate: %s, cost: %s' %
                                    (best_simplex_coord, best_simplex_perf_cost))
            
            # check if the time is not up yet
            if self.time_limit <= 0 or (time.time()-start_time) <= self.time_limit:

                # perform a local search on the best simplex coordinate
                (best_simplex_coord,
                 best_simplex_perf_cost) = self.searchBestNeighbor(best_simplex_coord,
                                                                   self.local_distance)
                
                # if the neighboring coordinate has a better performance cost
                if best_simplex_perf_cost < old_best_simplex_perf_cost:
                    if self.verbose: print ('---> better neighbor found: %s, cost: %s' %
                                            (best_simplex_coord, best_simplex_perf_cost))

            # compare to the global best coordinate and its performance cost
            if best_simplex_perf_cost < best_global_perf_cost:
                best_global_coord = best_simplex_coord
                best_global_perf_cost = best_simplex_perf_cost
                if self.verbose: print ('>>>> best coordinate found: %s, cost: %s' %
                                        (best_global_coord, best_global_perf_cost))

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
                                                                     
        if self.verbose: print '----- end simplex search -----'
        if self.verbose: print '----- begin summary -----'
        if self.verbose: print (' best coordinate: %s, cost: %s' %
                                (best_global_coord, best_global_perf_cost))
        if self.verbose: print ' total search time: %.2f seconds' % search_time
        if self.verbose: print ' total completed runs: %s' % runs
        if self.verbose: print '----- end summary -----'
 
        # return the best coordinate
        return best_global_coord

