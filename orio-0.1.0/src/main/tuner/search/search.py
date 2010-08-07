#
# The search engine used for search space exploration
#

import math, sys
from ...random import uniform, randint  # a workaround to resolve the 'random' package-name conflict

#-----------------------------------------------------

class Search:
    '''The search engine used to explore the search space '''

    MAXFLOAT = float('inf')

    #----------------------------------------------------------
    
    def __init__(self, cfrags, axis_names, axis_val_ranges, constraint, time_limit, total_runs,
                 search_opts, cmd_line_opts, ptcodegen, ptdriver, odriver):
        '''To instantiate a search engine'''

        # the class variables that are essential to know when developing a new search engine subclass
        self.time_limit = time_limit
        self.total_runs = total_runs
        self.search_opts = search_opts
        self.total_dims = len(axis_names)
        self.dim_uplimits = [len(r) for r in axis_val_ranges]
        self.space_size = 0
        if self.total_dims > 0:
            self.space_size = reduce(lambda x,y: x*y, self.dim_uplimits, 1)
        self.verbose = cmd_line_opts.verbose
        
        # the class variables that may be ignored when developing a new search engine subclass
        self.cfrags = cfrags
        self.axis_names = axis_names
        self.axis_val_ranges = axis_val_ranges
        self.constraint = constraint
        self.cmd_line_opts = cmd_line_opts
        self.ptcodegen = ptcodegen
        self.ptdriver = ptdriver
        self.odriver = odriver
        self.perf_cost_records = {}
        
    #----------------------------------------------------------

    def searchBestCoord(self):
        '''
        To explore the search space and retun the coordinate that yields the best performance
        (i.e. minimum performance cost).
        This is the function that needs to be created in each new search engine subclass.
        '''
        raise NotImplementedError('%s: unimplemented abstract function "searchBestCoord"' %
                                  self.__class__.__name__)
    
    #----------------------------------------------------------

    def search(self):
        '''To initiate the search process and return the best performance parameters'''

        # if the search space is empty
        if self.total_dims == 0:
            return {}

        # find the coordinate resulting in the best performance
        best_coord = self.searchBestCoord()

        # if no best coordinate can be found
        if best_coord == None:
            print ('error: the search cannot find a valid set of performance parameters. ' +
                   'the search time limit might be too short, or the performance parameter ' +
                   'constraints might prune out the entire search space.')
            sys.exit(1)

        # get the performance cost of the best parameters
        best_perf_cost = self.getPerfCost(best_coord)

        # convert the coordinate to the corresponding performance parameters
        best_perf_params = self.coordToPerfParams(best_coord)

        # return the best performance parameters
        return (best_perf_params, best_perf_cost)

    #----------------------------------------------------------

    def coordToPerfParams(self, coord):
        '''To convert the given coordinate to the corresponding performance parameters'''

        # get the performance parameters that correspond the given coordinate
        perf_params = {}
        for i in range(0, self.total_dims):
            ipoint = coord[i]
            perf_params[self.axis_names[i]] = self.axis_val_ranges[i][ipoint]

        # return the obtained performance parameters
        return perf_params

    #----------------------------------------------------------

    def getPerfCost(self, coord):
        '''
        To empirically evaluate the performance cost of the code variant generated using
        the given coordinate
        '''

        # if the given coordinate is out of the search space
        for i in range(0, self.total_dims):
            iuplimit = self.dim_uplimits[i]
            if coord[i] < 0 or coord[i] >= iuplimit:
                return self.MAXFLOAT

        # if the given coordinate has been computed before
        coord_key = str(coord)
        if coord_key in self.perf_cost_records:
            return self.perf_cost_records[coord_key]
        
        # get the performance parameters
        perf_params = self.coordToPerfParams(coord)
        
        # test if the performance parameters are valid
        try:
            is_valid = eval(self.constraint, perf_params)
        except Exception, e:
            print 'error: failed to evaluate the constraint expression: "%s"' % self.constraint
            print ' --> %s: %s' % (e.__class__.__name__, e)
            sys.exit(1)

        # if invalid performance parameters
        if not is_valid:
            perf_cost = self.MAXFLOAT
            
        # empirically evaluate the performance cost
        else:
            transformed_code_seq = self.odriver.optimizeCodeFrags(self.cfrags, perf_params)
            if len(transformed_code_seq) != 1:
                print 'internal error: the optimized annotation code cannot be multiple versions'
                sys.exit(1)
            transformed_code, _ = transformed_code_seq[0]
            test_code = self.ptcodegen.generate(transformed_code)
            perf_cost = self.ptdriver.run(test_code)

        # remember the performance cost of the given coordinate
        self.perf_cost_records[coord_key] = perf_cost
        
        # return the performance cost
        return perf_cost

    #----------------------------------------------------------

    def factorial(self, n):
        '''Return the factorial value of the given number'''
        return reduce(lambda x,y: x*y, range(1, n+1), 1)

    def roundInt(self, i):
        '''Proper rounding for integer'''
        return int(round(i))

    def getRandomInt(self, lbound, ubound):
        '''To generate a random integer N such that lbound <= N <= ubound'''
        if lbound > ubound:
            print ('internal error: the lower bound of genRandomInt must not be ' +
                   'greater than the upper bound')
            sys.exit()
        return randint(lbound, ubound)

    def getRandomReal(self, lbound, ubound):
        '''To generate a random real number N such that lbound <= N < ubound'''
        if lbound > ubound:
            print ('internal error: the lower bound of genRandomReal must not be ' +
                   'greater than the upper bound')
            sys.exit()
        return uniform(lbound, ubound)

    #----------------------------------------------------------

    def subCoords(self, coord1, coord2):
        '''coord1 - coord2'''
        return map(lambda x,y: x-y, coord1, coord2)
    
    def addCoords(self, coord1, coord2):
        '''coord1 + coord2'''
        return map(lambda x,y: x+y, coord1, coord2)

    def mulCoords(self, coef, coord):
        '''coef * coord'''
        return map(lambda x: self.roundInt(1.0*coef*x), coord)
    
    #----------------------------------------------------------

    def getCoordDistance(self, coord1, coord2):
        '''Return the distance between the given two coordinates'''

        d_sqr = 0
        for i in range(0, self.total_dims):
            d_sqr += (coord2[i] - coord1[i])**2
        d = math.sqrt(d_sqr)
        return d

    #----------------------------------------------------------

    def getRandomCoord(self):
        '''Randomly pick a coordinate within the search space'''

        random_coord = []
        for i in range(0, self.total_dims):
            iuplimit = self.dim_uplimits[i]
            ipoint = self.getRandomInt(0, iuplimit-1)
            random_coord.append(ipoint)
        return random_coord
                                                                     
    #----------------------------------------------------------

    def getNeighbors(self, coord, distance):
        '''Return all the neighboring coordinates (within the specified distance)'''
        
        # get all valid distances
        distances = [0] + range(1,distance+1,1) + range(-1,-distance-1,-1)

        # get all neighboring coordinates within the specified distance
        neigh_coords = [[]]
        for i in range(0, self.total_dims):
            iuplimit = self.dim_uplimits[i]
            cur_points = [coord[i]+d for d in distances]
            cur_points = filter(lambda x: 0 <= x < iuplimit, cur_points)
            n_neigh_coords = []
            for ncoord in neigh_coords:
                n_neigh_coords.extend([ncoord[:]+[p] for p in cur_points])
            neigh_coords = n_neigh_coords

        # remove the current coordinate from the neighboring coordinates list
        neigh_coords.remove(coord)
        
        # return all valid neighboring coordinates
        return neigh_coords

    #----------------------------------------------------------

    def searchBestNeighbor(self, coord, distance):
        '''
        Perform a local search by starting from the given coordinate then examining
        the performance costs of the neighboring coordinates (within the specified distance),
        then we perform this local search recursively once the neighbor with the best performance
        cost is found.
        '''

        # get all neighboring coordinates within the specified distance
        neigh_coords = self.getNeighbors(coord, distance)

        # record the best neighboring coordinate and its performance cost so far
        best_coord = coord
        best_perf_cost = self.getPerfCost(coord)

        # examine all neighboring coordinates
        for n in neigh_coords:
            perf_cost = self.getPerfCost(n)
            if perf_cost < best_perf_cost:
                best_coord = n
                best_perf_cost = perf_cost

        # recursively apply this local search, if new best neighboring coordinate is found
        if best_coord != coord:
            return self.searchBestNeighbor(best_coord, distance)
        
        # return the best neighboring coordinate and its performance cost
        return (best_coord, best_perf_cost)
    

