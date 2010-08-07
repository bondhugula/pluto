#
# To compile and execute the performance-testing code to get the performance cost
#

import os, sys, re

#-----------------------------------------------------

class PerfTestDriver:
    '''
    The performance-testing driver used to compile and execute the testing code
    to get the performance cost
    '''

    # the file names of the testing code (i.e. source and executable)
    __PTEST_SRC_FNAME = '_orio_perftest.c'

    # types of performance-counting methods
    __PCOUNT_BASIC = 'basic timer'
    
    #-----------------------------------------------------
    
    def __init__(self, build_cmd, build_opts, pcount_method, pcount_reps, driversrc=None, run_cmd=None, verbose=False):
        '''To instantiate the performance-testing driver'''

        self.build_cmd = build_cmd
        self.build_opts = build_opts
        self.pcount_method = pcount_method        # TODO
        self.pcount_reps = pcount_reps

        # TODO: add tuning spec options for the following two values
        self.batch = False
        if run_cmd: self.batch = True
        if driversrc:

            self.srcname = driversrc                  # if specified, use that file as the driver
            self.gendriver = False
        else: 
            self.srcname = self.__PTEST_SRC_FNAME
            if self.batch: 
                import random
                self.srcname = '_orio_perftest' + str(random.randint(1, 14141)) + '.c'
            self.gendriver = True
        self.execname = self.srcname[:self.srcname.find('.')] + '.exe'
        self.run_cmd = run_cmd
        self.verbose = verbose

        if self.pcount_method != self.__PCOUNT_BASIC:
            print 'error: unknown performance-counting method: "%s"' % self.pcount_method
            sys.exit(1)

    #-----------------------------------------------------

    def __write(self, test_code):
        '''Write the testing code into a file'''

        try:
            f = open(self.srcname, 'w')
            f.write(test_code)
            f.close()
        except:
            print 'error: cannot open file for writing: %s' % self.srcname
            sys.exit(1)

    #-----------------------------------------------------

    def __build(self):
        '''Compile the testing code'''
        #import distutils.sysconfig
        
        if not self.build_cmd.startswith('make'):
            # add extra options
            extra_compiler_opts = ''
            extra_compiler_opts += '-DREPS=%s' % self.pcount_reps
            
            # compile the testing code
            cmd = ('%s %s %s -o %s %s' %
                   (self.build_cmd, self.build_opts, extra_compiler_opts,
                    self.execname, self.srcname))
            if self.verbose: print ' compiling:\n\t' + cmd
            status = os.system(cmd)
            if status:
                print 'error: failed to compile the testing code: "%s"' % cmd
                sys.exit(1)
        else:
            from subprocess import Popen, PIPE
            
            # Using make
            # add extra defines
            extra_vars = ' AUTOTUNE_REPS=%s' % self.pcount_reps
            
            cmd = ('%s %s %s > build.log 2>&1' %
                   self.build_cmd, self.build_opts, extra_vars)
            if self.verbose: print ' compiling:\n\t' + cmd
            try:
                p1 = Popen([cmd], stdout=PIPE, stderr=PIPE)
            except Exception,e:
                print 'build command error: %s' % str(e)
                sys.exit(1)
                
            try:
                p2 = Popen(["tee", "build.log"], stdin=p1.stdout, stdout=PIPE)
            except Exception,e:
                print 'warning: could not save build output in build.log file: %s' % str(e)

            # save build output
            if p2:
                output = p2.communicate()[0]
                print output
            
            

    #-----------------------------------------------------

    def __execute(self):
        '''Execute the executable to get the performance cost'''

        # execute the executable and collect the performance cost
        perf_costs = {}
        if os.path.exists(self.execname):
            if self.run_cmd == None:
                cmd = './%s' % self.execname
                if self.verbose: print ' running:\n\t'  + cmd
                try:
                    f = os.popen(cmd)    
                    # The output is the list of times using python list syntax, e.g., {'[0,1]':0.2, '[1,1]': 0.3]
                    output = f.read()
                    if output: perf_costs = eval(output)
                    f.close()
                except Exception, e:
                    print 'error: failed to execute the testing code: "%s"' % cmd
                    print ' --> %s: %s' % (e.__class__.__name__, e)
                    sys.exit(1)
            else:
                cmd = self.run_cmd + ' ' + self.execname
                if self.verbose: print ' running:\n\t'  + cmd
                # TODO: redo this to take output file name
                try:
                    f = os.popen(cmd)    
                    # The output is the list of times using python list syntax, e.g., {'[0,1]':0.2, '[1,1]': 0.3]
                    output = f.read()
                    if self.batch and output:
                        import time
                        time.sleep(3);
                        status = '1'
                        while status == '1': 
                            # FIX TODO: very bad assumption that the last number out is the batch job name
                            jobid = output.strip().split('\n')[-1]
                            # TODO: make this an option in tuning spec (status command):
                            chkstatus_command = 'cqstat %s | grep %s | wc -l' % (jobid, jobid)
                            f2 = os.popen(chkstatus_command)
                            status = f2.read().strip()
                            f2.close()
                        outfile = jobid + '.output'
                        while not os.path.exists(outfile):
                            time.sleep(5)
                        try:
                            f3 = open(outfile)
                            output = f3.read()
                            if output: perf_costs = eval(output)
                            f3.close()
                        except:
                            perf_costs = {}
                    else:    
                        if output: perf_costs = eval(output)
                        f.close()
                except Exception, e:
                    print 'error: failed to execute the testing code: "%s"' % cmd
                    print ' --> %s: %s' % (e.__class__.__name__, e)
                    sys.exit(1)


        
        # check if the performance cost is already obtained
        if not perf_costs:
            print 'error: performance testing failed: "%s"' % cmd
            sys.exit(1)

        # return the performance costs dictionary (indexed by the string representation of the search coordinates)
        return perf_costs
            
    #-----------------------------------------------------

    def __cleanup(self):
        '''Delete all the generated files'''

        flist = [self.execname]
        if self.gendriver: flist.append(self.srcname)
        for fname in flist:
            try:
                if os.path.exists(fname):
                    os.unlink(fname)
            except:
                print 'error: cannot delete file: %s' % fname
                sys.exit(1)

    #-----------------------------------------------------

    def run(self, test_code):
        '''To compile and to execute the given testing code to get the performance cost
        @param test_code: the code for testing multiple coordinates in the search space
        @return: a dictionary of the times corresponding to each coordinate in the search space
        '''

        if self.gendriver: self.__write(test_code)
        self.__build()
        # Return a dictionary of performance codes indexed by the string represeentation of search coordinates
        perf_costs = self.__execute()
        self.__cleanup()
        return perf_costs


