#!/usr/bin/env python

import re,os,shutil

basefile = 'axpy4'
sizes = [10, 100,1000,10000,50000,100000,500000,1000000,2000000,5000000]

for s in sizes:
	f = open('axpy.bg.spec','r')
	contents = f.read()
	f.close()
	contents = re.sub('@THESIZE@',str(s),contents)
	fname = 'axpy' + str(s) + '.bg.spec'
	f = open(fname,'w')
	f.write(contents)
	f.close()
        srcfile = basefile + '_' + str(s) + '.c'	
	shutil.copyfile(basefile + '.c',srcfile)
	# run the tests
	cmd = 'ancc -v -s ' + fname + ' ' + srcfile + ' > axpy' + str(s) + '.bg.output.txt 2>&1 &'	
	print cmd
	os.system(cmd)
