#!/usr/bin/env python

# A Python script used for installing the Orio tool

#-----------------------------------------------------------

import os, sys
from distutils.core import setup

#-----------------------------------------------------------

# to traverse the source code directory to get all python packages
py_packages = []
cur_dir = os.getcwd()
src_dir = os.path.join(cur_dir, 'src')
for root, dirs, files in os.walk(src_dir, topdown=True):
    if '__init__.py' in files:
        rel_dir = root[len(cur_dir)+1:]
        dir_names = rel_dir.split(os.sep)
        py_packages.append('.'.join(['orio'] + dir_names[1:]))

#-----------------------------------------------------------

# make a call to the setup function
setup(name = 'orio',
      version = '0.0.1',
      description = 'ORIO -- An Annotation-Based Performance Tuning Tool',
      author = 'Albert Hartono',
      author_email = 'hartonoa@cse.ohio-state.edu',
      maintainer = 'Albert Hartono',
      maintainer_email = 'hartonoa@cse.ohio-state.edu',
      url = 'https://trac.mcs.anl.gov/projects/performance/wiki/Orio',
      packages = py_packages,
      package_dir = {'orio' : 'src'},
      package_data = {'orio' : ['tool/zestyparser/*']},
      scripts = ['orcc', 'orf'])


